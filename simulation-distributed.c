/**
 * Runs the organism-evolution simulator.
 * 
 * This version runs in parallel using distributed and shared memory. Compile with:
 *     mpicc -Wall -O3 -march=native -fopenmp simulation-distributed.c organism.c util.c -o simulation-distributed
 * 
 * Authors: Jonah Beers and Mark Morykan
 */
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>

#include "util.h"
#include "organism.h"

#include <stdint.h>
#include <limits.h>
#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "MPI_SIZE_T indeterminant"
#endif

/**
 * Spawn an organism in the environment with the given genome (already
 * converted) at the given location. The location is guaranteed to be within
 * the environment (from 0,0 to WORLD_SIZE,WORLD_SIZE). The genome must be
 * immediately used or copied as after calling this function that memory may
 * be overwritten.
 */
void organism_spawn(environment* env, uint8_t* genome, uint32_t size, uint32_t i, uint32_t j) {
    // Save the organisms to spawn at the end of the current iteration
    #pragma omp critical
    {
        size_t num_organisms = env->num_stored_organisms;
        env->organisms[num_organisms] = &env->grid[i*env->world_size+j];
        env->genomes[num_organisms] = genome;
        env->sizes[num_organisms] = size;
        env->num_stored_organisms++;
    }
}

/**
 * Records that an organism has completed a task. The task value is one of the
 * constants TASK_* such as TASK_OR and TASK_EQUAL.
 */
void organism_completed_task(organism* o, int task) {
    o->energy += ENERGY_PER_TASK[task];
    o->tasks_completed[task]++;
    #pragma omp atomic
    o->env->tasks_completed[task]++;
}


/**
 * Receives organism info from another process and sets them in the environment.
 */
void organism_recv_and_set(int from_rank, MPI_Datatype org_info, organism* env_grid, size_t world_size) {
    size_t amt_to_set;
    MPI_Recv(&amt_to_set, 1, MPI_SIZE_T, from_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (amt_to_set > 0) {
        work_to_send from_rank_orgs[amt_to_set];
        MPI_Recv(from_rank_orgs, amt_to_set, org_info, from_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (size_t k = 0; k < amt_to_set; k++) {
            work_to_send org = from_rank_orgs[k];
            organism_set(&env_grid[org.location_i * world_size + org.location_j], org.genome, org.size);
        }
    }
}


/**
 * Sends the number of organisms to send to another process and then the actual organisms to another process.
 */
void organism_send(int to_rank, MPI_Datatype org_info, work_to_send* orgs, size_t num_orgs) {
    MPI_Send(&num_orgs, 1, MPI_SIZE_T, to_rank, 0, MPI_COMM_WORLD);
    if (num_orgs > 0) MPI_Send(orgs, num_orgs, org_info, to_rank, 0, MPI_COMM_WORLD);
}


/**
 * Creates a struct containing organism info to send between processes.
 */
work_to_send create_org_to_send(uint8_t* genome, uint32_t size, size_t location_i, size_t location_j) {
    work_to_send work;
    memcpy(work.genome, genome, sizeof(work.genome));
    work.size = size, work.location_i = location_i, work.location_j = location_j;
    return work;
}


/**
 * Function for a worker processes. Receives a section of the grid to work
 * on from process 0, updates and sets organisms, communicates with neighboring 
 * processes if organisms spawn in a different section of the environment, and sends 
 * all organisms from its side back to process 0 after all iterations have finished.
 */
void organism_worker(MPI_Datatype org_info, MPI_Datatype final_org, double mutation_prob, size_t iterations, size_t world_size, 
                     size_t grid_size, const char* input_path) {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); // This process' rank
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size); // Number of processes
    int lower_rank = rank - 1, upper_rank = rank + 1; // Rank above/below this process
    size_t num_threads = get_num_cores_affinity(); 

    // Setup environment
    environment env; 
    env.world_size = world_size, env.mutation_prob = mutation_prob;
    if (!environment_load_from_path(&env, input_path)) { return; }
    environment_clear_tasks(&env);
    
    // Receive the number of organisms and index of grid to start at from process 0
    size_t process_info[2];
    MPI_Recv(&process_info, 2, MPI_SIZE_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    size_t num_organisms = process_info[0], beginning = process_info[1], ending = beginning + num_organisms;

    // Allocate memory for environment organisms, genomes, and sizes
    env.organisms = (organism**) malloc(num_organisms*sizeof(organism*));
    env.genomes = (u_int8_t**) malloc(num_organisms*sizeof(u_int8_t*));
    env.sizes = (uint32_t*) malloc(num_organisms*sizeof(uint32_t));
    env.num_stored_organisms = 0; // start at 0 organsims stored

    for (size_t t = 0; t < iterations; t++) {
        // Update this rank's section of the grid
        #pragma omp parallel for num_threads(num_threads) schedule(static, 32)
        for (size_t i = beginning; i < ending; i++) {
            organism_update(&env.grid[i]);
        }

        // Counters and arrays for data to send to other processes
        size_t upper_rank_send_count = 0, lower_rank_send_count = 0;
        work_to_send lower_rank_orgs[world_size], upper_rank_orgs[world_size];

        // Set all organisms that tried to spawn and package up organisms that tried to spawn outside of this processes's environment
        for (size_t j = 0; j < env.num_stored_organisms; j++) {
            organism* org = env.organisms[j];
            size_t location_i = org->location_i, location_j = org->location_j;
            size_t location = location_i * world_size + location_j;
            if (location < beginning || location >= ending) { // Organism spawned outside of rank's grid section
                work_to_send work = create_org_to_send(env.genomes[j], env.sizes[j], location_i, location_j);
                if (location < beginning) lower_rank_orgs[lower_rank_send_count++] = work; // Send to rank below
                else upper_rank_orgs[upper_rank_send_count++] = work; // Send to rank above
            } else { // Organism spawed in this rank's grid section
                organism_set(env.organisms[j], env.genomes[j], env.sizes[j]);
            }
        }
        env.num_stored_organisms = 0; // reset organisms stored
        
        // Sends/receives organism info to/from rank below
        organism_send(lower_rank, org_info, lower_rank_orgs, lower_rank_send_count); 
        organism_recv_and_set(lower_rank, org_info, env.grid, world_size);
        if (rank < size - 1) { // Sends/receives organism info to/from rank above (if exists)
            organism_send(upper_rank, org_info, upper_rank_orgs, upper_rank_send_count);
            organism_recv_and_set(upper_rank, org_info, env.grid, world_size);
        }
    }

    // Send back the environment organisms and tasks completed
    MPI_Send(env.grid + beginning, num_organisms, final_org, 0, 0, MPI_COMM_WORLD);
    MPI_Send(env.tasks_completed, NUM_TASKS, MPI_SIZE_T, 0, 0, MPI_COMM_WORLD);
}


int main(int argc, char*const argv[]) {
    // default command-line options
    size_t iterations = 4*1024*1024;
    unsigned int random_seed = time(NULL);
    double mutation_prob = 0.02;
    size_t num_threads = get_num_cores_affinity();

    // parse command line arguments
    int opt;
    while ((opt = getopt(argc, argv, "n:r:m:")) != -1) {
        char* end;
        switch (opt) {
        case 'n': iterations = strtoumax(optarg, &end, 10); break;
        case 'r': random_seed = strtoul(optarg, &end, 10); break;
        case 'm': mutation_prob = atof(optarg); break;
        default:
            fprintf(stderr, "usage: %s [-n num-iterations] [-r random-seed] [-m mutation-prob] input output\n", argv[0]);
            return 1;
        }
    }
    if (optind + 2 < argc || iterations == 0) {
        fprintf(stderr, "usage: %s [-n num-iterations] [-r random-seed] [-m mutation-prob] input output\n", argv[0]);
        return 1;
    }
    const char* input_path = argv[optind];
    const char* output_path = argv[optind+1];

    // load the environment
    environment env;
    if (!environment_load_from_path(&env, input_path)) { return 1; }
    environment_clear_tasks(&env);
    env.mutation_prob = mutation_prob;
    size_t world_size = env.world_size;
    size_t grid_size = world_size * world_size;

    // initialize MPI
    MPI_Init(NULL, NULL);
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);

    // The information needed to spawn an organism
    MPI_Datatype org_info;
    MPI_Type_create_struct(4, org_info_field_counts, org_info_field_offsets, org_info_field_types, &org_info);
    MPI_Type_commit(&org_info);

    // The information needed to set other processes' organisms into the master processe's environment
    MPI_Datatype final_org;
    MPI_Type_create_struct(3, final_org_field_counts, final_org_field_offsets, final_org_field_types, &final_org);
    MPI_Type_commit(&final_org);

    double rows_per_worker = world_size / (double) size;

    // start the timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);    

    // seed random number generator
    srand(random_seed+rank);

    // run the simulation
    if (rank != 0) {
        organism_worker(org_info, final_org, mutation_prob, iterations, world_size, grid_size, input_path);
        MPI_Finalize();
        return 0;
    }

    size_t beginning[size-1];  // Each processe's environment starting index 
    size_t num_organisms[size-1];  // The number of oranisms each process is updating
    size_t process_info[2];  // The starting index and number of organisms for each process
    for (size_t i = 1; i < size; i++) {
        size_t row_start = (size_t) (i * rows_per_worker);
        size_t amt_of_rows = (size_t) (((i + 1) * rows_per_worker) - row_start);
        num_organisms[i-1] = process_info[0] = amt_of_rows * world_size; // number of organisms
        beginning[i-1] = process_info[1] = row_start * world_size;
        MPI_Send(process_info, 2, MPI_SIZE_T, i, 0, MPI_COMM_WORLD);
    }

    // Set up fields to store organisms before setting them at end of iteration
    size_t num_of_organisms = (size_t) (rows_per_worker * world_size);
    env.organisms = (organism**) malloc(num_of_organisms*sizeof(organism*));
    env.genomes = (u_int8_t**) malloc(num_of_organisms*sizeof(u_int8_t*));
    env.sizes = (uint32_t*) malloc(num_of_organisms*sizeof(uint32_t));
    env.num_stored_organisms = 0;

    for (size_t t = 0; t < iterations; t++) {

        // Update process 0's section of the grid
        #pragma omp parallel for num_threads(num_threads) schedule(static, 32)
        for (size_t i = 0; i < num_of_organisms; i++) {
            organism_update(&env.grid[i]);
        }
        
        // Set all organisms that tried to spawn and package up organisms that tried to spawn outside of this processe's environment
        size_t upper_rank_send_count = 0;
        work_to_send upper_rank_orgs[world_size];
        for (size_t j = 0; j < env.num_stored_organisms; j++) {
            organism* org = env.organisms[j];
            size_t location_i = org->location_i, location_j = org->location_j;
            size_t location = location_i * world_size + location_j;
            if (location >= beginning[0]) {
                // Create organism to send to upper process
                upper_rank_orgs[upper_rank_send_count++] = create_org_to_send(env.genomes[j], env.sizes[j], location_i, location_j);
            } else {
                organism_set(env.organisms[j], env.genomes[j], env.sizes[j]);
            }
        }
        env.num_stored_organisms = 0;

        if (size > 1) {
            // Send organisms to upper process and receive and set upper processe's organisms that spawned in this environment
            organism_send(1, org_info, upper_rank_orgs, upper_rank_send_count);
            organism_recv_and_set(1, org_info, env.grid, world_size);
        }
    }

    // Receive other workers environments and numbers of tasks completed
    size_t tasks_completed[NUM_TASKS];
    for (size_t k = 1; k < size; k++) {
        MPI_Recv(env.grid + beginning[k-1], num_organisms[k-1], final_org, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(tasks_completed, NUM_TASKS, MPI_SIZE_T, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (size_t m = 0; m < NUM_TASKS; m++) {
            env.tasks_completed[m] += tasks_completed[m];
        }
    }
    
    // get the elapsed time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = end.tv_sec-start.tv_sec+(end.tv_nsec-start.tv_nsec)/1000000000.0;
    printf("Time: %g secs\n", time);

    // save the results
    environment_save_to_path(&env, output_path);

    // print the tasks completed
    printf("NOT:     %zu\n", env.tasks_completed[TASK_NOT]);
    printf("NAND:    %zu\n", env.tasks_completed[TASK_NAND]);
    printf("AND:     %zu\n", env.tasks_completed[TASK_AND]);
    printf("OR-NOT:  %zu\n", env.tasks_completed[TASK_OR_NOT]);
    printf("OR:      %zu\n", env.tasks_completed[TASK_OR]);
    printf("AND-NOT: %zu\n", env.tasks_completed[TASK_AND_NOT]);
    printf("NOR:     %zu\n", env.tasks_completed[TASK_NOR]);
    printf("XOR:     %zu\n", env.tasks_completed[TASK_XOR]);
    printf("EQUAL:   %zu\n", env.tasks_completed[TASK_EQUAL]);

    // cleanup
    free(env.grid); free(env.organisms); free(env.genomes); free(env.sizes);
    MPI_Type_free(&final_org); MPI_Type_free(&org_info);
    MPI_Finalize();
}
