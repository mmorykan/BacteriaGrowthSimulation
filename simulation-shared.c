/**
 * Runs the organism-evolution simulator.
 * 
 * This version runs in parallel using shared memory. Compile with:
 *     gcc -Wall -O3 -march=native -fopenmp simulation-shared.c organism.c util.c -o simulation-shared
 * 
 * Authors: Jonah Beers and Mark Morykan
 */
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <unistd.h>
#include <omp.h>

#include "util.h"
#include "organism.h"

/**
 * Spawn an organism in the environment with the given genome (already
 * converted) at the given location. The location is guaranteed to be within
 * the environment (from 0,0 to WORLD_SIZE,WORLD_SIZE). The genome must be
 * immediately used or copied as after calling this function that memory may
 * be overwritten.
 */
void organism_spawn(environment* env, uint8_t* genome, uint32_t size, uint32_t i, uint32_t j) {
    // Save the organism to set it at the end of the current iteration
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


int main(int argc, char*const argv[]) {
    // default command-line options
    size_t iterations = 4*1024*1024;
    unsigned int random_seed = time(NULL);
    double mutation_prob = 0.02;

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

    // Allocate space for organisms to spawn at the end of each time step
    size_t grid_size = world_size * world_size;
    env.organisms = malloc(grid_size*sizeof(organism*));
    env.genomes = malloc(grid_size*sizeof(u_int8_t*));
    env.sizes = malloc(grid_size*sizeof(uint32_t));
    env.num_stored_organisms = 0;

    // get number of threads
    size_t num_threads = get_num_cores_affinity();

    // start the timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // seed random number generator so all processes generate the same sequence of numbers 
    srand(random_seed);

    // run the simulation 
    for (size_t t = 0; t < iterations; t++) {

        // Update all organisms
        #pragma omp parallel for num_threads(num_threads) schedule(static, 32)
        for (size_t i = 0; i < grid_size; i++) {
            organism_update(&env.grid[i]);
        }

        // Set all organisms that tried to spawn at end of iteration
        for (size_t j = 0; j < env.num_stored_organisms; j++) {
            organism_set(env.organisms[j], env.genomes[j], env.sizes[j]);
        }
        env.num_stored_organisms = 0;
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
    free(env.grid);
    free(env.organisms);
    free(env.genomes);
    free(env.sizes);
}
