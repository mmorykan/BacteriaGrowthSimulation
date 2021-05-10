/**
 * Runs the organism-evolution simulator.
 * 
 * This version runs in serial. Compile with:
 *     gcc -Wall -O3 -march=native simulation-serial.c organism.c -o simulation-serial
 */
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <unistd.h>

#include "organism.h"

/**
 * Spawn an organism in the environment with the given genome (already
 * converted) at the given location. The location is guaranteed to be within
 * the environment (from 0,0 to WORLD_SIZE,WORLD_SIZE). The genome must be
 * immediately used or copied as after calling this function that memory may
 * be overwritten.
 * 
 * This function is *NOT* implemented in organism.c. You must implement this
 * method yourself in simulation.c. The serial version is provided there
 * and just calls organism_set(). In the MPI version you may sometimes have
 * to call an MPI function.
 */
void organism_spawn(environment* env, uint8_t* genome, uint32_t size, uint32_t i, uint32_t j) {
    organism_set(&env->grid[i*env->world_size+j], genome, size);
}

/**
 * Records that an organism has completed a task. The task value is one of the
 * constants TASK_* such as TASK_OR and TASK_EQUAL.
 * 
 * This function is *NOT* implemented in organism.c. You must implement this
 * method yourself in simulation.c. The serial version is provided and updates
 * array as appropriate along with granting the organism additional energy.
 */
void organism_completed_task(organism* o, int task) {
    o->energy += ENERGY_PER_TASK[task];
    o->tasks_completed[task]++;
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

    // start the timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // seed random number generator
    srand(random_seed);

    // run the simulation
    for (size_t t = 0; t < iterations; t++) {
        for (size_t i = 0; i < world_size; i++) {
            for (size_t j = 0; j < world_size; j++) {
                organism_update(&env.grid[i*world_size + j]);
            }
        }
        //if (t % (1024*1024) == 0) { environment_print(&env); } // for debugging
    }
    //if (iterations % (1024*1024) != 0) { environment_print(&env); } // for debugging

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
}
