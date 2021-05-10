/**
 * Defines the interface for an organism and the environment. Two of the
 * functions are not defined in organism.c and must be implemented elsewhere.
 * This file should mostly not need to be modifed with the exception being the
 * environment struct. 
 */
#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>


////////////////////////////
///// Global Constants /////
////////////////////////////

// Organism Properties and Constants
#define MAX_GENOME_LENGTH 256 // no organism will be created with a genome larger than this (official implementation has 2048, but our program is more restricted)
#define MAX_STACK_SIZE 8   // maximum size of the stacks (official implementation has 10)
#define NUM_STACKS 2       // number of stacks
#define MAX_TEMPLATE_LEN 8 // maximum length of templates (official implementation has 10), also applies to the amount of numebr of recent copied instructions should be remembered

// Definitions of the registers and heads in the organism
#define AX 0
#define BX 1
#define CX 2
#define NUM_REGISTERS 3
#define IP 0
#define READ_HEAD 1
#define WRITE_HEAD 2
#define FLOW_HEAD 3
#define NUM_HEADS 4

// Tasks to be completed
#define TASK_NOT     0
#define TASK_NAND    1
#define TASK_AND     2
#define TASK_OR_NOT  3
#define TASK_OR      4
#define TASK_AND_NOT 5
#define TASK_NOR     6
#define TASK_XOR     7
#define TASK_EQUAL   8
#define NUM_TASKS    9

// Energy they granted by tasks
static const uint32_t ENERGY_PER_TASK[] = {1024, 1024, 2048, 2048, 4096, 4096, 8192, 8192, 16384}; // last one is enough to fully replicate plus some


//////////////////////
///// Structures /////
//////////////////////

typedef struct _environment environment;


/**
 * An individual organism within the environment.
 * Each organism requires ~188 bytes more than just the double-length genome (total 512 with the MAX_GENOME_LENGTH of 256)
 */
typedef struct _organism {
    // NOTE: do NOT change this structure

    // information about the organism directly related to the CPU processing
    uint8_t memory[2*MAX_GENOME_LENGTH]; // cyclic memory buffer (CPU instructions and flags)
    uint32_t size; // used size of the memory buffer
    union {
        uint32_t heads[NUM_HEADS]; // "heads": instruction pointer, read, write, and flow heads
        struct { uint32_t ip, read_head, write_head, flow_head; }; // give names to the values in the array (sometimes array is easier to use)
    };
    union {
        int32_t registers[NUM_REGISTERS]; // registers AX, BX, and CX
        struct { int32_t ax, bx, cx; }; // give names to the values in the array (sometimes array is easier to use)
    };
    int32_t stacks[NUM_STACKS][MAX_STACK_SIZE]; // 2 stacks
    uint32_t stack_sizes[NUM_STACKS]; // the number of values currently in the stacks
    uint8_t stack_active; // which stack is active
    //void *in_buffer, *out_buffer; // these are not explicit within this implementation

    // additional information required for some instructions
    bool has_allocated; // if h-alloc has been called and h-divide hasn't been
    uint32_t energy; // the amount of energy the organism has, it takes 64 units to execute one instruction, gain 1 per update
    uint32_t genome_size; // original size of the genome
    uint32_t instructions_executed; // number of instructions executed (it's 'age')
    uint8_t recent_copies[MAX_TEMPLATE_LEN], recent_copies_offset, recent_copies_size; // recently copied values (cyclic buffer)
    uint32_t num_inputs; // number of inputs this organism has obtained
    uint32_t tasks_completed[NUM_TASKS]; // number of times this organism has completed each task

    // link to the environment
    environment* env; uint32_t location_i, location_j; // the environment this organism is in and its location within that environment
} organism;

/**
 * The environment containing all of the organisms.
 * Any organism with a size of 0 is non-existent.
 */
struct _environment {
    // NOTE: you may change this structure to save additional information as needed
    organism* grid; // 2D array of organisms; world_size-by-world_size
    size_t world_size; // the size of the grid
    double mutation_prob; // probability of mutation
    size_t tasks_completed[NUM_TASKS]; // number of times any organism has completed each task
    // Save organism locations, genomes and sizes
    size_t num_stored_organisms;
    organism** organisms;
    uint32_t* sizes;
    uint8_t** genomes;
};

/**
 * Organism information that needs to be sent when an 
 * organism spawns in another process' section of the grid.
 */
typedef struct _work_to_send {
    uint8_t genome[MAX_GENOME_LENGTH];
    uint32_t size, location_i, location_j;
} work_to_send;


/////////////////////////////////////////////////////
//// Organism Info Being Sent Between Processes /////
////////// for Creating Custom MPI Types ////////////
/////////////////////////////////////////////////////
#ifdef MPI_VERSION
static const int org_info_field_counts[] = {
    MAX_GENOME_LENGTH, 1, 1, 1
};
static const MPI_Aint org_info_field_offsets[] = {
    offsetof(work_to_send, genome), offsetof(work_to_send, size), offsetof(work_to_send, location_i), offsetof(work_to_send, location_j)
};
MPI_Datatype org_info_field_types[] = {
    MPI_BYTE, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED
};
#endif


#ifdef MPI_VERSION
static const int final_org_field_counts[] = {
    2*MAX_GENOME_LENGTH, 1, 1
};
static const MPI_Aint final_org_field_offsets[] = {
    offsetof(organism, memory), offsetof(organism, genome_size), offsetof(organism, size)
};
MPI_Datatype final_org_field_types[] = {
    MPI_BYTE, MPI_UNSIGNED, MPI_UNSIGNED
};
#endif


///////////////////////////////////////////////////////
///// Organism Info for Creating Custom MPI Types /////
///////////////////////////////////////////////////////
// To enable this section, you have to include mpi.h before you include this file
#ifdef MPI_VERSION
static const int organism_field_counts[] = {
    MAX_GENOME_LENGTH*2, 1, NUM_HEADS, NUM_REGISTERS, NUM_STACKS*MAX_STACK_SIZE, NUM_STACKS, 1,
    1, 1, 1, 1, MAX_TEMPLATE_LEN, 1, 1, 1, NUM_TASKS, 1, 1
};
static const MPI_Aint organism_field_offsets[] = {
    offsetof(organism, memory), offsetof(organism, size), offsetof(organism, heads), offsetof(organism, registers),
    offsetof(organism, stacks), offsetof(organism, stack_sizes), offsetof(organism, stack_active),
    offsetof(organism, has_allocated), offsetof(organism, energy), offsetof(organism, genome_size), offsetof(organism, instructions_executed),
    offsetof(organism, recent_copies), offsetof(organism, recent_copies_offset), offsetof(organism, recent_copies_size),
    offsetof(organism, num_inputs), offsetof(organism, tasks_completed), offsetof(organism, location_i), offsetof(organism, location_j)
};
static const MPI_Datatype organism_field_types[] = {
    MPI_BYTE, MPI_UNSIGNED, MPI_UNSIGNED, MPI_INT, MPI_INT, MPI_UNSIGNED, MPI_BYTE, MPI_BYTE, MPI_UNSIGNED, 
    MPI_UNSIGNED, MPI_UNSIGNED, MPI_BYTE, MPI_BYTE, MPI_BYTE, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED
};
#endif


///////////////////////////////////////////////
///// Core Environment/Organism Functions /////
///////////////////////////////////////////////

/**
 * Load the environment from the path to a given text file. See
 * environment_load_from_file() for details.
 */
bool environment_load_from_path(environment* env, const char* path);

/**
 * Load the environment from the given text file. The file must have a single
 * number on the first line for the world size and every line after that has
 * 2 numbers separated by whitespace for the location of an organism followed
 * by the genome for the organism (which is a sequence of lowercase letters).
 * This returns false if there are any issues reading the file or the data is
 * not a proper environment (it also prints out an error message).
 */
bool environment_load_from_file(environment* env, FILE* file);

/**
 * Save the environment to the path of the given text file. See
 * environment_save_to_file() for details.
 */
bool environment_save_to_path(environment* env, const char* path);

/**
 * Save the environment to the given text file. It will have the same format as
 * used by environment_load_from_file().
 */
bool environment_save_to_file(environment* env, FILE* file);

/**
 * Setup the grid for an environment. This allocates the memory, clears the
 * grid making all cells unoccupied / all organisms non-existent, and sets up a
 * few values so the organisms link back to the environment. The env->grid
 * field must be freed after this is called. Before calling this function, make
 * sure the world_size if set properly in the environment.
 */
void environment_alloc_grid(environment* env);

/**
 * Print out the environment. This displays it as a grid with spaces for empty
 * cells, single dots for occupied cells with weak (or the only) organisms, and
 * then blocks of various brightness for cells with stronger organisms.
 */
void environment_print(const environment* env);

/**
 * Set task counts for the environemnt to 0. Does not effect organisms.
 */
void environment_clear_tasks(environment* env);

/**
 * Check if a cell in the environment is occupied. Used like:
 *      cell_occupied(&env[i][j]);
 */
bool cell_occupied(organism* org);

/**
 * Converts a genome from an alphabetical sequence (e.g. "wzcagcccccccccc...")
 * to one usable by an organism.
 */
void convert_genome(const char* genome_alpha, uint8_t* genome, uint32_t size);

/**
 * Converts a genome from an alphabetical sequence (e.g. "wzcagcccccccccc...")
 * to one usable by an organism. The memory is copied and the copy is returned.
 * The returned memory must be freed.
 */
uint8_t* convert_genome_copy(const char* genome_alpha, uint32_t size);

/**
 * Converts a genome to an alphabetical sequence (e.g. "wzcagcccccccccc...").
 * The given buffer for genome_alpha must have room for the null terminator.
 */
void convert_genome_to_alpha(const uint8_t* genome, char* genome_alpha, uint32_t size);

/**
 * Converts a genome to an alphabetical sequence (e.g. "wzcagcccccccccc...").
 * to one usable by an organism. The memory is copied and the copy is returned.
 * The returned memory must be freed.
 */
char* convert_genome_to_alpha_copy(const uint8_t* genome, uint32_t size);

/**
 * Set an organism to have the given genome. The genome must be converted.
 * The genome is copied so after this is called it can be freed without
 * affecting the organism. If there is already an organism it is erased/cleared.
 */
void organism_set(organism* org, uint8_t* genome, uint32_t size);

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
void organism_spawn(environment* env, uint8_t* genome, uint32_t size, uint32_t i, uint32_t j);

/**
 * Records that an organism has completed a task. The task value is one of the
 * constants TASK_* such as TASK_OR and TASK_EQUAL.
 * 
 * This function is *NOT* implemented in organism.c. You must implement this
 * method yourself in simulation.c. The serial version is provided and updates
 * array as appropriate along with granting the organism additional energy.
 */
void organism_completed_task(organism* o, int task);

/**
 * Update an organism by having it execute an instruction (if it has enough
 * energy) and then giving it a small amount of energy. The organism may die
 * (if it has grown too old) and it may divide. If the organism is non-existent
 * (it is cleared/the cell is not occupied) this function does nothing but is
 * still safe to call on it.
 */
void organism_update(organism* o);

/**
 * Print out all of the details about an organism. This includes its genome (as
 * an alphebetical string), the locations of the instruction point and heads,
 * the values of the registers, the stacks, recent inputs, number of
 * instructions executed, and number of tasks completed. If the organism is
 * non-existent (it is cleared/the cell is not occupied) "Non-existent" is
 * printed.
 */
void organism_print(const organism* org);

/**
 * Reset an organism. It only keeps:
 *   - genome (however any additional flags are removed)
 *   - energy
 *   - number of tasks completed
 * Everything else is as if it had be cleared and set.
 */
void organism_reset(organism* org);

/**
 * Clear an organism causing it to become non-existent. If this is a cell
 * within an environment, that cell is now un-occupied.
 */
void organism_clear(organism* org);
