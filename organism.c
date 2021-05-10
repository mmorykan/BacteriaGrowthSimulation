/**
 * Defines the logic for an organism and the environment. Two of the organism
 * functions are not defined here and must be implemented elsewhere. You should
 * not have to modify this file at all.
 */
#include "organism.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>


//////////////////////////////
///// Internal Constants /////
//////////////////////////////

// Organism Properties and Constants
#define ENERGY_PER_INSTRUCTION 64 // how much energy an organism requires to execute a single instruction (note: energy is very different than the original implementation)
#define AGE_LIMIT_FACTOR 20       // once an organism has executed this many times its genome length without reproducing it dies
#define MIN_GENOME_LENGTH 8       // no organism will be created with a genome smaller than this
#define OFFSPRING_SIZE_RANGE 1.0  // range of "viable" offspring sizes
#define NUM_RECENT_INPUT 3        // number of recent inputs to consider when checking reccent tasks
const static int32_t INPUTS[] = {0x0FA8458D, 0x33B49485, 0x5590312A}; // fixed set of inputs that is cycled through

// Available CPU Instructions and memory flags
#define INSTRUCTION_NOP_A    0 // a
#define INSTRUCTION_NOP_B    1 // b
#define INSTRUCTION_NOP_C    2 // c
#define INSTRUCTION_IF_N_EQU 3 // d
#define INSTRUCTION_IF_LESS  4 // e
#define INSTRUCTION_IF_LABEL 5 // f
#define INSTRUCTION_MOV_HEAD 6 // g
#define INSTRUCTION_JMP_HEAD 7 // h
#define INSTRUCTION_GET_HEAD 8 // i
#define INSTRUCTION_SET_FLOW 9 // j
#define INSTRUCTION_SHIFT_R  10 // k
#define INSTRUCTION_SHIFT_L  11 // l
#define INSTRUCTION_INC      12 // m
#define INSTRUCTION_DEC      13 // n
#define INSTRUCTION_PUSH     14 // o
#define INSTRUCTION_POP      15 // p
#define INSTRUCTION_SWAP_STK 16 // q
#define INSTRUCTION_SWAP     17 // r
#define INSTRUCTION_ADD      18 // s
#define INSTRUCTION_SUB      19 // t
#define INSTRUCTION_NAND     20 // u
#define INSTRUCTION_H_COPY   21 // v
#define INSTRUCTION_H_ALLOC  22 // w
#define INSTRUCTION_H_DIVIDE 23 // x
#define INSTRUCTION_IO       24 // y
#define INSTRUCTION_H_SEARCH 25 // z
#define NUM_INSTRUCTIONS     26

#define INSTRUCTIONS_MASK    0x1F // instructions will only use these bits, flags can use all other bits
#define FLAG_EXECUTED        0x40
#define FLAG_COPIED          0x80


//////////////////////////
///// Utility Macros /////
//////////////////////////
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define NAND(a, b) ~((a)&(b))
#define ARRAY_SIZE(a) (sizeof(a)/sizeof(*(a)))
#define RECENT_INPUT(n, i) INPUTS[((n)+16*ARRAY_SIZE(INPUTS)-(i)-1)%ARRAY_SIZE(INPUTS)]
#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)


///////////////////////////////////////////////
///// Core Environment/Organism Functions /////
///////////////////////////////////////////////

bool environment_load_from_path(environment* env, const char* path) {
    FILE* file = fopen(path, "r");
    if (!file) { perror("error opening file for reading"); return false; }
    bool retval = environment_load_from_file(env, file);
    fclose(file);
    return retval;
}

bool environment_load_from_file(environment* env, FILE* file) {

    // read the world size and allocate the grid
    size_t world_size;
    if (fscanf(file, "%zu", &world_size) != 1) { perror("error reading world-size"); return false; }
    env->world_size = world_size;
    environment_alloc_grid(env);

    // read and set all of the organisms
    while (true) {
        size_t i, j;
        char genome_alpha[MAX_GENOME_LENGTH+1];
        uint8_t genome[MAX_GENOME_LENGTH];
        int matches = fscanf(file, "%zu %zu %" STRINGIFY(MAX_GENOME_LENGTH) "[abcdefghijklmnopqrstuvwxyz] ",
                            &i, &j, genome_alpha);
        if (matches != 3) {
            if (feof(file) && (matches == 0 || matches == EOF)) { return true; } // finished reading the file successfully!
            perror("error reading organism");
            break;
        }
        if (i >= world_size || j >= world_size) {
            fprintf(stderr, "organism location outside the world: %zu >= %zu || %zu >= %zu\n", i, world_size, j, world_size);
            break;
        }
        size_t size = strlen(genome_alpha);
        if (size < MIN_GENOME_LENGTH) {
            fprintf(stderr, "organism genome too small: %zu < %u\n", size, MIN_GENOME_LENGTH);
            break;
        }
        convert_genome(genome_alpha, genome, size);
        organism_set(&env->grid[i*world_size+j], genome, size);
    }

    // failed
    free(env->grid);
    return false;
}

bool environment_save_to_path(environment* env, const char* path) {
    FILE* file = fopen(path, "w");
    if (!file) { perror("error opening file for writing"); return false; }
    bool retval = environment_save_to_file(env, file);
    fclose(file);
    return retval;
}

bool environment_save_to_file(environment* env, FILE* file) {
    // write the world size
    size_t world_size = env->world_size;
    fprintf(file, "%zu\n", env->world_size);

    // write all of the organisms
    for (size_t i = 0; i < world_size; i++) {
        for (size_t j = 0; j < world_size; j++) {
            organism* org = &env->grid[i*world_size + j];
            if (cell_occupied(org)) {
                char genome_alpha[MAX_GENOME_LENGTH+1];
                convert_genome_to_alpha(org->memory, genome_alpha, MIN(org->genome_size, org->size));
                fprintf(file, "%zu %zu %s\n", i, j, genome_alpha);
            }
        }
    }
    return true;
}

void environment_alloc_grid(environment* env) {
    size_t size = env->world_size;
    size_t nbytes = size*size*sizeof(organism);
    organism* grid = env->grid = (organism*)memset(malloc(nbytes), 0, nbytes);
    for (uint32_t i = 0; i < size; i++) {
        for (uint32_t j = 0; j < size; j++) {
            grid[i*size+j].env = env;
            grid[i*size+j].location_i = i;
            grid[i*size+j].location_j = j;
        }
    }
}

void environment_clear_tasks(environment* env) {
    memset(env->tasks_completed, 0, sizeof(env->tasks_completed));
}

bool cell_occupied(organism* org) { return org->size != 0; }

void convert_genome(const char* genome_alpha, uint8_t* genome, uint32_t size) {
    for (uint32_t i = 0; i < size; i++) { genome[i] = genome_alpha[i] - 'a'; }
}

uint8_t* convert_genome_copy(const char* genome_alpha, uint32_t size) {
    uint8_t* output = (uint8_t*)malloc(size);
    convert_genome(genome_alpha, output, size);
    return output;
}

void convert_genome_to_alpha(const uint8_t* genome, char* genome_alpha, uint32_t size) {
    for (uint32_t i = 0; i < size; i++) { genome_alpha[i] = (genome[i] & INSTRUCTIONS_MASK) + 'a'; }
    genome_alpha[size] = 0; // null string terminator
}

char* convert_genome_to_alpha_copy(const uint8_t* genome, uint32_t size) {
    char* output = (char*)malloc(size+1);
    convert_genome_to_alpha(genome, output, size);
    return output;
}

void organism_clear(organism* org) {
    if (org->size) {
        environment* env = org->env;
        uint32_t i = org->location_i, j = org->location_j;
        memset(org, 0, sizeof(organism));
        org->env = env;
        org->location_i = i;
        org->location_j = j;
    }
}

void organism_set(organism* org, uint8_t* genome, uint32_t size) {
    organism_clear(org);
    for (uint32_t i = 0; i < size; i++) { org->memory[i] = genome[i] & INSTRUCTIONS_MASK; }
    org->size = org->genome_size = size;
}

void organism_reset(organism* org) {
    if (org->size) {
        memset(org->registers, 0, sizeof(org->registers));
        memset(org->heads, 0, sizeof(org->heads));
        memset(org->stacks, 0, sizeof(org->stacks));
        memset(org->stack_sizes, 0, sizeof(org->stack_sizes));
        memset(org->recent_copies, 0, sizeof(org->recent_copies));
        for (uint32_t i = 0; i < org->size; i++) { org->memory[i] &= INSTRUCTIONS_MASK; }
        org->genome_size = org->size;
        org->stack_active = 0;
        org->ip = org->size - 1;
        org->instructions_executed = 0;
        org->num_inputs = 0;
        org->has_allocated = false;
        org->recent_copies_offset = 0;
        org->recent_copies_size = 0;
    }
}

void organism_print(const organism* org) {
    if (!org->size) { printf("Non-existent\n"); return; }
    printf("Genome: ");
    for (uint32_t i = 0; i < org->size; i++) {
        if (i == org->ip) { printf("I"); }
        else if (i == org->flow_head) { printf("F"); }
        else if (i == org->write_head) { printf("W"); }
        else if (i == org->read_head) { printf("R"); }
        else { printf(" "); }
    }
    printf("\n");
    printf("IP: %u, Read Head: %u, Write Head: %u, Flow Head: %u\n", org->ip, org->read_head, org->write_head, org->flow_head);
    printf("Size: %u, Instructions Executed: %u\n", org->size, org->instructions_executed);
    for (int r = 0; r < NUM_REGISTERS; r++) { printf("%cX: %d (0x%x)%s", r + 'A', org->registers[r], org->registers[r], r != NUM_REGISTERS-1 ? ", " : "\n"); }
    for (int s = 0; s < NUM_STACKS; s++) {
        printf("%cStack %c:", org->stack_active == s ? '*' : ' ', s + 'A');
        for (int i = 0; i < org->stack_sizes[s]; i++) { printf(" %d (0x%x)", org->stacks[s][i], org->stacks[s][i]); }
        printf("\n");        
    }
    printf("Recent Input:");
    size_t ni = org->num_inputs;
    for (int i = 0; i < NUM_RECENT_INPUT; i++) { printf(" %d (0x%x)", RECENT_INPUT(ni, i), RECENT_INPUT(ni, i)); }
    printf("\n");
    printf("Tasks Completed:");
    for (int i = 0; i < NUM_TASKS; i++) { printf(" %u", org->tasks_completed[i]); }
    printf("\n");
}

void environment_print(const environment* env) {
    size_t size = env->world_size;
    uint32_t limit = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            organism* org = &env->grid[i*size+j];
            size_t n = 0;
            for (int k = 0; k < NUM_TASKS; k++) { n += org->tasks_completed[k]; }
            if (n > limit) { limit = n; }
        }
    }

    const static char* BLOCKS[] = {"·", "░", "▒", "▓", "█"};
    printf("┏");
    for (int j = 0; j < size; j++) { printf("━"); }
    printf("┓\n");
    for (int i = 0; i < size; i++) {
        printf("┃");
        for (int j = 0; j < size; j++) {
            organism* org = &env->grid[i*size+j];
            if (org->size) {
                uint32_t level = 0;
                if (limit != 0) {
                    size_t n = 0;
                    for (int k = 0; k < NUM_TASKS; k++) { n += org->tasks_completed[k]; }
                    level = (n*4+limit-1)/limit;
                }
                printf("%s", BLOCKS[level]);
            } else {
                printf(" ");
            }
        }
        printf("┃\n");
    }
    printf("┗");
    for (int j = 0; j < size; j++) { printf("━"); }
    printf("┛\n");
}


///////////////////////////////////////////
///// CPU Instruction Implementations /////
///////////////////////////////////////////main

///// Helper functions /////
static inline bool _memcmp_cyclic(uint8_t* buf, uint32_t off, uint32_t buf_size, uint8_t* match, uint32_t match_size) {
    if (match_size == 0) { return true; } // matches behavior of memcmp
    for (uint32_t i = 0; i < match_size; i++) {
        if ((buf[(i+off)%buf_size] & INSTRUCTIONS_MASK) != (match[i] & INSTRUCTIONS_MASK)) { return false; }
    }
    return true;
}
static inline uint8_t _get_next_instruction(organism* o) { return o->memory[(o->ip+1)%o->size] & INSTRUCTIONS_MASK; }
static inline void _advance(organism* o, uint32_t *head) { *head += 1; if (*head == o->size) { *head = 0; } }
static inline void _advance_ip(organism* o) { _advance(o, &o->ip); }
static inline void _retreat_ip(organism* o) { o->ip = (o->ip == 0 ? o->size : o->ip) - 1; }
static inline bool _is_nop(uint8_t instruction) { return instruction - INSTRUCTION_NOP_A < NUM_REGISTERS; }
static inline uint8_t _convert_nop_to_register(uint8_t instruction) { return instruction - INSTRUCTION_NOP_A + AX; }
static inline uint8_t _convert_nop_to_head(uint8_t instruction) { return instruction - INSTRUCTION_NOP_A + IP; }
static inline uint8_t _get_dynamic_register(organism* o, uint8_t def) {
    uint8_t next = _get_next_instruction(o);
    if (_is_nop(next)) { _advance_ip(o); return _convert_nop_to_register(next); }
    return def;
}
static inline uint8_t _get_register_complement(uint8_t reg) { return reg == NUM_REGISTERS-1 ? AX : reg + 1; }
static inline uint8_t _get_dynamic_head(organism* o, uint8_t def) {
    uint8_t next = _get_next_instruction(o);
    if (_is_nop(next)) { _advance_ip(o); return _convert_nop_to_head(next); }
    return def;
}
static inline uint32_t _get_template_complement(organism* o, uint8_t* template) {
    uint32_t size = 0;
    uint8_t next = _get_next_instruction(o);
    while (size < MAX_TEMPLATE_LEN && _is_nop(next)) {
        _advance_ip(o);
        uint8_t com = _get_register_complement(next - INSTRUCTION_NOP_A) + INSTRUCTION_NOP_A;
        template[size++] = com;
        next = _get_next_instruction(o);
    }
    return size;
}
static inline uint32_t _count_flagged_mem(organism* o, uint32_t start, uint32_t end, uint8_t flag) {
    uint32_t count = 0;
    for (uint32_t i = start; i < end; i++) { count += (o->memory[i] & flag) != 0; }
    return count;
}

#define REG_VAL_DYN o->registers[_get_dynamic_register(o, BX)]   // dynamic register, default BX

///// NOP Instructions /////
static inline void inst_nop(organism* o) { } // does nothing, works for a, b, and c

///// Conditional Instructions /////
static inline void inst_if_n_equ(organism* o) {
    uint8_t reg = _get_dynamic_register(o, BX);
    uint8_t com = _get_register_complement(reg);
    if (o->registers[reg] == o->registers[com]) { _advance_ip(o); } // skip next
}
static inline void inst_if_less(organism* o) {
    uint8_t reg = _get_dynamic_register(o, BX);
    uint8_t com = _get_register_complement(reg);
    if (o->registers[reg] >= o->registers[com]) { _advance_ip(o); } // skip next
}
static inline void inst_if_label(organism* o) {
    uint8_t template[MAX_TEMPLATE_LEN];
    uint32_t template_size = _get_template_complement(o, template);
    if (template_size != o->recent_copies_size ||
        !_memcmp_cyclic(o->recent_copies, o->recent_copies_offset, template_size, template, template_size)) { _advance_ip(o); } // skip next
}

///// Flow Control Instructions /////
static inline void inst_mov_head(organism* o) {
    uint8_t head = _get_dynamic_head(o, IP);
    o->heads[head] = o->flow_head;
    if (head == IP) { _retreat_ip(o); } // counteract the upcoming _advance_ip();
}
static inline void inst_jmp_head(organism* o) {
    uint8_t head = _get_dynamic_head(o, IP);
    int32_t pos = (o->heads[head] + o->cx) % o->size;
    if (pos < 0) { pos = 0; } // NOTE: this is the official implementation, but I would personally do while (pos < 0) { pos += o->size }
    o->heads[head] = pos;
    // if (head == IP) { _retreat_ip(o); } // counteract the upcoming _advance_ip(); // NOTE: official implementation doesn't have this but they also have a comment asking about if they should
}
static inline void inst_get_head(organism* o) { o->cx = o->heads[_get_dynamic_head(o, IP)]; }
static inline void inst_set_flow(organism* o) {
    int32_t pos = o->registers[_get_dynamic_register(o, CX)] % o->size;
    if (pos < 0) { pos += o->size; }
    o->flow_head = pos;
}

///// Arithmetic/Logic Instructions /////
static inline void inst_shift_r(organism* o) { REG_VAL_DYN >>= 1; }
static inline void inst_shift_l(organism* o) { REG_VAL_DYN <<= 1; }
static inline void inst_inc(organism* o) { REG_VAL_DYN++; }
static inline void inst_dec(organism* o) { REG_VAL_DYN--; }
static inline void inst_add(organism* o) { REG_VAL_DYN = o->bx + o->cx; }
static inline void inst_sub(organism* o) { REG_VAL_DYN = o->bx - o->cx; }
static inline void inst_nand(organism* o) { REG_VAL_DYN = NAND(o->bx, o->cx); }

///// Data Instructions /////
#define STACK_ACTIVE o->stacks[o->stack_active]
#define STACK_SIZE_ACTIVE o->stack_sizes[o->stack_active]
static inline void inst_push(organism* o) {
    STACK_ACTIVE[STACK_SIZE_ACTIVE++] = REG_VAL_DYN;
    if (STACK_SIZE_ACTIVE == MAX_STACK_SIZE) { STACK_SIZE_ACTIVE = MAX_STACK_SIZE - 1; }
}
static inline void inst_pop(organism* o) { REG_VAL_DYN = STACK_SIZE_ACTIVE ? STACK_ACTIVE[--STACK_SIZE_ACTIVE] : 0; }
static inline void inst_swap_stack(organism* o) { o->stack_active = (o->stack_active + 1) % NUM_STACKS; }
static inline void inst_swap(organism* o) {
    uint8_t reg = _get_dynamic_register(o, BX);
    uint8_t com = _get_register_complement(reg);
    int32_t temp = o->registers[reg];
    o->registers[reg] = o->registers[com];
    o->registers[com] = temp;
}

///// IO/Energy Instruction /////
// Note: this is somewhat different than the official implementation
static inline void inst_io(organism* o) {
    uint8_t reg = _get_dynamic_register(o, BX);

    // get the produced output and check it for a task that earns energy
    uint32_t out = o->registers[reg];
    uint32_t in2 = RECENT_INPUT(o->num_inputs, 0);
    for (int i = 0; i < NUM_RECENT_INPUT-1; i++) {
        uint32_t in1 = in2; in2 = RECENT_INPUT(o->num_inputs, i+1);
        if ((~(in1^in2)) == out) { organism_completed_task(o, TASK_EQUAL); break; }
        else if ((in1^in2) == out) { organism_completed_task(o, TASK_XOR); break; }
        else if ((~(in1|in2)) == out) { organism_completed_task(o, TASK_NOR); break; }
        else if ((in1&~in2) == out) { organism_completed_task(o, TASK_AND_NOT); break; }
        else if ((in1|in2) == out) { organism_completed_task(o, TASK_OR); break; }
        else if ((in1|~in2) == out) { organism_completed_task(o, TASK_OR_NOT); break; }
        else if ((in1&in2) == out) { organism_completed_task(o, TASK_AND); break; }
        else if (NAND(in1,in2) == out) { organism_completed_task(o, TASK_NAND); break; }
        else if (i == 0 && ~in1 == out) { organism_completed_task(o, TASK_NOT); break; }
    }

    // get a new input
    o->registers[reg] = RECENT_INPUT(++o->num_inputs, 0);
}

///// Head/Lifecycle Instructions /////
static inline void inst_h_copy(organism* o) {
    // copy (and possibly mutate)
    uint8_t orig = o->memory[o->read_head] & INSTRUCTIONS_MASK;
    uint8_t value = (rand()/(float)RAND_MAX <= o->env->mutation_prob) ? rand() % NUM_INSTRUCTIONS : orig;
    o->memory[o->write_head] = value | FLAG_COPIED;

    // save recent copies of NOPs (and reset for non-nops)
    //value = orig; // official implementation can have the recent copy data also "mutated" causing issues with some other instructions
    if (_is_nop(value)) {
        o->recent_copies[o->recent_copies_offset++] = value;
        if (o->recent_copies_offset == MAX_TEMPLATE_LEN) { o->recent_copies_offset = 0; }
        if (o->recent_copies_size < MAX_TEMPLATE_LEN) { o->recent_copies_size++; }
    } else { o->recent_copies_offset = o->recent_copies_size = 0; }

    // advance
    _advance(o, &o->read_head);
    _advance(o, &o->write_head);
}
static inline void inst_h_alloc(organism* o) {
    if (o->has_allocated) { return; } // h_alloc cannot be called twice before h_divide
    uint32_t old_size = o->size;
    uint32_t new_size = o->size = MIN(MAX_GENOME_LENGTH, old_size*2);
    memset(o->memory+old_size, 0, new_size-old_size);
    o->registers[_get_dynamic_register(o, AX)] = old_size;
    o->has_allocated = true;
}
static inline void inst_h_divide(organism* o) { // NOTE: this is the only CPU instruction that truly requires the environment or location...
    if (!o->has_allocated) { return; } // h_alloc must be called before h_divide

    // get the region that will become the offspring
    const uint32_t start = o->read_head;
    const uint32_t end = o->write_head ? o->write_head : o->size;
    const int offspring_size = end - start;

    // compute basic viability
    const int min_size = MAX(MIN_GENOME_LENGTH, o->genome_size/OFFSPRING_SIZE_RANGE);
    const int max_size = MIN(MAX_GENOME_LENGTH, o->genome_size*OFFSPRING_SIZE_RANGE);
    const int executed_size = _count_flagged_mem(o, 0, start, FLAG_EXECUTED);
    const int copied_size = _count_flagged_mem(o, start, end, FLAG_COPIED);
    if (offspring_size < min_size || start < min_size ||
        offspring_size > max_size || start > max_size ||
        copied_size < 0.5*offspring_size || executed_size < 0.5*start) { return; }

    // figure out where to place the offspring
    size_t world_size = o->env->world_size;
    uint32_t i = o->location_i, j = o->location_j;
    int32_t i_off, j_off;
    do {
        i_off = rand() % 3 - 1;
        j_off = rand() % 3 - 1;
    } while ((i_off == 0 && j_off == 0) ||
             (i == 0 && i_off == -1) || (i + i_off >= world_size) ||
             (j == 0 && j_off == -1) || (j + j_off >= world_size));

    // spawn the offspring
    organism_spawn(o->env, o->memory+start, offspring_size, i+i_off, j+j_off);

    // remove the sequence from the parent
    o->size = start;

    // reset the parent
    organism_reset(o);
}

static inline void inst_h_search(organism* o) {
    uint32_t orig_ip = o->ip;
    uint8_t template[MAX_TEMPLATE_LEN];
    uint32_t template_size = _get_template_complement(o, template);
    uint32_t next_idx = (o->ip + 1) % o->size;
    o->cx = template_size;
    if (template_size > 0) {
        // perform search
        uint32_t index = next_idx;
        while (index != orig_ip) {
            if (_memcmp_cyclic(o->memory, index, o->size, template, template_size)) {
                // matches!
                o->bx = index - orig_ip - 1;
                o->flow_head = (index + template_size) % o->size;
                return;
            }
            index = (index + 1) % o->size;
        }
    }
    // no match found
    o->bx = 0;
    o->flow_head = next_idx;
}
///// Execute the next CPU instruction for the organism /////
void execute_instruction(organism* o) {
    uint8_t instruction = o->memory[o->ip] & INSTRUCTIONS_MASK;
    if (instruction == INSTRUCTION_NOP_A) { inst_nop(o); }
    else if (instruction == INSTRUCTION_NOP_B) { inst_nop(o); }
    else if (instruction == INSTRUCTION_NOP_C) { inst_nop(o); }
    else if (instruction == INSTRUCTION_IF_N_EQU) { inst_if_n_equ(o); }
    else if (instruction == INSTRUCTION_IF_LESS) { inst_if_less(o); }
    else if (instruction == INSTRUCTION_IF_LABEL) { inst_if_label(o); }
    else if (instruction == INSTRUCTION_MOV_HEAD) { inst_mov_head(o); }
    else if (instruction == INSTRUCTION_JMP_HEAD) { inst_jmp_head(o); }
    else if (instruction == INSTRUCTION_GET_HEAD) { inst_get_head(o); }
    else if (instruction == INSTRUCTION_SET_FLOW) { inst_set_flow(o); }
    else if (instruction == INSTRUCTION_SHIFT_R) { inst_shift_r(o); }
    else if (instruction == INSTRUCTION_SHIFT_L) { inst_shift_l(o); }
    else if (instruction == INSTRUCTION_INC) { inst_inc(o); }
    else if (instruction == INSTRUCTION_DEC) { inst_dec(o); }
    else if (instruction == INSTRUCTION_PUSH) { inst_push(o); }
    else if (instruction == INSTRUCTION_POP) { inst_pop(o); }
    else if (instruction == INSTRUCTION_SWAP_STK) { inst_swap_stack(o); }
    else if (instruction == INSTRUCTION_SWAP) { inst_swap(o); }
    else if (instruction == INSTRUCTION_ADD) { inst_add(o); }
    else if (instruction == INSTRUCTION_SUB) { inst_sub(o); }
    else if (instruction == INSTRUCTION_NAND) { inst_nand(o); }
    else if (instruction == INSTRUCTION_H_COPY) { inst_h_copy(o); }
    else if (instruction == INSTRUCTION_H_ALLOC) { inst_h_alloc(o); }
    else if (instruction == INSTRUCTION_H_DIVIDE) { inst_h_divide(o); }
    else if (instruction == INSTRUCTION_IO) { inst_io(o); }
    else if (instruction == INSTRUCTION_H_SEARCH) { inst_h_search(o); }
    else { fprintf(stderr, "invalid instruction: 0x%02x\n", instruction); }
    o->memory[o->ip] |= FLAG_EXECUTED;
    _advance_ip(o);
    o->instructions_executed++;
}


///// Execute the next CPU instruction for the organism /////
void organism_update(organism* o) {
    if (!o->size) { return; } // non-existent
    if (o->energy >= ENERGY_PER_INSTRUCTION) {
        execute_instruction(o);
        o->energy -= ENERGY_PER_INSTRUCTION;
        if (o->instructions_executed >= AGE_LIMIT_FACTOR * o->genome_size) { organism_clear(o); }
    }
    o->energy++;
}
