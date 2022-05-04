//code largely repurposed from HW1
//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler
#ifdef _WIN32
#include <windows.h> // Include if under windows
#endif

#ifdef linux
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))
#endif

#ifndef WIN32
#include <sys/time.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>

#ifdef __x86_64__
#include "tsc_x86.h"
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2.7e9
#define CALIBRATE
#define ZERO_PROBABILITY 10 //1/ZERO_PROBABILITY is the probability of an entry in the bit matrix being zero

typedef union {
    double d;
    uint64_t u;
} double_bin_t;

double_bin_t *C;

void fw_max_min(int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j].d = max(C[i*n + j].d, min(C[i*n + k].d, C[k*n + j].d));
            }
        }
    }
}

void fw_min_plus(int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j].d = min(C[i*n + j].d, C[i*n + k].d + C[k*n + j].d);
            }
        }
    }
}

void fw_or_and(int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j].u = C[i*n + j].u | (C[i*n + k].u & C[k*n + j].u);
            }
        }
    }
}

void init_matrix(int n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            C[i * n + j].d = ((double )rand() + 1)/RAND_MAX;
        }
    }   
}

void init_bit_matrix(int n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            C[i * n + j].d = (rand() % ZERO_PROBABILITY);
        }
    }   
}

/* 
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
double rdtsc(int n, void (*compute)(int)) {
    int i, num_runs;
    myInt64 cycles;
    myInt64 start;
    num_runs = NUM_RUNS;

    /* 
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to 
     * avoid measurements bias due to the timing overhead.
     */
#ifdef CALIBRATE
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (i = 0; i < num_runs; ++i) {
            compute(n);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        compute(n);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
#endif

//If you want to use a different timer, simply call this function with a different timer function
//Written like this to enable future extensions
double benchmark(int n, void (*init_matrix) (int), double (*timer) (int, void (int)), void (*compute)(int)) {
    init_matrix(n);
    return timer(n, compute);
}

static void (*fw[3]) (int) = {fw_min_plus, fw_or_and, fw_max_min};
static void (*init[3])(int) = {init_matrix, init_bit_matrix, init_matrix};
static char msg[3][10] = {"(min, +)", "(and, or)", "(max, min)"};

int main(int argc, char **argv) {
    if (argc!=3) {printf("usage: FW <n> <fw> (fw = 0,1,2 = (min,plus), (or,and), (max, min))\n"); return -1;}
    int n = atoi(argv[1]);
    int fwi = atoi(argv[2]);
    printf("n=%d \n",n);
    C = (double_bin_t *)malloc(n*n*sizeof(double));

#ifdef __x86_64__
    double r = benchmark(n, init[fwi], rdtsc, fw[fwi]);
    printf("%s\n", msg[fwi]);
    printf(" FW : RDTSC instruction:\n %lf cycles measured\n\n", r);
#endif

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            printf("%lf", C[n*i + j].d);

    return 0;
}