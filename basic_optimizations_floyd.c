//code largely repurposed from HW1
//#error Please comment out the next two lines under linux, then comment this error
//#include "stdafx.h"  //Visual studio expects this line to be the first one, comment out if different compiler

#ifdef linux
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))
#endif

// #ifndef WIN32
// #include <sys/time.h>
// #include <windows.h> 
// #endif
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <immintrin.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#ifdef __x86_64__
#include "tsc_x86.h"
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2.7e9
#define CALIBRATE
#define ZERO_PROBABILITY 10 //1/ZERO_PROBABILITY is the probability of an entry in the bit matrix being zero
#define EPS  0.000001


//============================================================================================
// =============================== BASE IMPLEMENTATIONS ======================================
//============================================================================================

void fw_max_min(double *C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j] = max(C[i*n + j], min(C[i*n + k], C[k*n + j]));
            }
        }
    }
}

void fw_min_plus(double *C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j] = min(C[i*n + j], C[i*n + k] + C[k*n + j]);
            }
        }
    }
}

void fw_or_and_int(u_int64_t *C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                //printf("%lx \n", C[i*n + k]);
                C[i*n + j] = C[i*n + j] | (C[i*n + k] & C[k*n + j]);
            }
        }
    }
}

//============================================================================================
// =============================== BASIC OPTIMIZATIONS =======================================
//============================================================================================

void basic_optimization_to_plot_max_min(double *C, int n) {
    int in = 0;
    int kn = 0;
    int inpj = 0;
    int inpk = 0;
    int knpj = 0;

    double c_in_p_k = 0.0;
    double c_in_p_j = 0.0;
    double c_kn_p_j = 0.0;
    double min_c = 0.0;
    double max_c = 0.0;

    for (size_t k = 0; k < n; k++) {
        kn = k * n;
        for (size_t i = 0; i < n; i++) {
            in = i * n;
            inpk = in + k;
            c_in_p_k = C[inpk];
            for (size_t j = 0; j < n; j++) {
                inpj = in + j;
                knpj = kn + j;
                c_kn_p_j = C[knpj];
                c_in_p_j = C[inpj];
                min_c = min(c_in_p_k, c_kn_p_j);
                max_c =  max(c_in_p_j, min_c);
                C[inpj] = max_c;
            }
        }
    }
}

void basic_optimization_to_plot_or_and(u_int64_t *C, int n) {
    int in = 0;
    int kn = 0;
    int inpj = 0;
    int inpk = 0;
    int knpj = 0;

    u_int64_t c_in_p_k = 0.0;
    u_int64_t c_in_p_j = 0.0;
    u_int64_t c_kn_p_j = 0.0;
    u_int64_t or_c = 0.0;
    u_int64_t and_c = 0.0;

    for (size_t k = 0; k < n; k++) {
        kn = k * n;
        for (size_t i = 0; i < n; i++) {
            in = i * n;
            inpk = in + k;
            c_in_p_k = C[inpk];
            for (size_t j = 0; j < n; j++) {
                inpj = in + j;
                knpj = kn + j;
                c_kn_p_j = C[knpj];
                c_in_p_j = C[inpj];
                and_c = c_in_p_k & c_kn_p_j;
                or_c =  c_in_p_j | and_c;
                C[inpj] = or_c;
            }
        }
    }
}

void basic_optimization_to_plot_min_plus(double *C, int n) {
    int in = 0;
    int kn = 0;
    int inpj = 0;
    int inpk = 0;
    int knpj = 0;

    double c_in_p_c_kn = 0.0;
    double c_in_p_k = 0.0;
    double c_in_p_j = 0.0;
    double c_kn_p_j = 0.0;
    double min_c = 0.0;

    for (size_t k = 0; k < n; k++) {
        kn = k * n;
        for (size_t i = 0; i < n; i++) {
            in = i * n;
            inpk = in + k;
            c_in_p_k = C[inpk];
            for (size_t j = 0; j < n; j++) {
                inpj = in + j;
                knpj = kn + j;
                c_kn_p_j = C[knpj];
                c_in_p_j = C[inpj];
                c_in_p_c_kn = c_in_p_k + c_kn_p_j;
                min_c =  min(c_in_p_j, c_in_p_c_kn);
                C[inpj] = min_c;
            }
        }
    }
}

//====================================================================================================

void init_matrix(double *C, int n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            C[i * n + j] = ((double )rand() + 1)/RAND_MAX;
        }
    }   
}

void init_matrices(double *C1, double *C2, int n) {
    double x;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            x = ((double )rand() + 1)/RAND_MAX;
            C1[i * n + j] = x;
            C2[i * n + j] = x;
        }
    }
}

void init_bit_matrix(uint64_t *C, int n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            C[i * n + j] = (rand() % ZERO_PROBABILITY);
        }
    }
}

void init_bit_matrices(u_int64_t *C1, u_int64_t *C2, int n) {
    u_int64_t x;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            x = (rand() % ZERO_PROBABILITY);
            C1[i * n + j] = x;
            C2[i * n + j] = x;
        }
    }
}

/* 
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
double rdtsc(double *C, int n, void (*compute)(double*, int)) {
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
            compute(C, n);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        compute(C, n);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
#endif

#ifdef __x86_64__
double rdtsc_int(uint64_t *C, int n, void (*compute)(uint64_t*, int)) {
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
            compute(C, n);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        compute(C, n);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
#endif

static void (*fw[2]) (double*, int) = {fw_min_plus, fw_max_min};
static void (*fw_intermediate[2]) (double*, int) = {basic_optimization_to_plot_min_plus,
                                                    basic_optimization_to_plot_max_min};


double benchmark(int fw_idx, int n) {
    double time = 0;
    if(fw_idx == 2) {
        uint64_t *C = (uint64_t *)malloc(n*n*sizeof(uint64_t));
        init_bit_matrix(C, n);
        time = rdtsc_int(C, n, basic_optimization_to_plot_or_and);
        free(C);
    } else {
        double *C = (double *)malloc(n*n*sizeof(double));
        init_matrix(C, n);
        time = rdtsc(C, n, fw_intermediate[fw_idx]);
        free(C);
    }
    return time;
}

void test(int n, void (*baseline)(double*, int), void (*optimization)(double*, int)) {
    double *C_base = (double *)malloc(n*n*sizeof(double));
    double *C_opt = (double *)aligned_alloc(32, n*n*sizeof(double));
    init_matrices(C_base, C_opt, n);

//    printf("\nBEFORE\n");
//    for(int i = 0; i < n; i++) {
//        for(int j = 0; j < n; j++) {
//            printf("base[%d][%d] = %lf ", i, j, C_base[n*i + j]);
//            printf("opt[%d][%d] = %lf\n", i, j, C_opt[n*i + j]);
//        }
//    }

    // Run baseline function on C
    baseline(C_base, n);
    // Run optimized function on C
    optimization(C_opt, n);

//    printf("\nAFTER\n");
//    for(int i = 0; i < n; i++) {
//        for(int j = 0; j < n; j++) {
//            printf("base[%d][%d] = %lf ", i, j, C_base[n*i + j]);
//            printf("opt[%d][%d] = %lf\n", i, j, C_opt[n*i + j]);
//        }
//    }

    // Compare both
    for(int i = 0; i < n; ++i) {
        assert(fabs(C_opt[i] - C_base[i]) <= EPS);
    }

    free(C_base);
    free(C_opt);
}

void test_or_and(int n){
    u_int64_t *C_base = (u_int64_t *)malloc(n*n*sizeof(double));
    u_int64_t *C_opt = (u_int64_t *)malloc(n*n*sizeof(double));
    init_bit_matrices(C_base, C_opt, n);

//    printf("\nBEFORE\n");
//    for(int i = 0; i < n; i++) {
//        for(int j = 0; j < n; j++) {
//            printf("base[%d][%d] = %lu ", i, j, C_base[n*i + j]);
//            printf("opt[%d][%d] = %lu\n", i, j, C_opt[n*i + j]);
//        }
//    }

    // Run baseline function on C
    fw_or_and_int(C_base, n);
    // Run optimized function on C
    basic_optimization_to_plot_or_and(C_opt, n);

//    printf("\nAFTER\n");
//    for(int i = 0; i < n; i++) {
//        for(int j = 0; j < n; j++) {
//            printf("base[%d][%d] = %lu ", i, j, C_base[n*i + j]);
//            printf("opt[%d][%d] = %lu\n", i, j, C_opt[n*i + j]);
//        }
//    }

    // Compare both
    for(int i = 0; i < n; ++i) {
        assert(C_opt[i] == C_base[i]);
    }

    free(C_base);
    free(C_opt);
}


int main(int argc, char **argv) {
    if (argc!=3) {
        printf("usage: FW <n> <fw> (fw = 0,1,2 = (min,plus), (max, min), (or,and))\n");
        return -1;
    }

    int n = atoi(argv[1]);
    int fw_idx = atoi(argv[2]);
    if (fw_idx < 0 || fw_idx > 2) {
        printf("usage: FW <n> <fw> (fw = 0,1,2 = (min,plus), (max, min), (or,and))\n");
        return -1;
    }

#ifdef __x86_64__
    double r  = benchmark(fw_idx, n);
    printf(" %lf", r);
#endif

    // Compare our results with the base implementation
//    test(n, fw_min_plus, basic_optimization_to_plot_min_plus);
//    test_or_and(n);
//    test(n, fw_max_min, basic_optimization_to_plot_max_min);

    return 0;
}
