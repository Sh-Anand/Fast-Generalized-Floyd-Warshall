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

#ifdef __x86_64__
#include "tsc_x86.h"
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2.7e9
#define CALIBRATE
#define ZERO_PROBABILITY 10 //1/ZERO_PROBABILITY is the probability of an entry in the bit matrix being zero


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

void basic_opt_fw_min_plus(double *C, int n) {
    double c_ik;
    int j0, j1, j2, j3, i_n, k_n;

    for (size_t k = 0; k < n; k++) {
        k_n = k*n;
        for (size_t i = 0; i < n; i++) {
            i_n = i*n;
            c_ik = C[i_n + k];

            for (size_t j = 0; j < n-4; j+=4) {
                j1 = j+1, j2 = j+2, j3 = j+3;
                C[i_n + j0] = min(C[i_n + j0], c_ik + C[k_n + j0]);
                C[i_n + j1] = min(C[i_n + j1], c_ik + C[k_n + j1]);
                C[i_n + j2] = min(C[i_n + j2], c_ik + C[k_n + j2]);
                C[i_n + j3] = min(C[i_n + j3], c_ik + C[k_n + j3]);
            }
            switch(n % 4) {
                case 1:
                    C[i_n + (n-1)] = min(C[i_n + (n-1)], c_ik + C[k_n + (n-1)]);
                    break;
                case 2:
                    C[i_n + (n-2)] = min(C[i_n + (n-2)], c_ik + C[k_n + (n-2)]);
                    C[i_n + (n-1)] = min(C[i_n + (n-1)], c_ik + C[k_n + (n-1)]);
                    break;
                case 3:
                    C[i_n + (n-3)] = min(C[i_n + (n-3)], c_ik + C[k_n + (n-3)]);
                    C[i_n + (n-2)] = min(C[i_n + (n-2)], c_ik + C[k_n + (n-2)]);
                    C[i_n + (n-1)] = min(C[i_n + (n-1)], c_ik + C[k_n + (n-1)]);
                    break;
            }
        }
    }
}

void vect_fw_min_plus(double *C, int n) {
    int i_n = 0;

    double *addr_ij, *addr_ik, *addr_kj;
    __m256d c_ij, c_ik, c_kj, c2, cmp_lt, res;
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            addr_ik = &(C[i_n + k]);
            c_ik = _mm256_load_pd(addr_ik);

            for (size_t j = 0; j < n; j+=4) {
                addr_ij = &(C[i_n + j]);
                addr_kj = &(C[k*n + j]);
                
                c_ij = _mm256_load_pd(addr_ij);
                c_kj = _mm256_load_pd(addr_kj);

                c2 = _mm256_add_pd(c_ik, c_kj);
                // Compute min
                cmp_lt = _mm256_cmp_pd(c_ij, c2, _CMP_LT_OQ);
                res = _mm256_blendv_pd(c2, c_ij, cmp_lt);
                _mm256_store_pd(addr_ij, res);
            }
            i_n += n;
        }
    }
}

void opt_fw_max_min(double *C, int n) {
    int i_n = 0;

    double *c_addr = C, *addr_ij, *addr_ik, *addr_kj;
    __m256d c_ij, c_ik, c_kj, c2, cmp_lt, cmp_gt, res;

    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            size_t j = 0;
            for (; j < n; j+=4) {
                addr_ij = c_addr + sizeof(double) * (i_n + j);
                addr_ik = c_addr + sizeof(double) * (i_n + k);
                addr_kj = c_addr + sizeof(double) * (k*n + j);

                c_ij = _mm256_load_pd(addr_ij);
                c_kj = _mm256_load_pd(addr_kj);
                 
                // Compute min
                cmp_lt = _mm256_cmp_pd(c_ik, c_kj, _CMP_LT_OQ);
                c2 = _mm256_blendv_pd(c_kj, c_ik, cmp_lt);
                   
                // Compute min
                cmp_gt = _mm256_cmp_pd(c_ij, c2, _CMP_GT_OQ);
                res = _mm256_blendv_pd(c2, c_ij, cmp_gt);
                _mm256_store_pd(addr_ij, res);
            }
            for (; j < n; j++) {
                C[i*n + j] = max(C[i*n + j], min(C[i*n + k], C[k*n + j]));
            }
            i_n += n;
        }
    }
}

//extremely inefficient to cast at every iteration, find fix with minimal code duplication!
void fw_or_and(double *C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j] = (uint64_t)C[i*n + j] | ((uint64_t)C[i*n + k] & (uint64_t)C[k*n + j]);
            }
        }
    }
}

void opt_fw_or_and_256(double *C, int n) {
    int i_n = 0;

    double *c_addr = &C[0], *addr_ij, *addr_ik, *addr_kj;
    __m256d c_ij, c_ik, c_kj, c2, cmp_lt, res;
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j+=4) {
                addr_ij = c_addr + sizeof(double) * (i_n + j);
                addr_ik = c_addr + sizeof(double) * (i_n + k);
                addr_kj = c_addr + sizeof(double) * (k*n + j);
                
                c_ij = _mm256_load_pd(addr_ij);
                c_ik = _mm256_load_pd(addr_ik);
                c_kj = _mm256_load_pd(addr_kj);

                c2 = _mm256_add_pd(c_ik, c_kj);
                // Compute min
                cmp_lt = _mm256_cmp_pd(c_ij, c2, _CMP_LT_OQ);
                res = _mm256_blendv_pd(c2, c_ij, cmp_lt);
                _mm256_store_pd(addr_ij, res);
            }
            i_n += n;
        }
    }
}

// void opt_fw_or_and_512(double *C, int n) {
//     int i_n = 0;

//     double *c_addr = &C[0], *addr_ij, *addr_ik, *addr_kj;
//     __m512d c_ij, c_ik, c_kj, c2, cmp_lt, res;
//     for (size_t k = 0; k < n; k++) {
//         for (size_t i = 0; i < n; i++) {
//             size_t j = 0;
//             for (; j < n; j+=8) {
//                 addr_ij = c_addr + sizeof(double) * (i_n + j);
//                 addr_ik = c_addr + sizeof(double) * (i_n + k);
//                 addr_kj = c_addr + sizeof(double) * (k*n + j);
                
//                 c_ij = _mm512_load_pd(addr_ij);
//                 c_ik = _mm512_load_pd(addr_ik);
//                 c_kj = _mm512_load_pd(addr_kj);

//                 c2 = _mm512_and_pd(c_ik, c_kj);
//                 res = _mm512_or_pd(c_ij, c2);
//                 _mm512_store_pd(addr_ij, res);
//             }
//             for(;j < n; j++){
//                 u_int64_t c_ij_bits = *(u_int64_t *)(&C[i_n + j]);
//                 u_int64_t c_ik_bits = *(u_int64_t *)(&C[i_n + k]);
//                 u_int64_t c_kj_bits = *(u_int64_t *)(&C[k*n + j]);
//                 C[i_n + j] = c_ij_bits | (c_ik_bits & c_kj_bits); //maybe is wrong cause it gets converted to double
//             }
//             i_n += n;
//         }
//     }
// }

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

void init_bit_matrix(double *C, int n) {
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            C[i * n + j] = (rand() % ZERO_PROBABILITY);
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

//If you want to use a different timer, simply call this function with a different timer function
//Written like this to enable future extensions
double benchmark(double* C, int n, void (*init_matrix) (double*, int), double (*timer) (double*, int, void (double*, int)), void (*compute)(double*, int)) {
    init_matrix(C, n);
    return timer(C, n, compute);
}

void test(int n, void (*baseline)(double*, int), void (*optimization)(double*, int)) {
    double *C_base = (double *)malloc(n*n*sizeof(double));
    double *C_opt = (double *)malloc(n*n*sizeof(double));
    init_matrices(C_base, C_opt, n);
    // Run baseline function on C
    baseline(C_base, n);
    // Run optimized function on C
    optimization(C_opt, n);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("base[%d][%d] = %lf ", i, j, C_base[n*i + j]);
            printf("opt[%d][%d] = %lf\n", i, j, C_opt[n*i + j]);
        }
    }
    
    // Compare both
    for(int i = 0; i < n; ++i) {
        assert(C_opt[i] == C_base[i]);
    }

    free(C_base);
    free(C_opt);
}

static void (*fw[3]) (double*, int) = {fw_min_plus, fw_or_and, fw_max_min};
static void (*init[3])(double*, int) = {init_matrix, init_bit_matrix, init_matrix};
static char msg[3][10] = {"(min, +)", "(and, or)", "(max, min)"};

int main(int argc, char **argv) {
    if (argc!=3) {printf("usage: FW <n> <fw> (fw = 0,1,2 = (min,plus), (or,and), (max, min))\n"); return -1;}
    int n = atoi(argv[1]);
    int fwi = atoi(argv[2]);
    printf("n=%d \n",n);
    double *C = (double *)malloc(n*n*sizeof(double));

#ifdef __x86_64__
    double r = benchmark(C, n, init[fwi], rdtsc, fw[fwi]);
    printf("%s\n", msg[fwi]);
    printf(" FW : RDTSC instruction:\n %lf cycles measured\n\n", r);
#endif

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            printf("%lf", C[n*i + j]);

    // test(n, fw_min_plus, basic_opt_fw_min_plus);

    free(C);

    return 0;
}