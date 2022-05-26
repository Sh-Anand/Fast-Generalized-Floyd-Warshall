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

#ifdef __x86_64__
#include "tsc_x86.h"
#endif

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2.7e9
#define CALIBRATE
#define ZERO_PROBABILITY 10 //1/ZERO_PROBABILITY is the probability of an entry in the bit matrix being zero
#define EPS  0.000001

void print_bits(u_int64_t d){
    unsigned char * bits = (unsigned char *) & d;
    int i;

    for (i = sizeof (u_int64_t) - 1; i >= 0 ; i--) {
         printf ("%02X ", bits[i]);
    }
    printf ("\n");
}

void print_vector(__m256i v){
    u_int64_t values[4];
    memcpy(values, &v, sizeof(values));
    print_bits( values[3]);
    print_bits( values[2]);
    print_bits( values[1]);
    print_bits( values[0]);
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

void fw_or_and(u_int64_t *C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                //printf("%lx \n", C[i*n + k]);
                C[i*n + j] = C[i*n + j] | (C[i*n + k] & C[k*n + j]);
            }
        }
    }
}

void fw_max_min(double *C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j] = max(C[i*n + j], min(C[i*n + k], C[k*n + j]));
            }
        }
    }
}

//////////////////////VECT MIN PLUS//////////////////////

void vect_fw_min_plus(double *C, int n) {
    assert((u_int64_t) C % 32 == 0);
    double *addr_ij, *addr_ik, *addr_kj;
    __m256d c_ij, c_ik, c_kj, c2, cmp_lt, res;
    for (size_t k = 0; k < n; k++) {
        int i_n = 0;
        for (size_t i = 0; i < n; i++) {
            addr_ik = &(C[i_n + k]);

            c_ik = _mm256_set1_pd(*addr_ik);

            for (size_t j = 0; j < n - 4; j+=4) {
                addr_ij = &(C[i_n + j]);
                addr_kj = &(C[k*n + j]);
                c_ij = _mm256_load_pd(addr_ij);
                c_kj = _mm256_load_pd(addr_kj);

                c2 = _mm256_add_pd(c_ik, c_kj);

                cmp_lt = _mm256_cmp_pd(c_ij, c2, _CMP_LT_OQ);
                res = _mm256_blendv_pd(c2, c_ij, cmp_lt);
                _mm256_store_pd(addr_ij, res);
            }
            i_n += n;
        }
    }
}

void vect_fw_max_min(double *C, int n) {
    assert((u_int64_t) C % 32 == 0);
    double *addr_ij, *addr_ik, *addr_kj;
    __m256d c_ij, c_ik, c_kj, c2, cmp_lt, cmp_gt, res;
    for (size_t k = 0; k < n; k++) {
        int i_n = 0;
        for (size_t i = 0; i < n; i++) {
            addr_ik = &(C[i_n + k]);

            c_ik = _mm256_set1_pd(*addr_ik);

            for (size_t j = 0; j < n - 4; j+=4) {
                addr_ij = &(C[i_n + j]);
                addr_kj = &(C[k*n + j]);
                c_ij = _mm256_load_pd(addr_ij);
                c_kj = _mm256_load_pd(addr_kj);
                cmp_lt = _mm256_cmp_pd(c_ik, c_kj, _CMP_LT_OQ);
                c2 = _mm256_blendv_pd(c_kj, c_ik, cmp_lt);

                cmp_gt = _mm256_cmp_pd(c_ij, c2, _CMP_GT_OQ);
                res = _mm256_blendv_pd(c2, c_ij, cmp_gt);
                _mm256_store_pd(addr_ij, res);
            }
            i_n += n;
        }
    }
}

void vect_fw_or_and(u_int64_t *C,int n) {
    // with assert doesn't work
    //assert((u_int64_t) C % 32 == 0);
    u_int64_t *c_addr = C, *addr_ik, *addr_ij , *addr_kj;
    __m256i c_ij, c_ik, c_kj, c2, cmp_lt, res;
    for (size_t k = 0; k < n; k++) {
        int i_n = 0;
        for (size_t i = 0; i < n; i++) {
            size_t j = 0;
            addr_ik = c_addr + (i_n + k);
            
            c_ik = _mm256_set_epi64x (*addr_ik,
                                        *addr_ik,
                                        *addr_ik,
                                        *addr_ik);
            for (; j < n - 4; j+=4) {
                addr_ij = c_addr + (i_n + j);
                addr_kj = c_addr + (k*n + j);

                c_ij = _mm256_load_si256((__m256i *)addr_ij);
                
                c_kj = _mm256_load_si256((__m256i *)addr_kj);
  
                c2 = _mm256_and_si256(c_ik, c_kj);

                res = _mm256_or_si256(c_ij, c2);
                _mm256_store_si256((__m256i *) addr_ij ,res);
            }
            
            for(;j < n; j++){
                C[i_n + j] = C[i_n + j] | ( C[i_n + k] &  C[k*n + j]); 
            }
            
            i_n += n;
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
double rdtsc_or(uint64_t *C, int n, void (*compute)(uint64_t*, int)) {
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

void bench_test(int n, void (*baseline)(double*, int), void (*optimization)(double*, int)) {
    double *C_base = (double *)malloc(n*n*sizeof(double));
    double *C_opt = (double *)aligned_alloc(32, n*n*sizeof(double));
    init_matrices(C_base, C_opt, n);
    
    double base, opt;
    // Run baseline function on C
    base = rdtsc(C_base, n, baseline);
    printf("Time to run baseline %lf cycles\n", base);
    // Run optimized function on C
    opt = rdtsc(C_opt, n, optimization);
    printf("Time to run optimal  %lf cycles\n", opt);
    /*
    */
    // Compare both
    for(int i = 0; i < n; ++i) {
        assert(abs(C_opt[i] - C_base[i])  <= EPS);
    }

    free(C_base);
    free(C_opt);
}

void bench_test_or(int n, void (*baseline)(uint64_t*, int), void (*optimization)(uint64_t*, int)) {
    uint64_t *C_base = (uint64_t *)malloc(n*n*sizeof(uint64_t));
    uint64_t *C_opt = (uint64_t *)aligned_alloc(32, n*n*sizeof(uint64_t));
    init_bit_matrices(C_base, C_opt, n);
    
    double base, opt;
    // Run baseline function on C
    base = rdtsc_or(C_base, n, baseline);
    printf("Time to run baseline %lf cycles\n", base);
    // Run optimized function on C
    opt = rdtsc_or(C_opt, n, optimization);
    printf("Time to run optimal  %lf cycles\n", opt);
    /*
    */
    // Compare both
    for(int i = 0; i < n; ++i) {
        assert(C_opt[i] == C_base[i]);
    }

    free(C_base);
    free(C_opt);
}


int main(int argc, char **argv) {
    if (argc!=3) {printf("usage: FW <n> <fw> (fw = 0,1,2 = (min,plus), (or,and), (max, min))\n"); return -1;}
    int n = atoi(argv[1]);
    int fwi = atoi(argv[2]);
    printf("n=%d \n",n);

    switch (fwi)
    {
    case 0:
        bench_test(n, fw_min_plus, vect_fw_min_plus);
        break;
    case 1:
        bench_test_or(n, fw_or_and, vect_fw_or_and);
        break;
    case 2:
        bench_test(n, fw_max_min, vect_fw_max_min);
        break;
    default:
        break;
    }

    return 0;
}