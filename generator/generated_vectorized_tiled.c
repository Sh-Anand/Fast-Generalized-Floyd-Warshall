
#ifdef linux
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))
#endif

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

void fw_abc_or_and(uint64_t* A, uint64_t* B, uint64_t* C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j] = C[i*n + j] | (A[i*n + k] & B[k*n + j]);
            }
        }
    }
}

void opt_tiled(uint64_t* A, uint64_t* B, uint64_t* C, int L1, int n, int Bi, int Bj, int Bk) {
    int mm = n / L1;
    for(int k = 0; k < mm; ++k) {
        int l1 = k;
        int m1 = k;
        int sub_base_l = l1 * L1;
        int sub_base_m = m1 * L1;
        int ipln = 0;
        int kpbm = 0;
        int kpbln = 0;
        int iplnpkpbm = 0;
        uint64_t min_c = 0;
        __m256i a_v;
        uint64_t a = 0, c = 0, apb=0;
        int jp = 0;
        int jpm0; int iplnpjpm0; __m256i c_v0; __m256i b_v0; __m256i apb_v0; __m256i cmp_lt0; __m256i res0; 
         
        for(int k = 0; k < L1; ++k) {
            kpbm = k + sub_base_m;
            kpbln = ((k + sub_base_l) * n);
            for(int i = 0; i < L1; i += Bi) {
                for(int j = 0; j < L1; j += Bj) {
                    for(int ip = i; ip < i + Bi; ++ip) {
                        ipln = ((ip + sub_base_l) * n);
                        iplnpkpbm = ipln + kpbm;
                        a = A[iplnpkpbm];
                        a_v = _mm256_set_epi64x(a, a, a, a);
                        jp = j;
                        jpm0 = (jp + sub_base_m + 0);

                        for(; jp <= j + Bj - 4; jp += 4) {
                            iplnpjpm0 = ipln + jpm0;
                            c_v0 = _mm256_load_si256(C + iplnpjpm0);
                            b_v0 = _mm256_load_si256(B + kpbln + jpm0);
                            apb_v0 = _mm256_and_si256(a_v, b_v0);
                            res0 = _mm256_or_si256(c_v0, apb_v0);
                            _mm256_store_si256(C + iplnpjpm0, res0);
                            jpm0 += 4;
                        }
    
                        for(; jp < j + Bj; ++jp) {
                            jpm0 = (jp + sub_base_m);
                            iplnpjpm0 = ipln + jpm0;
                            c = C[iplnpjpm0];
                            apb = a & B[kpbln + jpm0];
                            min_c = c | apb;
                            C[iplnpjpm0] = min_c;
                        }
                    }
                }
            }
        }

        for(int j = 0; j < mm; ++j) {
            if(j != k) {
                //fwi_phase2_min_plus(A, B, C, n, k, j, L1, Bi, Bj);
                int l2 = k;
                int m2 = j;
                int sub_base_l = l2 * L1;
                int sub_base_m = m2 * L1;
                int ipln = 0;
                int kl = 0;
                int kln = 0;
                int iplnkl = 0;
                uint64_t apb = 0;
                uint64_t min_c = 0;
                uint64_t c = 0;
                uint64_t a = 0;

                __m256i a_v;
                int jp = 0;
                int jpm0; int iplnjpm0; __m256i c_v0; __m256i b_v0; __m256i apb_v0; __m256i cmp_lt0; __m256i res0; uint64_t *Bkj0; 

                for(int k = 0; k < L1; ++k) {
                    kl = (k + sub_base_l);
                    kln = kl * n;
                    for(int i = 0; i < L1; i += Bi) {
                        for(int j = 0; j < L1; j += Bj) {
                            for(int ip = i; ip < i + Bi; ++ip) {
                                ipln = ((ip + sub_base_l) * n);
                                iplnkl = ipln + kl;
                                a = A[iplnkl];
                                a_v = _mm256_set_epi64x(a, a, a, a);
                                jp = j;
                                jpm0 = (jp + sub_base_m + 0); iplnjpm0= ipln + jpm0; Bkj0=B + kln + jpm0;

                                for(; jp <= j + Bj - 4; jp += 4) {
                                    b_v0 = _mm256_load_si256(Bkj0);
                                    c_v0 = _mm256_load_si256(C + iplnjpm0);
                                    apb_v0 = _mm256_and_si256(a_v, b_v0);
                                    res0 = _mm256_or_si256(c_v0, apb_v0);
                                    _mm256_store_si256(C + iplnjpm0, res0);
                                    iplnjpm0 += 4;
                                    Bkj0 += 4;
                                }

                                for(; jp < j + Bj; ++jp) {
                                    jpm0 = (jp + sub_base_m);
                                    iplnjpm0 = ipln + jpm0; 
                                    apb = a + B[kln + jpm0];
                                    c = C[iplnjpm0];
                                    min_c = min(c, apb);
                                    C[iplnjpm0] = min_c;
                                }
                            }
                        }
                    }
                }
            }
        }

        for(int i = 0; i < mm; ++i) {
            if(i != k) {
                int l3 = k;
                int m3 = i;
                int sub_base_l = l3 * L1;
                int sub_base_m = m3 * L1;
                int ipmn = 0;
                int kl = 0;
                int kln = 0;
                int ipmnkl = 0;
                int jp = 0;
                uint64_t apb = 0;
                uint64_t c = 0;
                uint64_t min_c = 0;
                uint64_t a = 0;
                __m256i a_v;
                int jpl0; int ipmnjpl0; int klnjpl0; __m256i c_v0; __m256i b_v0; __m256i apb_v0; __m256i cmp_lt0; __m256i res0; 

                for(int k = 0; k < L1; ++k) {
                    kl = (k + sub_base_l);
                    kln = kl * n;
                    for(int i = 0; i < L1; i += Bi) {
                        for(int j = 0; j < L1; j += Bj) {
                            for(int ip = i; ip < i + Bi; ++ip) {
                                ipmn = ((ip + sub_base_m) * n);
                                ipmnkl = ipmn + kl; 
                                a = A[ipmnkl];
                                a_v = _mm256_set_epi64x(a, a, a, a);
                                jp = j;
                                jpl0 = (jp + sub_base_l + 0);

                                for(; jp <= j + Bj - 4; jp += 4) {
                                    ipmnjpl0 = ipmn + jpl0;
                                    klnjpl0 = kln + jpl0;
                                    b_v0 = _mm256_load_si256(B + klnjpl0);
                                    c_v0 = _mm256_load_si256(C + ipmnjpl0);
                                    apb_v0 = _mm256_and_si256(a_v, b_v0);
                                    res0 = _mm256_or_si256(c_v0, apb_v0);
                                    _mm256_store_si256(C + ipmnjpl0, res0);
                                    jpl0 += 4;
                                }

                                for(; jp < j + Bj; ++jp) {
                                    jpl0 = (jp + sub_base_l);
                                    ipmnjpl0 = ipmn + jpl0;
                                    klnjpl0 = kln + jpl0;
                                    apb = a & B[klnjpl0];
                                    c = C[ipmnjpl0];
                                    min_c = c | apb;
                                    C[ipmnjpl0] = min_c;
                                }
                            }
                        }
                    }
                }
            }
        }

        for(int i = 0; i < mm; ++i) {
            if(i != k){
                for(int j = 0; j < mm; ++j) {
                    if(j != k) {
                        int l4 = k, m4 = i, o4 = j;
                        int sub_base_l = l4 * L1;
                        int sub_base_m = m4 * L1;
                        int sub_base_o = o4 * L1;

                        int iBi = 0, jBj4 = 0, jBj = 0, kBk = 0;

                        int kpsubln = 0, ipsubmn = 0;
                        int kpsubl = 0;
                        
                        uint64_t *C_i, *B_k;
                        int jp0; int jpsubo0; __m256i c_v0; __m256i b_v0; __m256i apb_v0; __m256i cmp_lt0; __m256i res0; uint64_t *B_k_j0; uint64_t *C_i_j0; 

                        for(int i = 0; i < L1; i += Bi) {
                            iBi = i + Bi;
                            for(int j = 0; j < L1; j += Bj) {
                                jBj = j + Bj;
                                jBj4 = jBj - 4;
                                for(int k = 0; k < L1; k += Bk) {
                                    kBk = k + Bk;
                                    for(int kp = k; kp < kBk; ++kp) {
                                        kpsubl = kp + sub_base_l;
                                        kpsubln = (kpsubl) * n;
                                        B_k = B + kpsubln;
                                        for(int ip = i; ip < iBi; ++ip) {
                                            ipsubmn = (ip + sub_base_m) * n;
                                            C_i = C + ipsubmn; 
                                            jp = 0;
                                            uint64_t aa = A[ipsubmn + kpsubl];
                                            a_v = _mm256_set_epi64x(aa, aa, aa, aa);
                                            jpsubo0 = (sub_base_o + 0);

                                            for(; jp <= jBj4; jp += 4) {
                                                B_k_j0 = B_k + jpsubo0; C_i_j0 = C_i + jpsubo0;
                                                b_v0 = _mm256_load_si256(B_k_j0);
                                                c_v0 = _mm256_load_si256(C_i_j0);
                                                apb_v0 = _mm256_and_si256(a_v, b_v0);
                                                res0 = _mm256_or_si256(c_v0, apb_v0);
                                                _mm256_store_si256(C_i_j0, res0);
                                                jpsubo0+=4;
                                            }

                                            for(; jp < j + Bj; ++jp) {
                                                C[ipsubmn + (jp + sub_base_o)] =
                                                    C[ipsubmn + (jp + sub_base_o)] | (A[ipsubmn + (kp + sub_base_l)] & B[((kp + sub_base_l) * n) + (jp + sub_base_o)]);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void init_matrices(uint64_t *C1, uint64_t *C2, int n) {
    uint64_t x;
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
double rdtsc_generalized(uint64_t *A, uint64_t *B, uint64_t *C, int n,
        void (*compute)(uint64_t*, uint64_t*, uint64_t*, int)) {

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
            compute(A, B, C, n);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        compute(A, B, C, n);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
double rdtsc_tiled(uint64_t *A, uint64_t *B, uint64_t *C, int n, int L1, int Bi, int Bj, int Bk, 
        void (*compute)(uint64_t*, uint64_t*, uint64_t*, int, int, int, int, int)) {

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
            compute(A, B, C, L1, n, Bi, Bj, Bk);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for (i = 0; i < num_runs; ++i) {
        compute(A, B, C, L1, n, Bi, Bj, Bk);
    }

    cycles = stop_tsc(start)/num_runs;
    return (double) cycles;
}
#endif

#define epsilon 0.00000001
double benchmark_tiled_timed(int n, void (*baseline)(uint64_t*, uint64_t*, uint64_t*, int), 
    void (*compute)(uint64_t*, uint64_t*, uint64_t*, int, int, int, int, int),
    int L1, int Bi, int Bj, int Bk
) {

    uint64_t* C_base = (uint64_t*)malloc(n*n*sizeof(uint64_t));
    uint64_t* C_opt = (uint64_t*)aligned_alloc(32, n*n*sizeof(uint64_t));

    init_matrices(C_base, C_opt, n);

    double base = rdtsc_generalized(C_base, C_base, C_base, n, baseline);

    printf(" %f \n ", base);

    double time = rdtsc_tiled(C_opt, C_opt, C_opt, n, L1, Bi, Bj, Bk, compute);

    printf(" %f \n ", time);

    // Compare both 
    for(int i = 0; i < n*n; ++i) {
        assert(abs(C_opt[i] - C_base[i]) <= epsilon);
    }

    free(C_base);
    free(C_opt);

    return time;
}

int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    int L1 = atoi(argv[2]);
    int Bi,Bj,Bk;
    Bi = Bj = Bk = atoi(argv[3]);

#ifdef __x86_64__
    double r1 = benchmark_tiled_timed(n, fw_abc_or_and, opt_tiled, L1, Bi, Bj, Bk);
#endif

    return 0;
}
    