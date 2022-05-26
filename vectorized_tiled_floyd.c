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

void fw_abc_min_plus(double* A, double* B, double* C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j] = min(C[i*n + j], A[i*n + k] + B[k*n + j]);
            }
        }
    }
}

void fw_abc_max_min(double* A, double* B, double* C, int n) {
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                C[i*n + j] = max(C[i*n + j], min(A[i*n + k], B[k*n + j]));
            }
        }
    }
}

/**
 * @brief Tiled FW implementation for MIN_PLUS
 * @param A, the first operand matrix
 * @param B, the second operand matrix
 * @param C, the result matrix
 * @param L1, the tile size L1xL1 of the matrix
 * @param n, the size of the NxN matrices
 * @param Bi, tilling factor over i
 * @param Bj, tilling factor over j
 * @param Bk, tilling factor over k
 */
// NOTE All functions are inlined
void tiled_fw_min_plus(double* A, double* B, double* C, int L1, int n, int Bi, int Bj, int Bk) {
    int mm = n / L1;
    // printf("L1 : %d, Bi : %d, Bj : %d, Bk : %d, m : %d\n", L1, Bi, Bj, Bk, m);
    for(int k = 0; k < mm; ++k) {
        //Tilling phase 1 (update C_kk)
        int l1 = k;
        int m1 = k;
        int sub_base_l = l1 * L1;
        int sub_base_m = m1 * L1;
        int ipln = 0;
        int jpm = 0;
        int iplnpjpm = 0;
        int kpbm = 0;
        int kpbln = 0;
        int iplnpkpbm = 0;
        double min_c = 0.0;
        __m256d a_v,c_v,b_v,apb_v,cmp_lt,res;
        double a = 0.0, c = 0.0, apb=0.0;
        int jp = 0;
        for (int k = 0; k < L1; ++k) {
            kpbm = k + sub_base_m;
            kpbln = ((k + sub_base_l) * n);
            for (int i = 0; i < L1; i += Bi) {
                for (int j = 0; j < L1; j += Bj) {
                    for(int ip = i; ip < i + Bi; ++ip) {
                        ipln = ((ip + sub_base_l) * n);
                        iplnpkpbm = ipln + kpbm;
                        a = A[iplnpkpbm];
                        a_v = _mm256_set1_pd(a);
                        jp = j;
                        for(; jp <= j + Bj - 4; jp+=4) {
                            jpm = (jp + sub_base_m);
                            iplnpjpm = ipln + jpm;
                            c_v = _mm256_loadu_pd(C + iplnpjpm);
                            b_v = _mm256_loadu_pd(B + kpbln + jpm);
                            apb_v = _mm256_add_pd(a_v, b_v);
                            cmp_lt = _mm256_cmp_pd(c_v, apb_v, _CMP_LT_OQ);
                            res = _mm256_blendv_pd(apb_v, c_v, cmp_lt);
                            _mm256_storeu_pd(C + iplnpjpm, res);
                        }
                        for(; jp < j + Bj; ++jp) {
                            jpm = (jp + sub_base_m);
                            iplnpjpm = ipln + jpm;
                            c = C[iplnpjpm];
                            apb = a + B[kpbln + jpm];
                            min_c = min(c, apb);
                            C[iplnpjpm] = min_c;
                        }
                    }
                }
            }
        }

        //Tilling phase 2 (Update all tiles in same row as C_kk)
        for(int j = 0; j < mm; ++j) {
            if(j != k) {
                //fwi_phase2_min_plus(A, B, C, n, k, j, L1, Bi, Bj);
                int l2 = k;
                int m2 = j;
                int sub_base_l = l2 * L1;
                int sub_base_m = m2 * L1;
                int ipln = 0;
                int jpm = 0;
                int kl = 0;
                int kln = 0;
                int iplnkl = 0;
                int iplnjpm = 0;
                double apb = 0.0;
                double min_c = 0.0;
                double c = 0.0;
                double a = 0.0;
                __m256d a_v,c_v,b_v,apb_v,cmp_lt,res;
                int jp = 0;
                for (int k = 0; k < L1; ++k) {
                    kl = (k + sub_base_l);
                    kln = kl * n;
                    for (int i = 0; i < L1; i += Bi) {
                        for (int j = 0; j < L1; j += Bj) {
                            for(int ip = i; ip < i + Bi; ++ip) {
                                ipln = ((ip + sub_base_l) * n);
                                iplnkl = ipln + kl;
                                a = A[iplnkl];
                                a_v = _mm256_set1_pd(a);
                                jp = j;
                                for(; jp <= j + Bj - 4; jp+=4) {
                                    jpm = (jp + sub_base_m);
                                    iplnjpm = ipln + jpm;
                                    b_v = _mm256_loadu_pd(B + kln + jpm); 
                                    apb_v = _mm256_add_pd(a_v, b_v);
                                    c_v = _mm256_loadu_pd(C + iplnjpm);
                                    cmp_lt = _mm256_cmp_pd(c_v, apb_v, _CMP_LT_OQ);
                                    res = _mm256_blendv_pd(apb_v, c_v, cmp_lt);
                                    _mm256_storeu_pd(C + iplnjpm, res);
                                }
                                for(; jp < j + Bj; ++jp) {
                                    jpm = (jp + sub_base_m);
                                    iplnjpm = ipln + jpm; 
                                    apb = a + B[kln + jpm];
                                    c = C[iplnjpm];
                                    min_c = min(c, apb);
                                    C[iplnjpm] = min_c;
                                }
                            }
                        }
                    }
                }
            }
        }

        //Tilling phase 3 (Update all tiles in same column as C_kk)
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
                int jpl = 0;
                int ipmnjpl = 0;
                int klnjpl = 0;
                double apb = 0.0;
                double c = 0.0;
                double min_c = 0.0;
                double a = 0.0;
                int jp = 0;
                __m256d a_v,c_v,b_v,apb_v,cmp_lt,res;
                for (int k = 0; k < L1; ++k) {
                    kl = (k + sub_base_l);
                    kln = kl * n;
                    for (int i = 0; i < L1; i += Bi) {
                        for (int j = 0; j < L1; j += Bj) {
                            for(int ip = i; ip < i + Bi; ++ip) {
                                ipmn = ((ip + sub_base_m) * n);
                                ipmnkl = ipmn + kl; 
                                a = A[ipmnkl];
                                a_v = _mm256_set1_pd(a);
                                jp = j;
                                for(; jp < j + Bj - 4; jp+=4) {
                                    jpl = (jp + sub_base_l);
                                    ipmnjpl = ipmn + jpl;
                                    klnjpl = kln + jpl;
                                    b_v = _mm256_loadu_pd(B + klnjpl);
                                    apb_v = _mm256_add_pd(a_v, b_v);
                                    c_v = _mm256_loadu_pd(C + ipmnjpl);
                                    cmp_lt = _mm256_cmp_pd(c_v, apb_v, _CMP_LT_OQ);
                                    res = _mm256_blendv_pd(apb_v, c_v, cmp_lt);
                                    _mm256_storeu_pd(C + ipmnjpl, res);
                                }
                                for(; jp < j + Bj; ++jp) {
                                    jpl = (jp + sub_base_l);
                                    ipmnjpl = ipmn + jpl;
                                    klnjpl = kln + jpl;
                                    apb = a + B[klnjpl];
                                    c = C[ipmnjpl];
                                    min_c = min(c, apb);
                                    C[ipmnjpl] = min_c;
                                }
                            }
                        }
                    }
                }
            }
        }

        //Tilling phase 4 (Update all remaining tiles)
        for(int i = 0; i < mm; ++i) {
            if(i != k){
                for(int j = 0; j < mm; ++j) {
                    if(j != k) {
                        int l4 = k; 
                        int m4 = i;
                        int o4 = j;
                        int sub_base_l = l4 * L1;
                        int sub_base_m = m4 * L1;
                        int sub_base_o = o4 * L1;
                        int ipmn = 0;
                        int jpo = 0;
                        int ipmnpjpo = 0;
                        int kpl = 0;
                        int kplnpjpo = 0;
                        double apb = 0.0;
                        double c = 0.0;
                        double min_c = 0.0;
                        // int kp = 0;
                        // __m256d a_v,c_v,b_v,apb_v,cmp_lt,res;
                        for (int i = 0; i < L1; i += Bi) {
                            for (int j = 0; j < L1; j += Bj) {
                                for (int k = 0; k < L1; k += Bk) {
                                    for(int ip = i; ip < i + Bi; ++ip) {
                                        ipmn = (ip + sub_base_m) * n;
                                        for(int jp = j; jp < j + Bj; ++jp) {
                                            jpo = (jp + sub_base_o);
                                            ipmnpjpo = ipmn + jpo;
                                            c = C[ipmnpjpo];
                                            min_c = c;
                                            for(int kp = k; kp < k + Bk; ++kp) {
                                                kpl = kp + sub_base_l;
                                                kplnpjpo = (kpl * n) + jpo;
                                                apb = A[ipmn + kpl] + B[kplnpjpo];
                                                min_c = min(min_c, apb);
                                            }
                                            C[ipmnpjpo] = min_c;
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

void init_matrices(double *C1, double *C2, int n) {
    double x;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            x = ((double )rand() + 1);
            C1[i * n + j] = x;
            C2[i * n + j] = x;
        }
    }   
}

/* 
 * Timing function based on the TimeStep Counter of the CPU.
 */
#ifdef __x86_64__
double rdtsc_generalized(double *A, double *B, double *C, int n,
        void (*compute)(double*, double*, double*, int)) {

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
double rdtsc_tiled(double *A, double *B, double *C, int n, int L1, int Bi, int Bj, int Bk, 
        void (*compute)(double*, double*, double*, int, int, int, int, int)) {

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
double benchmark_tiled_timed(int n, void (*baseline)(double*, double*, double*, int), 
    void (*compute)(double*, double*, double*, int, int, int, int, int),
    int L1, int Bi, int Bj, int Bk
) {

    double *C_base = (double *)malloc(n*n*sizeof(double));
    double *C_opt = (double *)aligned_alloc(32, n*n*sizeof(double));

    init_matrices(C_base, C_opt, n);

    printf("RECEIVED L1 : %d\n", L1);

    // Run baseline function on C
    //baseline(C_base, C_base, C_base, n);

    // Run optimized function on C
    //compute(C_opt, C_opt, C_opt, L1, n, Bi, Bj, Bk);

    double base = rdtsc_generalized(C_base, C_base, C_base, n, baseline);

    printf("Time taken to run generalized : %f\n", base);

    double time = rdtsc_tiled(C_opt, C_opt, C_opt, n, L1, Bi, Bj, Bk, compute);

    printf("Time taken to run tiled : %f\n", time);

    // Compare both 
    for(int i = 0; i < n*n; ++i) {
        assert(abs(C_opt[i] - C_base[i]) <= epsilon);
    }

    printf("Comparison passed!\n");

    free(C_base);
    free(C_opt);

    return time;
}

int main(int argc, char **argv) {
    if (argc!=4) {printf("usage: FW <n> <L1> <B> (fw = 0,1,2 = (min,plus), (or,and), (max, min)), (tiled = 0,1)\n"); return -1;}
    int n = atoi(argv[1]);
    int L1 = atoi(argv[2]);
    int Bi,Bj,Bk;
    Bi = Bj = Bk = atoi(argv[3]);
    printf("n=%d \n",n);

#ifdef __x86_64__
    double r = benchmark_tiled_timed(n, fw_abc_min_plus, tiled_fw_min_plus, L1, Bi, Bj, Bk);
    printf(" FW : RDTSC instruction:\n %lf cycles measured\n\n", r);
#endif

    return 0;
}