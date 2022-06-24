
#ifdef linux
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))
#endif

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <immintrin.h>
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
    for(size_t k = 0; k < n; k++) {
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < n; j++) {
                C[i*n + j] = min(C[i*n + j], A[i*n + k] + B[k*n + j]);
            }
        }
    }
}

void fw_abc_max_min(double* A, double* B, double* C, int n) {
    for(size_t k = 0; k < n; k++) {
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < n; j++) {
                C[i*n + j] = max(C[i*n + j], min(A[i*n + k], B[k*n + j]));
            }
        }
    }
}

void opt_tiled(double* A, double* B, double* C, int L1, int n, int Bi, int Bj, int Bk) {
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
        double min_c = 0.0;
        __m256d a_v;
        double a = 0.0, c = 0.0, apb=0.0;
        int jp = 0;
        int jpm0; int iplnpjpm0; __m256d c_v0; __m256d b_v0; __m256d apb_v0; __m256d cmp_lt0; __m256d res0; 
        int jpm1; int iplnpjpm1; __m256d c_v1; __m256d b_v1; __m256d apb_v1; __m256d cmp_lt1; __m256d res1; 
        int jpm2; int iplnpjpm2; __m256d c_v2; __m256d b_v2; __m256d apb_v2; __m256d cmp_lt2; __m256d res2; 
         
        for(int k = 0; k < L1; ++k) {
            kpbm = k + sub_base_m;
            kpbln = ((k + sub_base_l) * n);
            for(int i = 0; i < L1; i += Bi) {
                for(int j = 0; j < L1; j += Bj) {
                    for(int ip = i; ip < i + Bi; ++ip) {
                        ipln = ((ip + sub_base_l) * n);
                        iplnpkpbm = ipln + kpbm;
                        a = A[iplnpkpbm];
                        a_v = _mm256_set1_pd(a);
                        jp = j;
                        jpm0 = (jp + sub_base_m + 0);
                        jpm1 = (jp + sub_base_m + 4);
                        jpm2 = (jp + sub_base_m + 8);

                        for(; jp <= j + Bj - 4; jp += 12) {
                            iplnpjpm0 = ipln + jpm0;
                            iplnpjpm1 = ipln + jpm1;
                            iplnpjpm2 = ipln + jpm2;
                            c_v0 = _mm256_load_pd(C + iplnpjpm0);
                            c_v1 = _mm256_load_pd(C + iplnpjpm1);
                            c_v2 = _mm256_load_pd(C + iplnpjpm2);
                            b_v0 = _mm256_load_pd(B + kpbln + jpm0);
                            b_v1 = _mm256_load_pd(B + kpbln + jpm1);
                            b_v2 = _mm256_load_pd(B + kpbln + jpm2);
                            apb_v0 = _mm256_min_pd(a_v, b_v0);
                            apb_v1 = _mm256_min_pd(a_v, b_v1);
                            apb_v2 = _mm256_min_pd(a_v, b_v2);
                            res0 = _mm256_max_pd(c_v0, apb_v0);
                            res1 = _mm256_max_pd(c_v1, apb_v1);
                            res2 = _mm256_max_pd(c_v2, apb_v2);
                            _mm256_store_pd(C + iplnpjpm0, res0);
                            _mm256_store_pd(C + iplnpjpm1, res1);
                            _mm256_store_pd(C + iplnpjpm2, res2);
                            jpm0 += 12;
                            jpm1 += 12;
                            jpm2 += 12;
                        }
    
                        for(; jp < j + Bj; ++jp) {
                            jpm0 = (jp + sub_base_m);
                            iplnpjpm0 = ipln + jpm0;
                            c = C[iplnpjpm0];
                            apb = min(a, B[kpbln + jpm0]);
                            min_c = max(c, apb);
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
                double apb = 0.0;
                double min_c = 0.0;
                double c = 0.0;
                double a = 0.0;

                __m256d a_v;
                int jp = 0;
                int jpm0; int iplnjpm0; __m256d c_v0; __m256d b_v0; __m256d apb_v0; __m256d cmp_lt0; __m256d res0; double *Bkj0; 
                int jpm1; int iplnjpm1; __m256d c_v1; __m256d b_v1; __m256d apb_v1; __m256d cmp_lt1; __m256d res1; double *Bkj1; 
                int jpm2; int iplnjpm2; __m256d c_v2; __m256d b_v2; __m256d apb_v2; __m256d cmp_lt2; __m256d res2; double *Bkj2; 

                for(int k = 0; k < L1; ++k) {
                    kl = (k + sub_base_l);
                    kln = kl * n;
                    for(int i = 0; i < L1; i += Bi) {
                        for(int j = 0; j < L1; j += Bj) {
                            for(int ip = i; ip < i + Bi; ++ip) {
                                ipln = ((ip + sub_base_l) * n);
                                iplnkl = ipln + kl;
                                a = A[iplnkl];
                                a_v = _mm256_set1_pd(a);
                                jp = j;
                                jpm0 = (jp + sub_base_m + 0); iplnjpm0= ipln + jpm0; Bkj0=B + kln + jpm0;
                                jpm1 = (jp + sub_base_m + 4); iplnjpm1= ipln + jpm1; Bkj1=B + kln + jpm1;
                                jpm2 = (jp + sub_base_m + 8); iplnjpm2= ipln + jpm2; Bkj2=B + kln + jpm2;

                                for(; jp <= j + Bj - 12; jp += 12) {
                                    b_v0 = _mm256_load_pd(Bkj0);
                                    b_v1 = _mm256_load_pd(Bkj1);
                                    b_v2 = _mm256_load_pd(Bkj2);
                                    c_v0 = _mm256_load_pd(C + iplnjpm0);
                                    c_v1 = _mm256_load_pd(C + iplnjpm1);
                                    c_v2 = _mm256_load_pd(C + iplnjpm2);
                                    apb_v0 = _mm256_min_pd(a_v, b_v0);
                                    apb_v1 = _mm256_min_pd(a_v, b_v1);
                                    apb_v2 = _mm256_min_pd(a_v, b_v2);
                                    res0 = _mm256_max_pd(c_v0, apb_v0);
                                    res1 = _mm256_max_pd(c_v1, apb_v1);
                                    res2 = _mm256_max_pd(c_v2, apb_v2);
                                    _mm256_store_pd(C + iplnjpm0, res0);
                                    _mm256_store_pd(C + iplnjpm1, res1);
                                    _mm256_store_pd(C + iplnjpm2, res2);
                                    iplnjpm0 += 12;
                                    iplnjpm1 += 12;
                                    iplnjpm2 += 12;
                                    Bkj0 += 12;
                                    Bkj1 += 12;
                                    Bkj2 += 12;
                                }

                                for(; jp < j + Bj; ++jp) {
                                    jpm0 = (jp + sub_base_m);
                                    iplnjpm0 = ipln + jpm0; 
                                    apb = min(a, B[kln + jpm0]);
                                    c = C[iplnjpm0];
                                    min_c = max(c, apb);
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
                double apb = 0.0;
                double c = 0.0;
                double min_c = 0.0;
                double a = 0.0;
                __m256d a_v;
                int jpl0; int ipmnjpl0; int klnjpl0; __m256d c_v0; __m256d b_v0; __m256d apb_v0; __m256d cmp_lt0; __m256d res0; 
                int jpl1; int ipmnjpl1; int klnjpl1; __m256d c_v1; __m256d b_v1; __m256d apb_v1; __m256d cmp_lt1; __m256d res1; 
                int jpl2; int ipmnjpl2; int klnjpl2; __m256d c_v2; __m256d b_v2; __m256d apb_v2; __m256d cmp_lt2; __m256d res2; 

                for(int k = 0; k < L1; ++k) {
                    kl = (k + sub_base_l);
                    kln = kl * n;
                    for(int i = 0; i < L1; i += Bi) {
                        for(int j = 0; j < L1; j += Bj) {
                            for(int ip = i; ip < i + Bi; ++ip) {
                                ipmn = ((ip + sub_base_m) * n);
                                ipmnkl = ipmn + kl; 
                                a = A[ipmnkl];
                                a_v = _mm256_set1_pd(a);
                                jp = j;
                                jpl0 = (jp + sub_base_l + 0);
                                jpl1 = (jp + sub_base_l + 4);
                                jpl2 = (jp + sub_base_l + 8);

                                for(; jp <= j + Bj - 12; jp += 12) {
                                    ipmnjpl0 = ipmn + jpl0;
                                    ipmnjpl1 = ipmn + jpl1;
                                    ipmnjpl2 = ipmn + jpl2;
                                    klnjpl0 = kln + jpl0;
                                    klnjpl1 = kln + jpl1;
                                    klnjpl2 = kln + jpl2;
                                    b_v0 = _mm256_load_pd(B + klnjpl0);
                                    b_v1 = _mm256_load_pd(B + klnjpl1);
                                    b_v2 = _mm256_load_pd(B + klnjpl2);
                                    c_v0 = _mm256_load_pd(C + ipmnjpl0);
                                    c_v1 = _mm256_load_pd(C + ipmnjpl1);
                                    c_v2 = _mm256_load_pd(C + ipmnjpl2);
                                    apb_v0 = _mm256_min_pd(a_v, b_v0);
                                    apb_v1 = _mm256_min_pd(a_v, b_v1);
                                    apb_v2 = _mm256_min_pd(a_v, b_v2);
                                    res0 = _mm256_max_pd(c_v0, apb_v0);
                                    res1 = _mm256_max_pd(c_v1, apb_v1);
                                    res2 = _mm256_max_pd(c_v2, apb_v2);
                                    _mm256_store_pd(C + ipmnjpl0, res0);
                                    _mm256_store_pd(C + ipmnjpl1, res1);
                                    _mm256_store_pd(C + ipmnjpl2, res2);
                                    jpl0 += 12;
                                    jpl1 += 12;
                                    jpl2 += 12;
                                }

                                for(; jp < j + Bj; ++jp) {
                                    jpl0 = (jp + sub_base_l);
                                    ipmnjpl0 = ipmn + jpl0;
                                    klnjpl0 = kln + jpl0;
                                    apb = min(a, B[klnjpl0]);
                                    c = C[ipmnjpl0];
                                    min_c = max(c, apb);
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
                        
                        double *C_i, *B_k;
                        int jp0; int jpsubo0; __m256d c_v0; __m256d b_v0; __m256d apb_v0; __m256d cmp_lt0; __m256d res0; double *B_k_j0; double *C_i_j0; 
                        int jp1; int jpsubo1; __m256d c_v1; __m256d b_v1; __m256d apb_v1; __m256d cmp_lt1; __m256d res1; double *B_k_j1; double *C_i_j1; 
                        int jp2; int jpsubo2; __m256d c_v2; __m256d b_v2; __m256d apb_v2; __m256d cmp_lt2; __m256d res2; double *B_k_j2; double *C_i_j2; 

                        for(int i = 0; i < L1; i += Bi) {
                            iBi = i + Bi;
                            for(int j = 0; j < L1; j += Bj) {
                                jBj = j + Bj;
                                jBj4 = jBj - 12;
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
                                            a_v = _mm256_set1_pd(A[ipsubmn + kpsubl]);
                                            jpsubo0 = (sub_base_o + 0);
                                            jpsubo1 = (sub_base_o + 4);
                                            jpsubo2 = (sub_base_o + 8);

                                            for(; jp <= jBj4; jp += 12) {
                                                B_k_j0 = B_k + jpsubo0; C_i_j0 = C_i + jpsubo0;
                                                B_k_j1 = B_k + jpsubo1; C_i_j1 = C_i + jpsubo1;
                                                B_k_j2 = B_k + jpsubo2; C_i_j2 = C_i + jpsubo2;
                                                b_v0 = _mm256_load_pd(B_k_j0);
                                                b_v1 = _mm256_load_pd(B_k_j1);
                                                b_v2 = _mm256_load_pd(B_k_j2);
                                                c_v0 = _mm256_load_pd(C_i_j0);
                                                c_v1 = _mm256_load_pd(C_i_j1);
                                                c_v2 = _mm256_load_pd(C_i_j2);
                                                apb_v0 = _mm256_min_pd(a_v, b_v0);
                                                apb_v1 = _mm256_min_pd(a_v, b_v1);
                                                apb_v2 = _mm256_min_pd(a_v, b_v2);
                                                res0 = _mm256_max_pd(c_v0, apb_v0);
                                                res1 = _mm256_max_pd(c_v1, apb_v1);
                                                res2 = _mm256_max_pd(c_v2, apb_v2);
                                                _mm256_store_pd(C_i_j0, res0);
                                                _mm256_store_pd(C_i_j1, res1);
                                                _mm256_store_pd(C_i_j2, res2);
                                                jpsubo0+=12;
                                                jpsubo1+=12;
                                                jpsubo2+=12;
                                            }

                                            for(; jp < j + Bj; ++jp) {
                                                C[ipsubmn + (jp + sub_base_o)] = max(
                                                    C[ipsubmn + (jp + sub_base_o)], 
                                                    min(A[ipsubmn + (kp + sub_base_l)], B[((kp + sub_base_l) * n) + (jp + sub_base_o)])
                                                );
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


void init_matrices(double *C1, double *C2, int n) {
    double x;
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < n; j++) {
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
        for(i = 0; i < num_runs; ++i) {
            compute(A, B, C, n);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for(i = 0; i < num_runs; ++i) {
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
        for(i = 0; i < num_runs; ++i) {
            compute(A, B, C, L1, n, Bi, Bj, Bk);
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
#endif

    start = start_tsc();
    for(i = 0; i < num_runs; ++i) {
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

    double base = rdtsc_generalized(C_base, C_base, C_base, n, baseline);

    double time = rdtsc_tiled(C_opt, C_opt, C_opt, n, L1, Bi, Bj, Bk, compute);

    printf("%f\n", time);
    printf("%f\n", base);

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
    double r1 = benchmark_tiled_timed(n, fw_abc_max_min, opt_tiled, L1, Bi, Bj, Bk);
#endif

    return 0;
}
    