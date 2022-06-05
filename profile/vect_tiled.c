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

#define NUM_RUNS 1
#define CYCLES_REQUIRED 1e8
#define FREQUENCY 2.7e9
#define CALIBRATE
#define ZERO_PROBABILITY 10 //1/ZERO_PROBABILITY is the probability of an entry in the bit matrix being zero
#define EPS  0.000001

    void opt_tiled_fw_min_plus(double* A, double* B, double* C, int L1, int n, int Bi, int Bj, int Bk) {
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
            int jpm0;int iplnpjpm0;__m256d c_v0;__m256d b_v0;__m256d apb_v0;__m256d cmp_lt0;__m256d res0;
int jpm1;int iplnpjpm1;__m256d c_v1;__m256d b_v1;__m256d apb_v1;__m256d cmp_lt1;__m256d res1;
         
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
                        jpm0 = (jp + sub_base_m + 0);
jpm1 = (jp + sub_base_m + 4);
for(; jp <= j + Bj - 4; jp += 8) {
iplnpjpm0 = ipln + jpm0;
c_v0 = _mm256_load_pd(C + iplnpjpm0);
b_v0 = _mm256_load_pd(B + kpbln + jpm0);
apb_v0 = _mm256_add_pd(a_v, b_v0);
res0 = _mm256_min_pd(c_v0, apb_v0);
_mm256_store_pd(C + iplnpjpm0, res0);
jpm0 += 8;
iplnpjpm1 = ipln + jpm1;
c_v1 = _mm256_load_pd(C + iplnpjpm1);
b_v1 = _mm256_load_pd(B + kpbln + jpm1);
apb_v1 = _mm256_add_pd(a_v, b_v1);
res1 = _mm256_min_pd(c_v1, apb_v1);
_mm256_store_pd(C + iplnpjpm1, res1);
jpm1 += 8;
}    
                            for(; jp < j + Bj; ++jp) {
                                jpm0 = (jp + sub_base_m);
                                iplnpjpm0 = ipln + jpm0;
                                c = C[iplnpjpm0];
                                apb = a + B[kpbln + jpm0];
                                min_c = min(c, apb);
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
                int jpm0;int iplnjpm0;__m256d c_v0;__m256d b_v0;__m256d apb_v0;__m256d cmp_lt0;__m256d res0;double *Bkj0;
int jpm1;int iplnjpm1;__m256d c_v1;__m256d b_v1;__m256d apb_v1;__m256d cmp_lt1;__m256d res1;double *Bkj1;

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
                                jpm0 = (jp + sub_base_m + 0); iplnjpm0= ipln + jpm0; Bkj0=B + kln + jpm0;
jpm1 = (jp + sub_base_m + 4); iplnjpm1= ipln + jpm1; Bkj1=B + kln + jpm1;
for(; jp <= j + Bj - 8; jp += 8) {
b_v0 = _mm256_load_pd(Bkj0)
;apb_v0 = _mm256_add_pd(a_v, b_v0)
; c_v0 = _mm256_load_pd(C + iplnjpm0)
;res0 = _mm256_min_pd(c_v0, apb_v0);
_mm256_store_pd(C + iplnjpm0, res0)
; iplnjpm0 += 8
; Bkj0 += 8
;b_v1 = _mm256_load_pd(Bkj1)
;apb_v1 = _mm256_add_pd(a_v, b_v1)
; c_v1 = _mm256_load_pd(C + iplnjpm1)
;res1 = _mm256_min_pd(c_v1, apb_v1);
_mm256_store_pd(C + iplnjpm1, res1)
; iplnjpm1 += 8
; Bkj1 += 8
;}

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
                    double apb = 0.0;
                    double c = 0.0;
                    double min_c = 0.0;
                    double a = 0.0;
                    __m256d a_v;

                    int jpl0;int ipmnjpl0;int klnjpl0;__m256d c_v0;__m256d b_v0;__m256d apb_v0;__m256d cmp_lt0;__m256d res0;
int jpl1;int ipmnjpl1;int klnjpl1;__m256d c_v1;__m256d b_v1;__m256d apb_v1;__m256d cmp_lt1;__m256d res1;

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
                                    jpl0 = (jp + sub_base_l + 0);
jpl1 = (jp + sub_base_l + 4);
for(; jp <= j + Bj - 8; jp += 8) {ipmnjpl0 = ipmn + jpl0;
klnjpl0 = kln + jpl0;
b_v0 = _mm256_load_pd(B + klnjpl0);
apb_v0 = _mm256_add_pd(a_v, b_v0);
c_v0 = _mm256_load_pd(C + ipmnjpl0);
res0 = _mm256_min_pd(c_v0, apb_v0);
_mm256_store_pd(C + ipmnjpl0, res0);
jpl0 += 8;
ipmnjpl1 = ipmn + jpl1;
klnjpl1 = kln + jpl1;
b_v1 = _mm256_load_pd(B + klnjpl1);
apb_v1 = _mm256_add_pd(a_v, b_v1);
c_v1 = _mm256_load_pd(C + ipmnjpl1);
res1 = _mm256_min_pd(c_v1, apb_v1);
_mm256_store_pd(C + ipmnjpl1, res1);
jpl1 += 8;
}
                                    for(; jp < j + Bj; ++jp) {
                                        jpl0 = (jp + sub_base_l);
                                        ipmnjpl0 = ipmn + jpl0;
                                        klnjpl0 = kln + jpl0;
                                        apb = a + B[klnjpl0];
                                        c = C[ipmnjpl0];
                                        min_c = min(c, apb);
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
                        int jp0;int jpsubo0;__m256d c_v0;__m256d b_v0;__m256d apb_v0;__m256d cmp_lt0;__m256d res0;double *B_k_j0;double *C_i_j0;
int jp1;int jpsubo1;__m256d c_v1;__m256d b_v1;__m256d apb_v1;__m256d cmp_lt1;__m256d res1;double *B_k_j1;double *C_i_j1;

                            for (int i = 0; i < L1; i += Bi) {
                                iBi = i + Bi;
                                for (int j = 0; j < L1; j += Bj) {
                                    jBj = j + Bj;
                                    jBj4 = jBj - 8;
                                    for (int k = 0; k < L1; k += Bk) {
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
for(; jp <= jBj4; jp += 8) {B_k_j0 = B_k + jpsubo0; C_i_j0 = C_i + jpsubo0;
b_v0 = _mm256_load_pd(B_k_j0);
c_v0 = _mm256_load_pd(C_i_j0);
apb_v0 = _mm256_add_pd(a_v, b_v0);
res0 = _mm256_min_pd(c_v0, apb_v0);
_mm256_store_pd(C_i_j0, res0);
jpsubo0+=8;
B_k_j1 = B_k + jpsubo1; C_i_j1 = C_i + jpsubo1;
b_v1 = _mm256_load_pd(B_k_j1);
c_v1 = _mm256_load_pd(C_i_j1);
apb_v1 = _mm256_add_pd(a_v, b_v1);
res1 = _mm256_min_pd(c_v1, apb_v1);
_mm256_store_pd(C_i_j1, res1);
jpsubo1+=8;
}
                                                for(; jp < j + Bj; ++jp) {
                                                    C[ipsubmn + (jp + sub_base_o)] = min(
                                                    C[ipsubmn + (jp + sub_base_o)], 
                                                    A[ipsubmn + (kp + sub_base_l)] + 
                                                        B[((kp + sub_base_l) * n) + (jp + sub_base_o)]
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

int main(int argc, char **argv) {
    if(argc != 4)
        printf("Wrong\n");

    int n = atoi(argv[1]);
    int L1 = atoi(argv[2]);
    int Bi, Bj, Bk;
    Bi = Bj = Bk = atoi(argv[3]);

    double *C = (double *)aligned_alloc(32, n*n*sizeof(double));

    /*for(int i = 0; i < n*n; i++)
        C[i] = i%3 == 0 ? 2.7 : (i%3 == 1 ? 3.2 : 0.3);*/

    printf("%d, %d, %d\n", n, L1, Bi);
    opt_tiled_fw_min_plus(C, C, C, L1, n , Bi, Bj, Bk);

    return 0;
}