
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
int jpm2;int iplnpjpm2;__m256d c_v2;__m256d b_v2;__m256d apb_v2;__m256d cmp_lt2;__m256d res2;
int jpm3;int iplnpjpm3;__m256d c_v3;__m256d b_v3;__m256d apb_v3;__m256d cmp_lt3;__m256d res3;
int jpm4;int iplnpjpm4;__m256d c_v4;__m256d b_v4;__m256d apb_v4;__m256d cmp_lt4;__m256d res4;
int jpm5;int iplnpjpm5;__m256d c_v5;__m256d b_v5;__m256d apb_v5;__m256d cmp_lt5;__m256d res5;
int jpm6;int iplnpjpm6;__m256d c_v6;__m256d b_v6;__m256d apb_v6;__m256d cmp_lt6;__m256d res6;
int jpm7;int iplnpjpm7;__m256d c_v7;__m256d b_v7;__m256d apb_v7;__m256d cmp_lt7;__m256d res7;
int jpm8;int iplnpjpm8;__m256d c_v8;__m256d b_v8;__m256d apb_v8;__m256d cmp_lt8;__m256d res8;
int jpm9;int iplnpjpm9;__m256d c_v9;__m256d b_v9;__m256d apb_v9;__m256d cmp_lt9;__m256d res9;
int jpm10;int iplnpjpm10;__m256d c_v10;__m256d b_v10;__m256d apb_v10;__m256d cmp_lt10;__m256d res10;
int jpm11;int iplnpjpm11;__m256d c_v11;__m256d b_v11;__m256d apb_v11;__m256d cmp_lt11;__m256d res11;
int jpm12;int iplnpjpm12;__m256d c_v12;__m256d b_v12;__m256d apb_v12;__m256d cmp_lt12;__m256d res12;
int jpm13;int iplnpjpm13;__m256d c_v13;__m256d b_v13;__m256d apb_v13;__m256d cmp_lt13;__m256d res13;
int jpm14;int iplnpjpm14;__m256d c_v14;__m256d b_v14;__m256d apb_v14;__m256d cmp_lt14;__m256d res14;
int jpm15;int iplnpjpm15;__m256d c_v15;__m256d b_v15;__m256d apb_v15;__m256d cmp_lt15;__m256d res15;
         
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
jpm2 = (jp + sub_base_m + 8);
jpm3 = (jp + sub_base_m + 12);
jpm4 = (jp + sub_base_m + 16);
jpm5 = (jp + sub_base_m + 20);
jpm6 = (jp + sub_base_m + 24);
jpm7 = (jp + sub_base_m + 28);
jpm8 = (jp + sub_base_m + 32);
jpm9 = (jp + sub_base_m + 36);
jpm10 = (jp + sub_base_m + 40);
jpm11 = (jp + sub_base_m + 44);
jpm12 = (jp + sub_base_m + 48);
jpm13 = (jp + sub_base_m + 52);
jpm14 = (jp + sub_base_m + 56);
jpm15 = (jp + sub_base_m + 60);
for(; jp <= j + Bj - 4; jp += 64) {
iplnpjpm0 = ipln + jpm0;
c_v0 = _mm256_load_pd(C + iplnpjpm0);
b_v0 = _mm256_load_pd(B + kpbln + jpm0);
apb_v0 = _mm256_add_pd(a_v, b_v0);
res0 = _mm256_min_pd(c_v0, apb_v0);
_mm256_store_pd(C + iplnpjpm0, res0);
jpm0 += 64;
iplnpjpm1 = ipln + jpm1;
c_v1 = _mm256_load_pd(C + iplnpjpm1);
b_v1 = _mm256_load_pd(B + kpbln + jpm1);
apb_v1 = _mm256_add_pd(a_v, b_v1);
res1 = _mm256_min_pd(c_v1, apb_v1);
_mm256_store_pd(C + iplnpjpm1, res1);
jpm1 += 64;
iplnpjpm2 = ipln + jpm2;
c_v2 = _mm256_load_pd(C + iplnpjpm2);
b_v2 = _mm256_load_pd(B + kpbln + jpm2);
apb_v2 = _mm256_add_pd(a_v, b_v2);
res2 = _mm256_min_pd(c_v2, apb_v2);
_mm256_store_pd(C + iplnpjpm2, res2);
jpm2 += 64;
iplnpjpm3 = ipln + jpm3;
c_v3 = _mm256_load_pd(C + iplnpjpm3);
b_v3 = _mm256_load_pd(B + kpbln + jpm3);
apb_v3 = _mm256_add_pd(a_v, b_v3);
res3 = _mm256_min_pd(c_v3, apb_v3);
_mm256_store_pd(C + iplnpjpm3, res3);
jpm3 += 64;
iplnpjpm4 = ipln + jpm4;
c_v4 = _mm256_load_pd(C + iplnpjpm4);
b_v4 = _mm256_load_pd(B + kpbln + jpm4);
apb_v4 = _mm256_add_pd(a_v, b_v4);
res4 = _mm256_min_pd(c_v4, apb_v4);
_mm256_store_pd(C + iplnpjpm4, res4);
jpm4 += 64;
iplnpjpm5 = ipln + jpm5;
c_v5 = _mm256_load_pd(C + iplnpjpm5);
b_v5 = _mm256_load_pd(B + kpbln + jpm5);
apb_v5 = _mm256_add_pd(a_v, b_v5);
res5 = _mm256_min_pd(c_v5, apb_v5);
_mm256_store_pd(C + iplnpjpm5, res5);
jpm5 += 64;
iplnpjpm6 = ipln + jpm6;
c_v6 = _mm256_load_pd(C + iplnpjpm6);
b_v6 = _mm256_load_pd(B + kpbln + jpm6);
apb_v6 = _mm256_add_pd(a_v, b_v6);
res6 = _mm256_min_pd(c_v6, apb_v6);
_mm256_store_pd(C + iplnpjpm6, res6);
jpm6 += 64;
iplnpjpm7 = ipln + jpm7;
c_v7 = _mm256_load_pd(C + iplnpjpm7);
b_v7 = _mm256_load_pd(B + kpbln + jpm7);
apb_v7 = _mm256_add_pd(a_v, b_v7);
res7 = _mm256_min_pd(c_v7, apb_v7);
_mm256_store_pd(C + iplnpjpm7, res7);
jpm7 += 64;
iplnpjpm8 = ipln + jpm8;
c_v8 = _mm256_load_pd(C + iplnpjpm8);
b_v8 = _mm256_load_pd(B + kpbln + jpm8);
apb_v8 = _mm256_add_pd(a_v, b_v8);
res8 = _mm256_min_pd(c_v8, apb_v8);
_mm256_store_pd(C + iplnpjpm8, res8);
jpm8 += 64;
iplnpjpm9 = ipln + jpm9;
c_v9 = _mm256_load_pd(C + iplnpjpm9);
b_v9 = _mm256_load_pd(B + kpbln + jpm9);
apb_v9 = _mm256_add_pd(a_v, b_v9);
res9 = _mm256_min_pd(c_v9, apb_v9);
_mm256_store_pd(C + iplnpjpm9, res9);
jpm9 += 64;
iplnpjpm10 = ipln + jpm10;
c_v10 = _mm256_load_pd(C + iplnpjpm10);
b_v10 = _mm256_load_pd(B + kpbln + jpm10);
apb_v10 = _mm256_add_pd(a_v, b_v10);
res10 = _mm256_min_pd(c_v10, apb_v10);
_mm256_store_pd(C + iplnpjpm10, res10);
jpm10 += 64;
iplnpjpm11 = ipln + jpm11;
c_v11 = _mm256_load_pd(C + iplnpjpm11);
b_v11 = _mm256_load_pd(B + kpbln + jpm11);
apb_v11 = _mm256_add_pd(a_v, b_v11);
res11 = _mm256_min_pd(c_v11, apb_v11);
_mm256_store_pd(C + iplnpjpm11, res11);
jpm11 += 64;
iplnpjpm12 = ipln + jpm12;
c_v12 = _mm256_load_pd(C + iplnpjpm12);
b_v12 = _mm256_load_pd(B + kpbln + jpm12);
apb_v12 = _mm256_add_pd(a_v, b_v12);
res12 = _mm256_min_pd(c_v12, apb_v12);
_mm256_store_pd(C + iplnpjpm12, res12);
jpm12 += 64;
iplnpjpm13 = ipln + jpm13;
c_v13 = _mm256_load_pd(C + iplnpjpm13);
b_v13 = _mm256_load_pd(B + kpbln + jpm13);
apb_v13 = _mm256_add_pd(a_v, b_v13);
res13 = _mm256_min_pd(c_v13, apb_v13);
_mm256_store_pd(C + iplnpjpm13, res13);
jpm13 += 64;
iplnpjpm14 = ipln + jpm14;
c_v14 = _mm256_load_pd(C + iplnpjpm14);
b_v14 = _mm256_load_pd(B + kpbln + jpm14);
apb_v14 = _mm256_add_pd(a_v, b_v14);
res14 = _mm256_min_pd(c_v14, apb_v14);
_mm256_store_pd(C + iplnpjpm14, res14);
jpm14 += 64;
iplnpjpm15 = ipln + jpm15;
c_v15 = _mm256_load_pd(C + iplnpjpm15);
b_v15 = _mm256_load_pd(B + kpbln + jpm15);
apb_v15 = _mm256_add_pd(a_v, b_v15);
res15 = _mm256_min_pd(c_v15, apb_v15);
_mm256_store_pd(C + iplnpjpm15, res15);
jpm15 += 64;
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
int jpm2;int iplnjpm2;__m256d c_v2;__m256d b_v2;__m256d apb_v2;__m256d cmp_lt2;__m256d res2;double *Bkj2;
int jpm3;int iplnjpm3;__m256d c_v3;__m256d b_v3;__m256d apb_v3;__m256d cmp_lt3;__m256d res3;double *Bkj3;
int jpm4;int iplnjpm4;__m256d c_v4;__m256d b_v4;__m256d apb_v4;__m256d cmp_lt4;__m256d res4;double *Bkj4;
int jpm5;int iplnjpm5;__m256d c_v5;__m256d b_v5;__m256d apb_v5;__m256d cmp_lt5;__m256d res5;double *Bkj5;
int jpm6;int iplnjpm6;__m256d c_v6;__m256d b_v6;__m256d apb_v6;__m256d cmp_lt6;__m256d res6;double *Bkj6;
int jpm7;int iplnjpm7;__m256d c_v7;__m256d b_v7;__m256d apb_v7;__m256d cmp_lt7;__m256d res7;double *Bkj7;
int jpm8;int iplnjpm8;__m256d c_v8;__m256d b_v8;__m256d apb_v8;__m256d cmp_lt8;__m256d res8;double *Bkj8;
int jpm9;int iplnjpm9;__m256d c_v9;__m256d b_v9;__m256d apb_v9;__m256d cmp_lt9;__m256d res9;double *Bkj9;
int jpm10;int iplnjpm10;__m256d c_v10;__m256d b_v10;__m256d apb_v10;__m256d cmp_lt10;__m256d res10;double *Bkj10;
int jpm11;int iplnjpm11;__m256d c_v11;__m256d b_v11;__m256d apb_v11;__m256d cmp_lt11;__m256d res11;double *Bkj11;
int jpm12;int iplnjpm12;__m256d c_v12;__m256d b_v12;__m256d apb_v12;__m256d cmp_lt12;__m256d res12;double *Bkj12;
int jpm13;int iplnjpm13;__m256d c_v13;__m256d b_v13;__m256d apb_v13;__m256d cmp_lt13;__m256d res13;double *Bkj13;
int jpm14;int iplnjpm14;__m256d c_v14;__m256d b_v14;__m256d apb_v14;__m256d cmp_lt14;__m256d res14;double *Bkj14;
int jpm15;int iplnjpm15;__m256d c_v15;__m256d b_v15;__m256d apb_v15;__m256d cmp_lt15;__m256d res15;double *Bkj15;

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
jpm2 = (jp + sub_base_m + 8); iplnjpm2= ipln + jpm2; Bkj2=B + kln + jpm2;
jpm3 = (jp + sub_base_m + 12); iplnjpm3= ipln + jpm3; Bkj3=B + kln + jpm3;
jpm4 = (jp + sub_base_m + 16); iplnjpm4= ipln + jpm4; Bkj4=B + kln + jpm4;
jpm5 = (jp + sub_base_m + 20); iplnjpm5= ipln + jpm5; Bkj5=B + kln + jpm5;
jpm6 = (jp + sub_base_m + 24); iplnjpm6= ipln + jpm6; Bkj6=B + kln + jpm6;
jpm7 = (jp + sub_base_m + 28); iplnjpm7= ipln + jpm7; Bkj7=B + kln + jpm7;
jpm8 = (jp + sub_base_m + 32); iplnjpm8= ipln + jpm8; Bkj8=B + kln + jpm8;
jpm9 = (jp + sub_base_m + 36); iplnjpm9= ipln + jpm9; Bkj9=B + kln + jpm9;
jpm10 = (jp + sub_base_m + 40); iplnjpm10= ipln + jpm10; Bkj10=B + kln + jpm10;
jpm11 = (jp + sub_base_m + 44); iplnjpm11= ipln + jpm11; Bkj11=B + kln + jpm11;
jpm12 = (jp + sub_base_m + 48); iplnjpm12= ipln + jpm12; Bkj12=B + kln + jpm12;
jpm13 = (jp + sub_base_m + 52); iplnjpm13= ipln + jpm13; Bkj13=B + kln + jpm13;
jpm14 = (jp + sub_base_m + 56); iplnjpm14= ipln + jpm14; Bkj14=B + kln + jpm14;
jpm15 = (jp + sub_base_m + 60); iplnjpm15= ipln + jpm15; Bkj15=B + kln + jpm15;
for(; jp <= j + Bj - 64; jp += 64) {
b_v0 = _mm256_load_pd(Bkj0)
;apb_v0 = _mm256_add_pd(a_v, b_v0)
; c_v0 = _mm256_load_pd(C + iplnjpm0)
;res0 = _mm256_min_pd(c_v0, apb_v0);
_mm256_store_pd(C + iplnjpm0, res0)
; iplnjpm0 += 64
; Bkj0 += 64
;b_v1 = _mm256_load_pd(Bkj1)
;apb_v1 = _mm256_add_pd(a_v, b_v1)
; c_v1 = _mm256_load_pd(C + iplnjpm1)
;res1 = _mm256_min_pd(c_v1, apb_v1);
_mm256_store_pd(C + iplnjpm1, res1)
; iplnjpm1 += 64
; Bkj1 += 64
;b_v2 = _mm256_load_pd(Bkj2)
;apb_v2 = _mm256_add_pd(a_v, b_v2)
; c_v2 = _mm256_load_pd(C + iplnjpm2)
;res2 = _mm256_min_pd(c_v2, apb_v2);
_mm256_store_pd(C + iplnjpm2, res2)
; iplnjpm2 += 64
; Bkj2 += 64
;b_v3 = _mm256_load_pd(Bkj3)
;apb_v3 = _mm256_add_pd(a_v, b_v3)
; c_v3 = _mm256_load_pd(C + iplnjpm3)
;res3 = _mm256_min_pd(c_v3, apb_v3);
_mm256_store_pd(C + iplnjpm3, res3)
; iplnjpm3 += 64
; Bkj3 += 64
;b_v4 = _mm256_load_pd(Bkj4)
;apb_v4 = _mm256_add_pd(a_v, b_v4)
; c_v4 = _mm256_load_pd(C + iplnjpm4)
;res4 = _mm256_min_pd(c_v4, apb_v4);
_mm256_store_pd(C + iplnjpm4, res4)
; iplnjpm4 += 64
; Bkj4 += 64
;b_v5 = _mm256_load_pd(Bkj5)
;apb_v5 = _mm256_add_pd(a_v, b_v5)
; c_v5 = _mm256_load_pd(C + iplnjpm5)
;res5 = _mm256_min_pd(c_v5, apb_v5);
_mm256_store_pd(C + iplnjpm5, res5)
; iplnjpm5 += 64
; Bkj5 += 64
;b_v6 = _mm256_load_pd(Bkj6)
;apb_v6 = _mm256_add_pd(a_v, b_v6)
; c_v6 = _mm256_load_pd(C + iplnjpm6)
;res6 = _mm256_min_pd(c_v6, apb_v6);
_mm256_store_pd(C + iplnjpm6, res6)
; iplnjpm6 += 64
; Bkj6 += 64
;b_v7 = _mm256_load_pd(Bkj7)
;apb_v7 = _mm256_add_pd(a_v, b_v7)
; c_v7 = _mm256_load_pd(C + iplnjpm7)
;res7 = _mm256_min_pd(c_v7, apb_v7);
_mm256_store_pd(C + iplnjpm7, res7)
; iplnjpm7 += 64
; Bkj7 += 64
;b_v8 = _mm256_load_pd(Bkj8)
;apb_v8 = _mm256_add_pd(a_v, b_v8)
; c_v8 = _mm256_load_pd(C + iplnjpm8)
;res8 = _mm256_min_pd(c_v8, apb_v8);
_mm256_store_pd(C + iplnjpm8, res8)
; iplnjpm8 += 64
; Bkj8 += 64
;b_v9 = _mm256_load_pd(Bkj9)
;apb_v9 = _mm256_add_pd(a_v, b_v9)
; c_v9 = _mm256_load_pd(C + iplnjpm9)
;res9 = _mm256_min_pd(c_v9, apb_v9);
_mm256_store_pd(C + iplnjpm9, res9)
; iplnjpm9 += 64
; Bkj9 += 64
;b_v10 = _mm256_load_pd(Bkj10)
;apb_v10 = _mm256_add_pd(a_v, b_v10)
; c_v10 = _mm256_load_pd(C + iplnjpm10)
;res10 = _mm256_min_pd(c_v10, apb_v10);
_mm256_store_pd(C + iplnjpm10, res10)
; iplnjpm10 += 64
; Bkj10 += 64
;b_v11 = _mm256_load_pd(Bkj11)
;apb_v11 = _mm256_add_pd(a_v, b_v11)
; c_v11 = _mm256_load_pd(C + iplnjpm11)
;res11 = _mm256_min_pd(c_v11, apb_v11);
_mm256_store_pd(C + iplnjpm11, res11)
; iplnjpm11 += 64
; Bkj11 += 64
;b_v12 = _mm256_load_pd(Bkj12)
;apb_v12 = _mm256_add_pd(a_v, b_v12)
; c_v12 = _mm256_load_pd(C + iplnjpm12)
;res12 = _mm256_min_pd(c_v12, apb_v12);
_mm256_store_pd(C + iplnjpm12, res12)
; iplnjpm12 += 64
; Bkj12 += 64
;b_v13 = _mm256_load_pd(Bkj13)
;apb_v13 = _mm256_add_pd(a_v, b_v13)
; c_v13 = _mm256_load_pd(C + iplnjpm13)
;res13 = _mm256_min_pd(c_v13, apb_v13);
_mm256_store_pd(C + iplnjpm13, res13)
; iplnjpm13 += 64
; Bkj13 += 64
;b_v14 = _mm256_load_pd(Bkj14)
;apb_v14 = _mm256_add_pd(a_v, b_v14)
; c_v14 = _mm256_load_pd(C + iplnjpm14)
;res14 = _mm256_min_pd(c_v14, apb_v14);
_mm256_store_pd(C + iplnjpm14, res14)
; iplnjpm14 += 64
; Bkj14 += 64
;b_v15 = _mm256_load_pd(Bkj15)
;apb_v15 = _mm256_add_pd(a_v, b_v15)
; c_v15 = _mm256_load_pd(C + iplnjpm15)
;res15 = _mm256_min_pd(c_v15, apb_v15);
_mm256_store_pd(C + iplnjpm15, res15)
; iplnjpm15 += 64
; Bkj15 += 64
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
int jpl2;int ipmnjpl2;int klnjpl2;__m256d c_v2;__m256d b_v2;__m256d apb_v2;__m256d cmp_lt2;__m256d res2;
int jpl3;int ipmnjpl3;int klnjpl3;__m256d c_v3;__m256d b_v3;__m256d apb_v3;__m256d cmp_lt3;__m256d res3;
int jpl4;int ipmnjpl4;int klnjpl4;__m256d c_v4;__m256d b_v4;__m256d apb_v4;__m256d cmp_lt4;__m256d res4;
int jpl5;int ipmnjpl5;int klnjpl5;__m256d c_v5;__m256d b_v5;__m256d apb_v5;__m256d cmp_lt5;__m256d res5;
int jpl6;int ipmnjpl6;int klnjpl6;__m256d c_v6;__m256d b_v6;__m256d apb_v6;__m256d cmp_lt6;__m256d res6;
int jpl7;int ipmnjpl7;int klnjpl7;__m256d c_v7;__m256d b_v7;__m256d apb_v7;__m256d cmp_lt7;__m256d res7;
int jpl8;int ipmnjpl8;int klnjpl8;__m256d c_v8;__m256d b_v8;__m256d apb_v8;__m256d cmp_lt8;__m256d res8;
int jpl9;int ipmnjpl9;int klnjpl9;__m256d c_v9;__m256d b_v9;__m256d apb_v9;__m256d cmp_lt9;__m256d res9;
int jpl10;int ipmnjpl10;int klnjpl10;__m256d c_v10;__m256d b_v10;__m256d apb_v10;__m256d cmp_lt10;__m256d res10;
int jpl11;int ipmnjpl11;int klnjpl11;__m256d c_v11;__m256d b_v11;__m256d apb_v11;__m256d cmp_lt11;__m256d res11;
int jpl12;int ipmnjpl12;int klnjpl12;__m256d c_v12;__m256d b_v12;__m256d apb_v12;__m256d cmp_lt12;__m256d res12;
int jpl13;int ipmnjpl13;int klnjpl13;__m256d c_v13;__m256d b_v13;__m256d apb_v13;__m256d cmp_lt13;__m256d res13;
int jpl14;int ipmnjpl14;int klnjpl14;__m256d c_v14;__m256d b_v14;__m256d apb_v14;__m256d cmp_lt14;__m256d res14;
int jpl15;int ipmnjpl15;int klnjpl15;__m256d c_v15;__m256d b_v15;__m256d apb_v15;__m256d cmp_lt15;__m256d res15;

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
jpl2 = (jp + sub_base_l + 8);
jpl3 = (jp + sub_base_l + 12);
jpl4 = (jp + sub_base_l + 16);
jpl5 = (jp + sub_base_l + 20);
jpl6 = (jp + sub_base_l + 24);
jpl7 = (jp + sub_base_l + 28);
jpl8 = (jp + sub_base_l + 32);
jpl9 = (jp + sub_base_l + 36);
jpl10 = (jp + sub_base_l + 40);
jpl11 = (jp + sub_base_l + 44);
jpl12 = (jp + sub_base_l + 48);
jpl13 = (jp + sub_base_l + 52);
jpl14 = (jp + sub_base_l + 56);
jpl15 = (jp + sub_base_l + 60);
for(; jp <= j + Bj - 64; jp += 64) {ipmnjpl0 = ipmn + jpl0;
klnjpl0 = kln + jpl0;
b_v0 = _mm256_load_pd(B + klnjpl0);
apb_v0 = _mm256_add_pd(a_v, b_v0);
c_v0 = _mm256_load_pd(C + ipmnjpl0);
res0 = _mm256_min_pd(c_v0, apb_v0);
_mm256_store_pd(C + ipmnjpl0, res0);
jpl0 += 64;
ipmnjpl1 = ipmn + jpl1;
klnjpl1 = kln + jpl1;
b_v1 = _mm256_load_pd(B + klnjpl1);
apb_v1 = _mm256_add_pd(a_v, b_v1);
c_v1 = _mm256_load_pd(C + ipmnjpl1);
res1 = _mm256_min_pd(c_v1, apb_v1);
_mm256_store_pd(C + ipmnjpl1, res1);
jpl1 += 64;
ipmnjpl2 = ipmn + jpl2;
klnjpl2 = kln + jpl2;
b_v2 = _mm256_load_pd(B + klnjpl2);
apb_v2 = _mm256_add_pd(a_v, b_v2);
c_v2 = _mm256_load_pd(C + ipmnjpl2);
res2 = _mm256_min_pd(c_v2, apb_v2);
_mm256_store_pd(C + ipmnjpl2, res2);
jpl2 += 64;
ipmnjpl3 = ipmn + jpl3;
klnjpl3 = kln + jpl3;
b_v3 = _mm256_load_pd(B + klnjpl3);
apb_v3 = _mm256_add_pd(a_v, b_v3);
c_v3 = _mm256_load_pd(C + ipmnjpl3);
res3 = _mm256_min_pd(c_v3, apb_v3);
_mm256_store_pd(C + ipmnjpl3, res3);
jpl3 += 64;
ipmnjpl4 = ipmn + jpl4;
klnjpl4 = kln + jpl4;
b_v4 = _mm256_load_pd(B + klnjpl4);
apb_v4 = _mm256_add_pd(a_v, b_v4);
c_v4 = _mm256_load_pd(C + ipmnjpl4);
res4 = _mm256_min_pd(c_v4, apb_v4);
_mm256_store_pd(C + ipmnjpl4, res4);
jpl4 += 64;
ipmnjpl5 = ipmn + jpl5;
klnjpl5 = kln + jpl5;
b_v5 = _mm256_load_pd(B + klnjpl5);
apb_v5 = _mm256_add_pd(a_v, b_v5);
c_v5 = _mm256_load_pd(C + ipmnjpl5);
res5 = _mm256_min_pd(c_v5, apb_v5);
_mm256_store_pd(C + ipmnjpl5, res5);
jpl5 += 64;
ipmnjpl6 = ipmn + jpl6;
klnjpl6 = kln + jpl6;
b_v6 = _mm256_load_pd(B + klnjpl6);
apb_v6 = _mm256_add_pd(a_v, b_v6);
c_v6 = _mm256_load_pd(C + ipmnjpl6);
res6 = _mm256_min_pd(c_v6, apb_v6);
_mm256_store_pd(C + ipmnjpl6, res6);
jpl6 += 64;
ipmnjpl7 = ipmn + jpl7;
klnjpl7 = kln + jpl7;
b_v7 = _mm256_load_pd(B + klnjpl7);
apb_v7 = _mm256_add_pd(a_v, b_v7);
c_v7 = _mm256_load_pd(C + ipmnjpl7);
res7 = _mm256_min_pd(c_v7, apb_v7);
_mm256_store_pd(C + ipmnjpl7, res7);
jpl7 += 64;
ipmnjpl8 = ipmn + jpl8;
klnjpl8 = kln + jpl8;
b_v8 = _mm256_load_pd(B + klnjpl8);
apb_v8 = _mm256_add_pd(a_v, b_v8);
c_v8 = _mm256_load_pd(C + ipmnjpl8);
res8 = _mm256_min_pd(c_v8, apb_v8);
_mm256_store_pd(C + ipmnjpl8, res8);
jpl8 += 64;
ipmnjpl9 = ipmn + jpl9;
klnjpl9 = kln + jpl9;
b_v9 = _mm256_load_pd(B + klnjpl9);
apb_v9 = _mm256_add_pd(a_v, b_v9);
c_v9 = _mm256_load_pd(C + ipmnjpl9);
res9 = _mm256_min_pd(c_v9, apb_v9);
_mm256_store_pd(C + ipmnjpl9, res9);
jpl9 += 64;
ipmnjpl10 = ipmn + jpl10;
klnjpl10 = kln + jpl10;
b_v10 = _mm256_load_pd(B + klnjpl10);
apb_v10 = _mm256_add_pd(a_v, b_v10);
c_v10 = _mm256_load_pd(C + ipmnjpl10);
res10 = _mm256_min_pd(c_v10, apb_v10);
_mm256_store_pd(C + ipmnjpl10, res10);
jpl10 += 64;
ipmnjpl11 = ipmn + jpl11;
klnjpl11 = kln + jpl11;
b_v11 = _mm256_load_pd(B + klnjpl11);
apb_v11 = _mm256_add_pd(a_v, b_v11);
c_v11 = _mm256_load_pd(C + ipmnjpl11);
res11 = _mm256_min_pd(c_v11, apb_v11);
_mm256_store_pd(C + ipmnjpl11, res11);
jpl11 += 64;
ipmnjpl12 = ipmn + jpl12;
klnjpl12 = kln + jpl12;
b_v12 = _mm256_load_pd(B + klnjpl12);
apb_v12 = _mm256_add_pd(a_v, b_v12);
c_v12 = _mm256_load_pd(C + ipmnjpl12);
res12 = _mm256_min_pd(c_v12, apb_v12);
_mm256_store_pd(C + ipmnjpl12, res12);
jpl12 += 64;
ipmnjpl13 = ipmn + jpl13;
klnjpl13 = kln + jpl13;
b_v13 = _mm256_load_pd(B + klnjpl13);
apb_v13 = _mm256_add_pd(a_v, b_v13);
c_v13 = _mm256_load_pd(C + ipmnjpl13);
res13 = _mm256_min_pd(c_v13, apb_v13);
_mm256_store_pd(C + ipmnjpl13, res13);
jpl13 += 64;
ipmnjpl14 = ipmn + jpl14;
klnjpl14 = kln + jpl14;
b_v14 = _mm256_load_pd(B + klnjpl14);
apb_v14 = _mm256_add_pd(a_v, b_v14);
c_v14 = _mm256_load_pd(C + ipmnjpl14);
res14 = _mm256_min_pd(c_v14, apb_v14);
_mm256_store_pd(C + ipmnjpl14, res14);
jpl14 += 64;
ipmnjpl15 = ipmn + jpl15;
klnjpl15 = kln + jpl15;
b_v15 = _mm256_load_pd(B + klnjpl15);
apb_v15 = _mm256_add_pd(a_v, b_v15);
c_v15 = _mm256_load_pd(C + ipmnjpl15);
res15 = _mm256_min_pd(c_v15, apb_v15);
_mm256_store_pd(C + ipmnjpl15, res15);
jpl15 += 64;
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
int jp2;int jpsubo2;__m256d c_v2;__m256d b_v2;__m256d apb_v2;__m256d cmp_lt2;__m256d res2;double *B_k_j2;double *C_i_j2;
int jp3;int jpsubo3;__m256d c_v3;__m256d b_v3;__m256d apb_v3;__m256d cmp_lt3;__m256d res3;double *B_k_j3;double *C_i_j3;
int jp4;int jpsubo4;__m256d c_v4;__m256d b_v4;__m256d apb_v4;__m256d cmp_lt4;__m256d res4;double *B_k_j4;double *C_i_j4;
int jp5;int jpsubo5;__m256d c_v5;__m256d b_v5;__m256d apb_v5;__m256d cmp_lt5;__m256d res5;double *B_k_j5;double *C_i_j5;
int jp6;int jpsubo6;__m256d c_v6;__m256d b_v6;__m256d apb_v6;__m256d cmp_lt6;__m256d res6;double *B_k_j6;double *C_i_j6;
int jp7;int jpsubo7;__m256d c_v7;__m256d b_v7;__m256d apb_v7;__m256d cmp_lt7;__m256d res7;double *B_k_j7;double *C_i_j7;
int jp8;int jpsubo8;__m256d c_v8;__m256d b_v8;__m256d apb_v8;__m256d cmp_lt8;__m256d res8;double *B_k_j8;double *C_i_j8;
int jp9;int jpsubo9;__m256d c_v9;__m256d b_v9;__m256d apb_v9;__m256d cmp_lt9;__m256d res9;double *B_k_j9;double *C_i_j9;
int jp10;int jpsubo10;__m256d c_v10;__m256d b_v10;__m256d apb_v10;__m256d cmp_lt10;__m256d res10;double *B_k_j10;double *C_i_j10;
int jp11;int jpsubo11;__m256d c_v11;__m256d b_v11;__m256d apb_v11;__m256d cmp_lt11;__m256d res11;double *B_k_j11;double *C_i_j11;
int jp12;int jpsubo12;__m256d c_v12;__m256d b_v12;__m256d apb_v12;__m256d cmp_lt12;__m256d res12;double *B_k_j12;double *C_i_j12;
int jp13;int jpsubo13;__m256d c_v13;__m256d b_v13;__m256d apb_v13;__m256d cmp_lt13;__m256d res13;double *B_k_j13;double *C_i_j13;
int jp14;int jpsubo14;__m256d c_v14;__m256d b_v14;__m256d apb_v14;__m256d cmp_lt14;__m256d res14;double *B_k_j14;double *C_i_j14;
int jp15;int jpsubo15;__m256d c_v15;__m256d b_v15;__m256d apb_v15;__m256d cmp_lt15;__m256d res15;double *B_k_j15;double *C_i_j15;

                            for (int i = 0; i < L1; i += Bi) {
                                iBi = i + Bi;
                                for (int j = 0; j < L1; j += Bj) {
                                    jBj = j + Bj;
                                    jBj4 = jBj - 64;
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
jpsubo2 = (sub_base_o + 8);
jpsubo3 = (sub_base_o + 12);
jpsubo4 = (sub_base_o + 16);
jpsubo5 = (sub_base_o + 20);
jpsubo6 = (sub_base_o + 24);
jpsubo7 = (sub_base_o + 28);
jpsubo8 = (sub_base_o + 32);
jpsubo9 = (sub_base_o + 36);
jpsubo10 = (sub_base_o + 40);
jpsubo11 = (sub_base_o + 44);
jpsubo12 = (sub_base_o + 48);
jpsubo13 = (sub_base_o + 52);
jpsubo14 = (sub_base_o + 56);
jpsubo15 = (sub_base_o + 60);
for(; jp <= jBj4; jp += 64) {B_k_j0 = B_k + jpsubo0; C_i_j0 = C_i + jpsubo0;
b_v0 = _mm256_load_pd(B_k_j0);
c_v0 = _mm256_load_pd(C_i_j0);
apb_v0 = _mm256_add_pd(a_v, b_v0);
res0 = _mm256_min_pd(c_v0, apb_v0);
_mm256_store_pd(C_i_j0, res0);
jpsubo0+=64;
B_k_j1 = B_k + jpsubo1; C_i_j1 = C_i + jpsubo1;
b_v1 = _mm256_load_pd(B_k_j1);
c_v1 = _mm256_load_pd(C_i_j1);
apb_v1 = _mm256_add_pd(a_v, b_v1);
res1 = _mm256_min_pd(c_v1, apb_v1);
_mm256_store_pd(C_i_j1, res1);
jpsubo1+=64;
B_k_j2 = B_k + jpsubo2; C_i_j2 = C_i + jpsubo2;
b_v2 = _mm256_load_pd(B_k_j2);
c_v2 = _mm256_load_pd(C_i_j2);
apb_v2 = _mm256_add_pd(a_v, b_v2);
res2 = _mm256_min_pd(c_v2, apb_v2);
_mm256_store_pd(C_i_j2, res2);
jpsubo2+=64;
B_k_j3 = B_k + jpsubo3; C_i_j3 = C_i + jpsubo3;
b_v3 = _mm256_load_pd(B_k_j3);
c_v3 = _mm256_load_pd(C_i_j3);
apb_v3 = _mm256_add_pd(a_v, b_v3);
res3 = _mm256_min_pd(c_v3, apb_v3);
_mm256_store_pd(C_i_j3, res3);
jpsubo3+=64;
B_k_j4 = B_k + jpsubo4; C_i_j4 = C_i + jpsubo4;
b_v4 = _mm256_load_pd(B_k_j4);
c_v4 = _mm256_load_pd(C_i_j4);
apb_v4 = _mm256_add_pd(a_v, b_v4);
res4 = _mm256_min_pd(c_v4, apb_v4);
_mm256_store_pd(C_i_j4, res4);
jpsubo4+=64;
B_k_j5 = B_k + jpsubo5; C_i_j5 = C_i + jpsubo5;
b_v5 = _mm256_load_pd(B_k_j5);
c_v5 = _mm256_load_pd(C_i_j5);
apb_v5 = _mm256_add_pd(a_v, b_v5);
res5 = _mm256_min_pd(c_v5, apb_v5);
_mm256_store_pd(C_i_j5, res5);
jpsubo5+=64;
B_k_j6 = B_k + jpsubo6; C_i_j6 = C_i + jpsubo6;
b_v6 = _mm256_load_pd(B_k_j6);
c_v6 = _mm256_load_pd(C_i_j6);
apb_v6 = _mm256_add_pd(a_v, b_v6);
res6 = _mm256_min_pd(c_v6, apb_v6);
_mm256_store_pd(C_i_j6, res6);
jpsubo6+=64;
B_k_j7 = B_k + jpsubo7; C_i_j7 = C_i + jpsubo7;
b_v7 = _mm256_load_pd(B_k_j7);
c_v7 = _mm256_load_pd(C_i_j7);
apb_v7 = _mm256_add_pd(a_v, b_v7);
res7 = _mm256_min_pd(c_v7, apb_v7);
_mm256_store_pd(C_i_j7, res7);
jpsubo7+=64;
B_k_j8 = B_k + jpsubo8; C_i_j8 = C_i + jpsubo8;
b_v8 = _mm256_load_pd(B_k_j8);
c_v8 = _mm256_load_pd(C_i_j8);
apb_v8 = _mm256_add_pd(a_v, b_v8);
res8 = _mm256_min_pd(c_v8, apb_v8);
_mm256_store_pd(C_i_j8, res8);
jpsubo8+=64;
B_k_j9 = B_k + jpsubo9; C_i_j9 = C_i + jpsubo9;
b_v9 = _mm256_load_pd(B_k_j9);
c_v9 = _mm256_load_pd(C_i_j9);
apb_v9 = _mm256_add_pd(a_v, b_v9);
res9 = _mm256_min_pd(c_v9, apb_v9);
_mm256_store_pd(C_i_j9, res9);
jpsubo9+=64;
B_k_j10 = B_k + jpsubo10; C_i_j10 = C_i + jpsubo10;
b_v10 = _mm256_load_pd(B_k_j10);
c_v10 = _mm256_load_pd(C_i_j10);
apb_v10 = _mm256_add_pd(a_v, b_v10);
res10 = _mm256_min_pd(c_v10, apb_v10);
_mm256_store_pd(C_i_j10, res10);
jpsubo10+=64;
B_k_j11 = B_k + jpsubo11; C_i_j11 = C_i + jpsubo11;
b_v11 = _mm256_load_pd(B_k_j11);
c_v11 = _mm256_load_pd(C_i_j11);
apb_v11 = _mm256_add_pd(a_v, b_v11);
res11 = _mm256_min_pd(c_v11, apb_v11);
_mm256_store_pd(C_i_j11, res11);
jpsubo11+=64;
B_k_j12 = B_k + jpsubo12; C_i_j12 = C_i + jpsubo12;
b_v12 = _mm256_load_pd(B_k_j12);
c_v12 = _mm256_load_pd(C_i_j12);
apb_v12 = _mm256_add_pd(a_v, b_v12);
res12 = _mm256_min_pd(c_v12, apb_v12);
_mm256_store_pd(C_i_j12, res12);
jpsubo12+=64;
B_k_j13 = B_k + jpsubo13; C_i_j13 = C_i + jpsubo13;
b_v13 = _mm256_load_pd(B_k_j13);
c_v13 = _mm256_load_pd(C_i_j13);
apb_v13 = _mm256_add_pd(a_v, b_v13);
res13 = _mm256_min_pd(c_v13, apb_v13);
_mm256_store_pd(C_i_j13, res13);
jpsubo13+=64;
B_k_j14 = B_k + jpsubo14; C_i_j14 = C_i + jpsubo14;
b_v14 = _mm256_load_pd(B_k_j14);
c_v14 = _mm256_load_pd(C_i_j14);
apb_v14 = _mm256_add_pd(a_v, b_v14);
res14 = _mm256_min_pd(c_v14, apb_v14);
_mm256_store_pd(C_i_j14, res14);
jpsubo14+=64;
B_k_j15 = B_k + jpsubo15; C_i_j15 = C_i + jpsubo15;
b_v15 = _mm256_load_pd(B_k_j15);
c_v15 = _mm256_load_pd(C_i_j15);
apb_v15 = _mm256_add_pd(a_v, b_v15);
res15 = _mm256_min_pd(c_v15, apb_v15);
_mm256_store_pd(C_i_j15, res15);
jpsubo15+=64;
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
        double r1 = benchmark_tiled_timed(n, fw_abc_min_plus, opt_tiled_fw_min_plus, L1, Bi, Bj, Bk);
    #endif

        return 0;
    }
    