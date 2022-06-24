
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
        int jpm3; int iplnpjpm3; __m256d c_v3; __m256d b_v3; __m256d apb_v3; __m256d cmp_lt3; __m256d res3; 
        int jpm4; int iplnpjpm4; __m256d c_v4; __m256d b_v4; __m256d apb_v4; __m256d cmp_lt4; __m256d res4; 
        int jpm5; int iplnpjpm5; __m256d c_v5; __m256d b_v5; __m256d apb_v5; __m256d cmp_lt5; __m256d res5; 
        int jpm6; int iplnpjpm6; __m256d c_v6; __m256d b_v6; __m256d apb_v6; __m256d cmp_lt6; __m256d res6; 
        int jpm7; int iplnpjpm7; __m256d c_v7; __m256d b_v7; __m256d apb_v7; __m256d cmp_lt7; __m256d res7; 
        int jpm8; int iplnpjpm8; __m256d c_v8; __m256d b_v8; __m256d apb_v8; __m256d cmp_lt8; __m256d res8; 
         
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
                        jpm3 = (jp + sub_base_m + 12);
                        jpm4 = (jp + sub_base_m + 16);
                        jpm5 = (jp + sub_base_m + 20);
                        jpm6 = (jp + sub_base_m + 24);
                        jpm7 = (jp + sub_base_m + 28);
                        jpm8 = (jp + sub_base_m + 32);

                        for(; jp <= j + Bj - 4; jp += 36) {
                            iplnpjpm0 = ipln + jpm0;
                            iplnpjpm1 = ipln + jpm1;
                            iplnpjpm2 = ipln + jpm2;
                            iplnpjpm3 = ipln + jpm3;
                            iplnpjpm4 = ipln + jpm4;
                            iplnpjpm5 = ipln + jpm5;
                            iplnpjpm6 = ipln + jpm6;
                            iplnpjpm7 = ipln + jpm7;
                            iplnpjpm8 = ipln + jpm8;
                            c_v0 = _mm256_load_pd(C + iplnpjpm0);
                            c_v1 = _mm256_load_pd(C + iplnpjpm1);
                            c_v2 = _mm256_load_pd(C + iplnpjpm2);
                            c_v3 = _mm256_load_pd(C + iplnpjpm3);
                            c_v4 = _mm256_load_pd(C + iplnpjpm4);
                            c_v5 = _mm256_load_pd(C + iplnpjpm5);
                            c_v6 = _mm256_load_pd(C + iplnpjpm6);
                            c_v7 = _mm256_load_pd(C + iplnpjpm7);
                            c_v8 = _mm256_load_pd(C + iplnpjpm8);
                            b_v0 = _mm256_load_pd(B + kpbln + jpm0);
                            b_v1 = _mm256_load_pd(B + kpbln + jpm1);
                            b_v2 = _mm256_load_pd(B + kpbln + jpm2);
                            b_v3 = _mm256_load_pd(B + kpbln + jpm3);
                            b_v4 = _mm256_load_pd(B + kpbln + jpm4);
                            b_v5 = _mm256_load_pd(B + kpbln + jpm5);
                            b_v6 = _mm256_load_pd(B + kpbln + jpm6);
                            b_v7 = _mm256_load_pd(B + kpbln + jpm7);
                            b_v8 = _mm256_load_pd(B + kpbln + jpm8);
                            apb_v0 = _mm256_add_pd(a_v, b_v0);
                            apb_v1 = _mm256_add_pd(a_v, b_v1);
                            apb_v2 = _mm256_add_pd(a_v, b_v2);
                            apb_v3 = _mm256_add_pd(a_v, b_v3);
                            apb_v4 = _mm256_add_pd(a_v, b_v4);
                            apb_v5 = _mm256_add_pd(a_v, b_v5);
                            apb_v6 = _mm256_add_pd(a_v, b_v6);
                            apb_v7 = _mm256_add_pd(a_v, b_v7);
                            apb_v8 = _mm256_add_pd(a_v, b_v8);
                            res0 = _mm256_min_pd(c_v0, apb_v0);
                            res1 = _mm256_min_pd(c_v1, apb_v1);
                            res2 = _mm256_min_pd(c_v2, apb_v2);
                            res3 = _mm256_min_pd(c_v3, apb_v3);
                            res4 = _mm256_min_pd(c_v4, apb_v4);
                            res5 = _mm256_min_pd(c_v5, apb_v5);
                            res6 = _mm256_min_pd(c_v6, apb_v6);
                            res7 = _mm256_min_pd(c_v7, apb_v7);
                            res8 = _mm256_min_pd(c_v8, apb_v8);
                            _mm256_store_pd(C + iplnpjpm0, res0);
                            _mm256_store_pd(C + iplnpjpm1, res1);
                            _mm256_store_pd(C + iplnpjpm2, res2);
                            _mm256_store_pd(C + iplnpjpm3, res3);
                            _mm256_store_pd(C + iplnpjpm4, res4);
                            _mm256_store_pd(C + iplnpjpm5, res5);
                            _mm256_store_pd(C + iplnpjpm6, res6);
                            _mm256_store_pd(C + iplnpjpm7, res7);
                            _mm256_store_pd(C + iplnpjpm8, res8);
                            jpm0 += 36;
                            jpm1 += 36;
                            jpm2 += 36;
                            jpm3 += 36;
                            jpm4 += 36;
                            jpm5 += 36;
                            jpm6 += 36;
                            jpm7 += 36;
                            jpm8 += 36;
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
                int jpm0; int iplnjpm0; __m256d c_v0; __m256d b_v0; __m256d apb_v0; __m256d cmp_lt0; __m256d res0; double *Bkj0; 
                int jpm1; int iplnjpm1; __m256d c_v1; __m256d b_v1; __m256d apb_v1; __m256d cmp_lt1; __m256d res1; double *Bkj1; 
                int jpm2; int iplnjpm2; __m256d c_v2; __m256d b_v2; __m256d apb_v2; __m256d cmp_lt2; __m256d res2; double *Bkj2; 
                int jpm3; int iplnjpm3; __m256d c_v3; __m256d b_v3; __m256d apb_v3; __m256d cmp_lt3; __m256d res3; double *Bkj3; 
                int jpm4; int iplnjpm4; __m256d c_v4; __m256d b_v4; __m256d apb_v4; __m256d cmp_lt4; __m256d res4; double *Bkj4; 
                int jpm5; int iplnjpm5; __m256d c_v5; __m256d b_v5; __m256d apb_v5; __m256d cmp_lt5; __m256d res5; double *Bkj5; 
                int jpm6; int iplnjpm6; __m256d c_v6; __m256d b_v6; __m256d apb_v6; __m256d cmp_lt6; __m256d res6; double *Bkj6; 
                int jpm7; int iplnjpm7; __m256d c_v7; __m256d b_v7; __m256d apb_v7; __m256d cmp_lt7; __m256d res7; double *Bkj7; 
                int jpm8; int iplnjpm8; __m256d c_v8; __m256d b_v8; __m256d apb_v8; __m256d cmp_lt8; __m256d res8; double *Bkj8; 

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
                                jpm3 = (jp + sub_base_m + 12); iplnjpm3= ipln + jpm3; Bkj3=B + kln + jpm3;
                                jpm4 = (jp + sub_base_m + 16); iplnjpm4= ipln + jpm4; Bkj4=B + kln + jpm4;
                                jpm5 = (jp + sub_base_m + 20); iplnjpm5= ipln + jpm5; Bkj5=B + kln + jpm5;
                                jpm6 = (jp + sub_base_m + 24); iplnjpm6= ipln + jpm6; Bkj6=B + kln + jpm6;
                                jpm7 = (jp + sub_base_m + 28); iplnjpm7= ipln + jpm7; Bkj7=B + kln + jpm7;
                                jpm8 = (jp + sub_base_m + 32); iplnjpm8= ipln + jpm8; Bkj8=B + kln + jpm8;

                                for(; jp <= j + Bj - 36; jp += 36) {
                                    b_v0 = _mm256_load_pd(Bkj0);
                                    b_v1 = _mm256_load_pd(Bkj1);
                                    b_v2 = _mm256_load_pd(Bkj2);
                                    b_v3 = _mm256_load_pd(Bkj3);
                                    b_v4 = _mm256_load_pd(Bkj4);
                                    b_v5 = _mm256_load_pd(Bkj5);
                                    b_v6 = _mm256_load_pd(Bkj6);
                                    b_v7 = _mm256_load_pd(Bkj7);
                                    b_v8 = _mm256_load_pd(Bkj8);
                                    c_v0 = _mm256_load_pd(C + iplnjpm0);
                                    c_v1 = _mm256_load_pd(C + iplnjpm1);
                                    c_v2 = _mm256_load_pd(C + iplnjpm2);
                                    c_v3 = _mm256_load_pd(C + iplnjpm3);
                                    c_v4 = _mm256_load_pd(C + iplnjpm4);
                                    c_v5 = _mm256_load_pd(C + iplnjpm5);
                                    c_v6 = _mm256_load_pd(C + iplnjpm6);
                                    c_v7 = _mm256_load_pd(C + iplnjpm7);
                                    c_v8 = _mm256_load_pd(C + iplnjpm8);
                                    apb_v0 = _mm256_add_pd(a_v, b_v0);
                                    apb_v1 = _mm256_add_pd(a_v, b_v1);
                                    apb_v2 = _mm256_add_pd(a_v, b_v2);
                                    apb_v3 = _mm256_add_pd(a_v, b_v3);
                                    apb_v4 = _mm256_add_pd(a_v, b_v4);
                                    apb_v5 = _mm256_add_pd(a_v, b_v5);
                                    apb_v6 = _mm256_add_pd(a_v, b_v6);
                                    apb_v7 = _mm256_add_pd(a_v, b_v7);
                                    apb_v8 = _mm256_add_pd(a_v, b_v8);
                                    res0 = _mm256_min_pd(c_v0, apb_v0);
                                    res1 = _mm256_min_pd(c_v1, apb_v1);
                                    res2 = _mm256_min_pd(c_v2, apb_v2);
                                    res3 = _mm256_min_pd(c_v3, apb_v3);
                                    res4 = _mm256_min_pd(c_v4, apb_v4);
                                    res5 = _mm256_min_pd(c_v5, apb_v5);
                                    res6 = _mm256_min_pd(c_v6, apb_v6);
                                    res7 = _mm256_min_pd(c_v7, apb_v7);
                                    res8 = _mm256_min_pd(c_v8, apb_v8);
                                    _mm256_store_pd(C + iplnjpm0, res0);
                                    _mm256_store_pd(C + iplnjpm1, res1);
                                    _mm256_store_pd(C + iplnjpm2, res2);
                                    _mm256_store_pd(C + iplnjpm3, res3);
                                    _mm256_store_pd(C + iplnjpm4, res4);
                                    _mm256_store_pd(C + iplnjpm5, res5);
                                    _mm256_store_pd(C + iplnjpm6, res6);
                                    _mm256_store_pd(C + iplnjpm7, res7);
                                    _mm256_store_pd(C + iplnjpm8, res8);
                                    iplnjpm0 += 36;
                                    iplnjpm1 += 36;
                                    iplnjpm2 += 36;
                                    iplnjpm3 += 36;
                                    iplnjpm4 += 36;
                                    iplnjpm5 += 36;
                                    iplnjpm6 += 36;
                                    iplnjpm7 += 36;
                                    iplnjpm8 += 36;
                                    Bkj0 += 36;
                                    Bkj1 += 36;
                                    Bkj2 += 36;
                                    Bkj3 += 36;
                                    Bkj4 += 36;
                                    Bkj5 += 36;
                                    Bkj6 += 36;
                                    Bkj7 += 36;
                                    Bkj8 += 36;
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
                double apb = 0.0;
                double c = 0.0;
                double min_c = 0.0;
                double a = 0.0;
                __m256d a_v;
                int jpl0; int ipmnjpl0; int klnjpl0; __m256d c_v0; __m256d b_v0; __m256d apb_v0; __m256d cmp_lt0; __m256d res0; 
                int jpl1; int ipmnjpl1; int klnjpl1; __m256d c_v1; __m256d b_v1; __m256d apb_v1; __m256d cmp_lt1; __m256d res1; 
                int jpl2; int ipmnjpl2; int klnjpl2; __m256d c_v2; __m256d b_v2; __m256d apb_v2; __m256d cmp_lt2; __m256d res2; 
                int jpl3; int ipmnjpl3; int klnjpl3; __m256d c_v3; __m256d b_v3; __m256d apb_v3; __m256d cmp_lt3; __m256d res3; 
                int jpl4; int ipmnjpl4; int klnjpl4; __m256d c_v4; __m256d b_v4; __m256d apb_v4; __m256d cmp_lt4; __m256d res4; 
                int jpl5; int ipmnjpl5; int klnjpl5; __m256d c_v5; __m256d b_v5; __m256d apb_v5; __m256d cmp_lt5; __m256d res5; 
                int jpl6; int ipmnjpl6; int klnjpl6; __m256d c_v6; __m256d b_v6; __m256d apb_v6; __m256d cmp_lt6; __m256d res6; 
                int jpl7; int ipmnjpl7; int klnjpl7; __m256d c_v7; __m256d b_v7; __m256d apb_v7; __m256d cmp_lt7; __m256d res7; 
                int jpl8; int ipmnjpl8; int klnjpl8; __m256d c_v8; __m256d b_v8; __m256d apb_v8; __m256d cmp_lt8; __m256d res8; 

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
                                jpl3 = (jp + sub_base_l + 12);
                                jpl4 = (jp + sub_base_l + 16);
                                jpl5 = (jp + sub_base_l + 20);
                                jpl6 = (jp + sub_base_l + 24);
                                jpl7 = (jp + sub_base_l + 28);
                                jpl8 = (jp + sub_base_l + 32);

                                for(; jp <= j + Bj - 36; jp += 36) {
                                    ipmnjpl0 = ipmn + jpl0;
                                    ipmnjpl1 = ipmn + jpl1;
                                    ipmnjpl2 = ipmn + jpl2;
                                    ipmnjpl3 = ipmn + jpl3;
                                    ipmnjpl4 = ipmn + jpl4;
                                    ipmnjpl5 = ipmn + jpl5;
                                    ipmnjpl6 = ipmn + jpl6;
                                    ipmnjpl7 = ipmn + jpl7;
                                    ipmnjpl8 = ipmn + jpl8;
                                    klnjpl0 = kln + jpl0;
                                    klnjpl1 = kln + jpl1;
                                    klnjpl2 = kln + jpl2;
                                    klnjpl3 = kln + jpl3;
                                    klnjpl4 = kln + jpl4;
                                    klnjpl5 = kln + jpl5;
                                    klnjpl6 = kln + jpl6;
                                    klnjpl7 = kln + jpl7;
                                    klnjpl8 = kln + jpl8;
                                    b_v0 = _mm256_load_pd(B + klnjpl0);
                                    b_v1 = _mm256_load_pd(B + klnjpl1);
                                    b_v2 = _mm256_load_pd(B + klnjpl2);
                                    b_v3 = _mm256_load_pd(B + klnjpl3);
                                    b_v4 = _mm256_load_pd(B + klnjpl4);
                                    b_v5 = _mm256_load_pd(B + klnjpl5);
                                    b_v6 = _mm256_load_pd(B + klnjpl6);
                                    b_v7 = _mm256_load_pd(B + klnjpl7);
                                    b_v8 = _mm256_load_pd(B + klnjpl8);
                                    c_v0 = _mm256_load_pd(C + ipmnjpl0);
                                    c_v1 = _mm256_load_pd(C + ipmnjpl1);
                                    c_v2 = _mm256_load_pd(C + ipmnjpl2);
                                    c_v3 = _mm256_load_pd(C + ipmnjpl3);
                                    c_v4 = _mm256_load_pd(C + ipmnjpl4);
                                    c_v5 = _mm256_load_pd(C + ipmnjpl5);
                                    c_v6 = _mm256_load_pd(C + ipmnjpl6);
                                    c_v7 = _mm256_load_pd(C + ipmnjpl7);
                                    c_v8 = _mm256_load_pd(C + ipmnjpl8);
                                    apb_v0 = _mm256_add_pd(a_v, b_v0);
                                    apb_v1 = _mm256_add_pd(a_v, b_v1);
                                    apb_v2 = _mm256_add_pd(a_v, b_v2);
                                    apb_v3 = _mm256_add_pd(a_v, b_v3);
                                    apb_v4 = _mm256_add_pd(a_v, b_v4);
                                    apb_v5 = _mm256_add_pd(a_v, b_v5);
                                    apb_v6 = _mm256_add_pd(a_v, b_v6);
                                    apb_v7 = _mm256_add_pd(a_v, b_v7);
                                    apb_v8 = _mm256_add_pd(a_v, b_v8);
                                    res0 = _mm256_min_pd(c_v0, apb_v0);
                                    res1 = _mm256_min_pd(c_v1, apb_v1);
                                    res2 = _mm256_min_pd(c_v2, apb_v2);
                                    res3 = _mm256_min_pd(c_v3, apb_v3);
                                    res4 = _mm256_min_pd(c_v4, apb_v4);
                                    res5 = _mm256_min_pd(c_v5, apb_v5);
                                    res6 = _mm256_min_pd(c_v6, apb_v6);
                                    res7 = _mm256_min_pd(c_v7, apb_v7);
                                    res8 = _mm256_min_pd(c_v8, apb_v8);
                                    _mm256_store_pd(C + ipmnjpl0, res0);
                                    _mm256_store_pd(C + ipmnjpl1, res1);
                                    _mm256_store_pd(C + ipmnjpl2, res2);
                                    _mm256_store_pd(C + ipmnjpl3, res3);
                                    _mm256_store_pd(C + ipmnjpl4, res4);
                                    _mm256_store_pd(C + ipmnjpl5, res5);
                                    _mm256_store_pd(C + ipmnjpl6, res6);
                                    _mm256_store_pd(C + ipmnjpl7, res7);
                                    _mm256_store_pd(C + ipmnjpl8, res8);
                                    jpl0 += 36;
                                    jpl1 += 36;
                                    jpl2 += 36;
                                    jpl3 += 36;
                                    jpl4 += 36;
                                    jpl5 += 36;
                                    jpl6 += 36;
                                    jpl7 += 36;
                                    jpl8 += 36;
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
                        int jp0; int jpsubo0; __m256d c_v0; __m256d b_v0; __m256d apb_v0; __m256d cmp_lt0; __m256d res0; double *B_k_j0; double *C_i_j0; 
                        int jp1; int jpsubo1; __m256d c_v1; __m256d b_v1; __m256d apb_v1; __m256d cmp_lt1; __m256d res1; double *B_k_j1; double *C_i_j1; 
                        int jp2; int jpsubo2; __m256d c_v2; __m256d b_v2; __m256d apb_v2; __m256d cmp_lt2; __m256d res2; double *B_k_j2; double *C_i_j2; 
                        int jp3; int jpsubo3; __m256d c_v3; __m256d b_v3; __m256d apb_v3; __m256d cmp_lt3; __m256d res3; double *B_k_j3; double *C_i_j3; 
                        int jp4; int jpsubo4; __m256d c_v4; __m256d b_v4; __m256d apb_v4; __m256d cmp_lt4; __m256d res4; double *B_k_j4; double *C_i_j4; 
                        int jp5; int jpsubo5; __m256d c_v5; __m256d b_v5; __m256d apb_v5; __m256d cmp_lt5; __m256d res5; double *B_k_j5; double *C_i_j5; 
                        int jp6; int jpsubo6; __m256d c_v6; __m256d b_v6; __m256d apb_v6; __m256d cmp_lt6; __m256d res6; double *B_k_j6; double *C_i_j6; 
                        int jp7; int jpsubo7; __m256d c_v7; __m256d b_v7; __m256d apb_v7; __m256d cmp_lt7; __m256d res7; double *B_k_j7; double *C_i_j7; 
                        int jp8; int jpsubo8; __m256d c_v8; __m256d b_v8; __m256d apb_v8; __m256d cmp_lt8; __m256d res8; double *B_k_j8; double *C_i_j8; 

                        for(int i = 0; i < L1; i += Bi) {
                            iBi = i + Bi;
                            for(int j = 0; j < L1; j += Bj) {
                                jBj = j + Bj;
                                jBj4 = jBj - 36;
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
                                            jpsubo3 = (sub_base_o + 12);
                                            jpsubo4 = (sub_base_o + 16);
                                            jpsubo5 = (sub_base_o + 20);
                                            jpsubo6 = (sub_base_o + 24);
                                            jpsubo7 = (sub_base_o + 28);
                                            jpsubo8 = (sub_base_o + 32);

                                            for(; jp <= jBj4; jp += 36) {
                                                B_k_j0 = B_k + jpsubo0; C_i_j0 = C_i + jpsubo0;
                                                B_k_j1 = B_k + jpsubo1; C_i_j1 = C_i + jpsubo1;
                                                B_k_j2 = B_k + jpsubo2; C_i_j2 = C_i + jpsubo2;
                                                B_k_j3 = B_k + jpsubo3; C_i_j3 = C_i + jpsubo3;
                                                B_k_j4 = B_k + jpsubo4; C_i_j4 = C_i + jpsubo4;
                                                B_k_j5 = B_k + jpsubo5; C_i_j5 = C_i + jpsubo5;
                                                B_k_j6 = B_k + jpsubo6; C_i_j6 = C_i + jpsubo6;
                                                B_k_j7 = B_k + jpsubo7; C_i_j7 = C_i + jpsubo7;
                                                B_k_j8 = B_k + jpsubo8; C_i_j8 = C_i + jpsubo8;
                                                b_v0 = _mm256_load_pd(B_k_j0);
                                                b_v1 = _mm256_load_pd(B_k_j1);
                                                b_v2 = _mm256_load_pd(B_k_j2);
                                                b_v3 = _mm256_load_pd(B_k_j3);
                                                b_v4 = _mm256_load_pd(B_k_j4);
                                                b_v5 = _mm256_load_pd(B_k_j5);
                                                b_v6 = _mm256_load_pd(B_k_j6);
                                                b_v7 = _mm256_load_pd(B_k_j7);
                                                b_v8 = _mm256_load_pd(B_k_j8);
                                                c_v0 = _mm256_load_pd(C_i_j0);
                                                c_v1 = _mm256_load_pd(C_i_j1);
                                                c_v2 = _mm256_load_pd(C_i_j2);
                                                c_v3 = _mm256_load_pd(C_i_j3);
                                                c_v4 = _mm256_load_pd(C_i_j4);
                                                c_v5 = _mm256_load_pd(C_i_j5);
                                                c_v6 = _mm256_load_pd(C_i_j6);
                                                c_v7 = _mm256_load_pd(C_i_j7);
                                                c_v8 = _mm256_load_pd(C_i_j8);
                                                apb_v0 = _mm256_add_pd(a_v, b_v0);
                                                apb_v1 = _mm256_add_pd(a_v, b_v1);
                                                apb_v2 = _mm256_add_pd(a_v, b_v2);
                                                apb_v3 = _mm256_add_pd(a_v, b_v3);
                                                apb_v4 = _mm256_add_pd(a_v, b_v4);
                                                apb_v5 = _mm256_add_pd(a_v, b_v5);
                                                apb_v6 = _mm256_add_pd(a_v, b_v6);
                                                apb_v7 = _mm256_add_pd(a_v, b_v7);
                                                apb_v8 = _mm256_add_pd(a_v, b_v8);
                                                res0 = _mm256_min_pd(c_v0, apb_v0);
                                                res1 = _mm256_min_pd(c_v1, apb_v1);
                                                res2 = _mm256_min_pd(c_v2, apb_v2);
                                                res3 = _mm256_min_pd(c_v3, apb_v3);
                                                res4 = _mm256_min_pd(c_v4, apb_v4);
                                                res5 = _mm256_min_pd(c_v5, apb_v5);
                                                res6 = _mm256_min_pd(c_v6, apb_v6);
                                                res7 = _mm256_min_pd(c_v7, apb_v7);
                                                res8 = _mm256_min_pd(c_v8, apb_v8);
                                                _mm256_store_pd(C_i_j0, res0);
                                                _mm256_store_pd(C_i_j1, res1);
                                                _mm256_store_pd(C_i_j2, res2);
                                                _mm256_store_pd(C_i_j3, res3);
                                                _mm256_store_pd(C_i_j4, res4);
                                                _mm256_store_pd(C_i_j5, res5);
                                                _mm256_store_pd(C_i_j6, res6);
                                                _mm256_store_pd(C_i_j7, res7);
                                                _mm256_store_pd(C_i_j8, res8);
                                                jpsubo0+=36;
                                                jpsubo1+=36;
                                                jpsubo2+=36;
                                                jpsubo3+=36;
                                                jpsubo4+=36;
                                                jpsubo5+=36;
                                                jpsubo6+=36;
                                                jpsubo7+=36;
                                                jpsubo8+=36;
                                            }

                                            for(; jp < j + Bj; ++jp) {
                                                C[ipsubmn + (jp + sub_base_o)] = min(
                                                    C[ipsubmn + (jp + sub_base_o)], 
                                                    A[ipsubmn + (kp + sub_base_l)] + B[((kp + sub_base_l) * n) + (jp + sub_base_o)]
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
    double r1 = benchmark_tiled_timed(n, fw_abc_min_plus, opt_tiled, L1, Bi, Bj, Bk);
#endif

    return 0;
}
    