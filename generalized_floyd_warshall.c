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

///////////////////////////BLOCKED MIN PLUS///////////////////////////

void basic_blocked_fw_min_plus(double *C, int n, int Bi, int Bj, int Bk) {
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; i += Bi) {
            for (int j = 0; j < n; j += Bj) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    for(int jp = j; jp < j + Bj; ++jp) {
                        C[(ip * n) + jp] =  min(C[ip * n + jp], C[ip * n + k] + C[k * n + jp]);
                    }
                }
            }
        }
    }
}

void opt_blocked_fw_min_plus(double *C, int n, int Bi, int Bj, int Bk) {
    int ipn = 0;
    int kn = 0;
    int ipnplusjp = 0;
    double c_ipnpluskp = 0.0;
    double c_kpnplusjp = 0.0;
    double c_0pc_1 = 0.0;
    double min_c = 0.0;
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; i += Bi) {
            for (int j = 0; j < n; j += Bj) {
                kn = k * n;
                for(int ip = i; ip < i + Bi; ++ip) {
                    ipn = ip * n;
                    c_ipnpluskp = C[ipn + k];
                    for(int jp = j; jp < j + Bj; ++jp) {
                        ipnplusjp = ipn + jp;
                        c_kpnplusjp = C[kn + jp];
                        c_0pc_1 = c_ipnpluskp + c_kpnplusjp;
                        min_c =  min(C[ipnplusjp], c_0pc_1);
                        C[ipnplusjp] = min_c;
                    }
                }
            }
        }
    }
}

///////////////////////////TILED MIN PLUS///////////////////////////
//Used as a basis for the sumatrix fwis
void basic_fwi_min_plus(double* A, double* B, double *C, int n, int Bi, int Bj) {
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; i += Bi) {
            for (int j = 0; j < n; j += Bj) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    for(int jp = j; jp < j + Bj; ++jp) {
                        C[(ip * n) + jp] =  min(C[ip * n + jp], A[ip * n + k] + B[k * n + jp]);
                    }
                }
            }
        }
    }
}

void base_fwi_phase1_min_plus(double* A, double* B, double* C, int n, int l, int m, int L1, int Bi, int Bj) {
    int sub_base_l = l * L1;
    int sub_base_m = m * L1;
    for (int k = 0; k < L1; ++k) {
        for (int i = 0; i < L1; i += Bi) {
            for (int j = 0; j < L1; j += Bj) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    for(int jp = j; jp < j + Bj; ++jp) {
                        C[((ip + sub_base_l) * n) + (jp + sub_base_m)] = min(
                            C[((ip + sub_base_l) * n) + (jp + sub_base_m)], 
                            (A[((ip + sub_base_l) * n) + (k + sub_base_m)] + B[((k + sub_base_l) * n) + (jp + sub_base_m)])
                        );
                    }
                }
            }
        }
    }
}

void base_fwi_phase2_min_plus(double* A, double* B, double* C, int n, int l, int m, int L1, int Bi, int Bj) {
    int sub_base_l = l * L1;
    int sub_base_m = m * L1;
    for (int k = 0; k < L1; ++k) {
        for (int i = 0; i < L1; i += Bi) {
            for (int j = 0; j < L1; j += Bj) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    for(int jp = j; jp < j + Bj; ++jp) {
                        C[((ip + sub_base_l) * n) + (jp + sub_base_m)] = min(
                            C[((ip + sub_base_l) * n) + (jp + sub_base_m)], 
                            A[((ip + sub_base_l) * n) + (k + sub_base_l)] + 
                                B[((k + sub_base_l) * n) + (jp + sub_base_m)]
                        );
                    }
                }
            }
        }
    }
}

void base_fwi_phase3_min_plus(double* A, double* B, double* C, int n, int l, int m, int L1, int Bi, int Bj) {
    int sub_base_l = l * L1;
    int sub_base_m = m * L1;
    for (int k = 0; k < L1; ++k) {
        for (int i = 0; i < L1; i += Bi) {
            for (int j = 0; j < L1; j += Bj) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    for(int jp = j; jp < j + Bj; ++jp) {
                        C[((ip + sub_base_m) * n) + (jp + sub_base_l)] = min(
                            C[((ip + sub_base_m) * n) + (jp + sub_base_l)], 
                            A[((ip + sub_base_m) * n) + (k + sub_base_l)] + 
                                B[((k + sub_base_l) * n) + (jp + sub_base_l)]
                        );
                    }
                }
            }
        }
    }
}

/**
 * @brief Final phase of tilling process, here A, B & C are distinct, and we need 3 submatrix parameters
 * @param [A, B, C] FW matrices
 * @param n original size of the NxN matrices
 * @param l the first submatrix parameter (k in the paper)
 * @param m the second submatrix parameter (i in the paper)
 * @param o the last submatrix parameter (j in the paper)
 * @param L1 the tilling size
 * @param [Bi, Bj, Bk] tilling factors 
 */
void base_fwi_abc_min_plus(double* A, double* B, double* C, int n, int l, int m, int o, int L1, int Bi, int Bj, int Bk) {
    int sub_base_l = l * L1;
    int sub_base_m = m * L1;
    int sub_base_o = o * L1;
    for (int i = 0; i < L1; i += Bi) {
        for (int j = 0; j < L1; j += Bj) {
            for (int k = 0; k < L1; k += Bk) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    for(int jp = j; jp < j + Bj; ++jp) {
                        for(int kp = k; kp < k + Bk; ++kp) {
                            C[((ip + sub_base_m) * n) + (jp + sub_base_o)] = min(
                            C[((ip + sub_base_m) * n) + (jp + sub_base_o)], 
                            A[((ip + sub_base_m) * n) + (kp + sub_base_l)] + 
                                B[((kp + sub_base_l) * n) + (jp + sub_base_o)]
                            );
                        }
                    }
                }
            }
        }
    }
}

void fwi_phase1_min_plus(double* A, double* B, double* C, int n, int l, int m, int L1, int Bi, int Bj) {
    int sub_base_l = l * L1;
    int sub_base_m = m * L1;
    int ipln = 0;
    int jpm = 0;
    int iplnpjpm = 0;
    int kpbm = 0;
    int kpbln = 0;
    int iplnpkpbm = 0;
    double apb = 0.0;
    double min_c = 0.0;
    double c = 0.0;
    double a = 0.0;
    for (int k = 0; k < L1; ++k) {
        kpbm = k + sub_base_m;
        kpbln = ((k + sub_base_l) * n);
        for (int i = 0; i < L1; i += Bi) {
            for (int j = 0; j < L1; j += Bj) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    ipln = ((ip + sub_base_l) * n);
                    iplnpkpbm = ipln + kpbm;
                    a = A[iplnpkpbm];
                    for(int jp = j; jp < j + Bj; ++jp) {
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
}

void fwi_phase2_min_plus(double* A, double* B, double* C, int n, int l, int m, int L1, int Bi, int Bj) {
    int sub_base_l = l * L1;
    int sub_base_m = m * L1;
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
    for (int k = 0; k < L1; ++k) {
        kl = (k + sub_base_l);
        kln = kl * n;
        for (int i = 0; i < L1; i += Bi) {
            for (int j = 0; j < L1; j += Bj) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    ipln = ((ip + sub_base_l) * n);
                    iplnkl = ipln + kl;
                    a = A[iplnkl];
                    for(int jp = j; jp < j + Bj; ++jp) {
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

void fwi_phase3_min_plus(double* A, double* B, double* C, int n, int l, int m, int L1, int Bi, int Bj) {
    int sub_base_l = l * L1;
    int sub_base_m = m * L1;
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
    for (int k = 0; k < L1; ++k) {
        kl = (k + sub_base_l);
        kln = kl * n;
        for (int i = 0; i < L1; i += Bi) {
            for (int j = 0; j < L1; j += Bj) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    ipmn = ((ip + sub_base_m) * n);
                    ipmnkl = ipmn + kl; 
                    a = A[ipmnkl];
                    for(int jp = j; jp < j + Bj; ++jp) {
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

/**
 * @brief Final phase of tilling process, here A, B & C are distinct, and we need 3 submatrix parameters
 * @param [A, B, C] FW matrices
 * @param n original size of the NxN matrices
 * @param l the first submatrix parameter (k in the paper)
 * @param m the second submatrix parameter (i in the paper)
 * @param o the last submatrix parameter (j in the paper)
 * @param L1 the tilling size
 * @param [Bi, Bj, Bk] tilling factors 
 */
void fwi_abc_min_plus(double* A, double* B, double* C, int n, int l, int m, int o, int L1, int Bi, int Bj, int Bk) {
    int sub_base_l = l * L1;
    int sub_base_m = m * L1;
    int sub_base_o = o * L1;
    int ipmn = 0;
    int jpo = 0;
    int ipmnpjpo = 0;
    int kpl = 0;
    int kplnpjpo = 0;
    double apb = 0.0;
    double c = 0.0;
    double min_c = 0.0;
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

/**
 * @brief FW iteration for MIN_PLUS
 * @param A, the first operand matrix
 * @param B, the second operand matrix
 * @param C, the result matrix 
 * @param n, the size of the NxN matrices
 * @param Bi, tilling factor over i
 * @param Bj, tilling factor over j
 */
void fwi_min_plus(double* A, double* B, double* C, int n, int Bi, int Bj) {
    int ipn = 0;
    int kn = 0;
    int ipnplusjp = 0;
    double a = 0.0;
    double b = 0.0;
    double apb = 0.0;
    double min_c = 0.0;
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; i += Bi) {
            for (int j = 0; j < n; j += Bj) {
                kn = k * n;
                for(int ip = i; ip < i + Bi; ++ip) {
                    ipn = ip * n;
                    a = A[ipn + k];
                    for(int jp = j; jp < j + Bj; ++jp) {
                        ipnplusjp = ipn + jp;
                        b = B[kn + jp];
                        apb = a + b;
                        min_c =  min(C[ipnplusjp], apb);
                        C[ipnplusjp] = min_c;
                    }
                }
            }
        }
    }
}

/**
 * @brief FW iteration for MIN_PLUS using three distinct matrices A, B, C
 * @param A, the first operand matrix
 * @param B, the second operand matrix
 * @param C, the result matrix 
 * @param n, the size of the NxN matrices
 * @param Bi, tilling factor over i
 * @param Bj, tilling factor over j
 * @param Bk, tilling factor over k
 */
void fwi_abc(double* A, double* B, double* C, int n, int Bi, int Bj, int Bk) {
    int ipn = 0;
    int ipnpjp = 0;
    int ipnpk = 0;
    int kpnpjp = 0;
    double apb = 0.0;
    double min_c = 0.0;
    for (int i = 0; i < n; i += Bi) {
        for (int j = 0; j < n; j += Bj) {
            for (int k = 0; k < n; k += Bk) {
                for(int ip = i; ip < i + Bi; ++ip) {
                    ipn = ip * n;
                    for(int jp = j; jp < j + Bj; ++jp) {
                        ipnpjp = ipn + jp;
                        for(int kp = k; kp < k + Bk; ++kp) {
                            ipnpk = ipn + k;
                            kpnpjp = (kp * n) + jp;
                            apb = A[ipnpk] + B[kpnpjp];
                            min_c = min(C[ipnpjp], apb);
                            C[ipnpjp] = min_c;
                        }
                    }
                }
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
void tiled_fw_min_plus(double* A, double* B, double* C, int L1, int n, int Bi, int Bj, int Bk) {
    int m = n / L1;
    for(int k = 0; k < m; ++k) {
        //Tilling phase 1 (update C_kk)
        fwi_phase1_min_plus(A, B, C, n, k, k, L1, Bi, Bj);

        //Tilling phase 2 (Update all tiles in same row as C_kk)
        for(int j = 0; j < m; ++j) {
            if(j != k) {
                fwi_phase2_min_plus(A, B, C, n, k, j, L1, Bi, Bj);
            }
        }

        //Tilling phase 3 (Update all tiles in same column as C_kk)
        for(int i = 0; i < m; ++i) {
            if(i != k) {
                fwi_phase3_min_plus(A, B, C, n, k, i, L1, Bi, Bj);
            }
        }

        //Tilling phase 4 (Update all remaining tiles)
        for(int i = 0; i < m; ++i) {
            if(i != k){
                for(int j = 0; j < m; ++j) {
                    if(j != k) {
                        fwi_abc_min_plus(A, B, C, n, k, i, j, L1, Bi, Bj, Bk);
                    }
                }
            }
        }
    }
}

//////////////////////VECT MIN PLUS//////////////////////

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

//////////////////////MAX MIN//////////////////////

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
                C[i*n + j] = (u_int64_t)C[i*n + j] | ((u_int64_t)C[i*n + k] & (u_int64_t)C[k*n + j]);
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

void opt_fw_or_and_256(u_int64_t *C,int n) {
    assert(C == (u_int64_t *)(__m256i*)C);
    u_int64_t *c_addr = C, *addr_ik, *addr_ij , *addr_kj;
    __m256i c_ij, c_ik, c_kj, c2, cmp_lt, res;
    for (size_t k = 0; k < n; k++) {
        int i_n = 0;
        for (size_t i = 0; i < n; i++) {
            size_t j = 0;
            addr_ik = c_addr + (i_n + k);
            //printf("access ik");
            
            c_ik = _mm256_set_epi64x (*addr_ik,
                                        *addr_ik,
                                        *addr_ik,
                                        *addr_ik);
            //printf("%04llx\n", (unsigned long long)*addr_ik);
            /* printf("cik : ");
            print_vector(c_ik); */
            for (; j <= n - 4; j+=4) {
                addr_ij = c_addr + (i_n + j);
                addr_kj = c_addr + (k*n + j);
                //printf("access ij");
                c_ij = _mm256_load_si256((__m256i *)addr_ij);
                
                /* printf("cij : ");
                print_vector(c_ij); */
                //printf("access kj");
                c_kj = _mm256_load_si256((__m256i *)addr_kj);
 /*                
                printf("ckj : ");
                print_vector(c_kj);
       */          
                c2 = _mm256_and_si256(c_ik, c_kj);

               /*  printf("and : ");
                print_vector(c2); */

                res = _mm256_or_si256(c_ij, c2);

                /* printf("or : ");
                print_vector(res); */
                //printf("store ij");
                _mm256_store_si256((__m256i *) addr_ij ,res);
            }
            
            for(;j < n; j++){
                C[i_n + j] = C[i_n + j] | ( C[i_n + k] &  C[k*n + j]); 
                printf("%d %d %d\n", (int)(i_n+j), (int)(i_n+k),(int)(k*n + j));
            }
            
            i_n += n;
        }
    }
}


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

void test_or_and(int n){
    u_int64_t *C_base = (u_int64_t *)malloc(n*n*sizeof(double));
    u_int64_t *C_opt = (u_int64_t *)malloc(n*n*sizeof(double));
    init_bit_matrices(C_base, C_opt, n);
    
     for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%lu ", C_base[n*i + j]);
            printf("%lu\n", C_opt[n*i + j]);
        }
    }

    // Run baseline function on C
    fw_or_and_int(C_base, n);
    // Run optimized function on C
    opt_fw_or_and_256(C_opt, n);
    /*
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            //printf("base[%d][%d] = %lf ", i, j, C_base[n*i + j]);
            //printf("opt[%d][%d] = %lf\n", i, j, C_opt[n*i + j]);
            print_bits(C_base[n*i + j]);
            print_bits(C_opt[n*i + j]);
        }
    }
    */
    // Compare both
    for(int i = 0; i < n; ++i) {
        assert(C_opt[i] == C_base[i]);
    }

    free(C_base);
    free(C_opt);
}

void test_blocked(int n, void (*baseline)(double*, int), void (*optimization)(double*, int, int, int, int)) {
    double *C_base = (double *)malloc(n*n*sizeof(double));
    double *C_opt = (double *)malloc(n*n*sizeof(double));
    init_matrices(C_base, C_opt, n);
    // Run baseline function on C
    baseline(C_base, n);

    int Bi, Bj, Bk;
    Bi = n/4;
    Bj = n/4;
    Bk = n/4;
    // Run optimized function on C
    optimization(C_opt, n, Bi, Bj, Bk);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("base[%d][%d] = %lf ", i, j, C_base[n*i + j]);
            printf("opt[%d][%d] = %lf\n", i, j, C_opt[n*i + j]);
        }
    }
    
    // Compare both
    for(int i = 0; i < n*n; ++i) {
        assert(C_opt[i] == C_base[i]);
    }

    free(C_base);
    free(C_opt);
}

void test_tiled(int n, void (*baseline)(double*, double*, double*, int), 
    void (*optimization)(double*, double*, double*, int, int, int, int, int)
) {
    double *A_base = (double *)malloc(n*n*sizeof(double));
    double *A_opt = (double *)malloc(n*n*sizeof(double));
    double *B_base = (double *)malloc(n*n*sizeof(double));
    double *B_opt = (double *)malloc(n*n*sizeof(double));
    double *C_base = (double *)malloc(n*n*sizeof(double));
    double *C_opt = (double *)malloc(n*n*sizeof(double));
    init_matrices(A_base, A_opt, n);
    init_matrices(B_base, B_opt, n);
    init_matrices(C_base, C_opt, n);
    // Run baseline function on C
    baseline(A_base, B_base, C_base, n);

    int L1 = n/4;

    int Bi, Bj, Bk;
    Bi = L1/4;
    Bj = L1/4;
    Bk = L1/4;

    printf("L1 = %d\n", L1);
    
    // Run optimized function on C
    optimization(A_opt, B_opt, C_opt, L1, n, Bi, Bj, Bk);
    /*
    printf("\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("base[%d][%d] = %lf ", i, j, C_base[n*i + j]);
            printf("opt[%d][%d] = %lf\n", i, j, C_opt[n*i + j]);
        }
    }
    printf("\n");
    */
    // Compare both 
    for(int i = 0; i < n*n; ++i) {
        assert(C_opt[i] == C_base[i]);
    }

    free(A_base);
    free(A_opt);
    free(B_base);
    free(B_opt);
    free(C_base);
    free(C_opt);
}

void test(int n, void (*baseline)(double*, int), void (*optimization)(double*, int)) {
    double *C_base = (double *)malloc(n*n*sizeof(double));
    double *C_opt = (double *)malloc(n*n*sizeof(double));
        init_matrices(C_base, C_opt, n);
    
/*     for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            printf("%lf", C_base[n*i + j]);
            printf("%lf", C_opt[n*i + j]);
        }
    }
 */
    // Run baseline function on C
    baseline(C_base, n);
    // Run optimized function on C
    optimization(C_opt, n);
    /*
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            //printf("base[%d][%d] = %lf ", i, j, C_base[n*i + j]);
            //printf("opt[%d][%d] = %lf\n", i, j, C_opt[n*i + j]);
            print_bits(C_base[n*i + j]);
            print_bits(C_opt[n*i + j]);
        }
    }
    */
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
    /*
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            printf("%lf", C[n*i + j]);
    */
    //printf("\n");
    //test(n, fw_or_and, opt_fw_or_and_256);
    //test_blocked(n, fw_min_plus, opt_blocked_fw_min_plus);
    //test_tiled(n, fw_abc_min_plus, tiled_fw_min_plus);
    //test_or_and(n);

    free(C);

    return 0;
}
