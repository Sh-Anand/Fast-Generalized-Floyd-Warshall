
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

int main(int argc, char **argv) {
    if (argc!=3) {printf("usage: FW <n> <fw> (fw = 0,1,2 = (min,plus), (or,and), (max, min))\n"); return -1;}
    int n = atoi(argv[1]);
    int fwi = atoi(argv[2]);
    printf("n=%d \n",n);
    double *C = (double *)malloc(n*n*sizeof(double));
    fw_min_plus(C, n);

    free(C);

    return 0;
}