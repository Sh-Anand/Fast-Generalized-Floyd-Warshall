#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

void initialize_matrix(double* adj, int n){
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            double r = ((double )rand() + 1);
            //printf("%f\n", r/RAND_MAX);
            adj[i * n + j] = r/RAND_MAX;

        }
        
    }
    
}

int main(){
    int n = 10;
    double* adj = (double *)malloc(n*n*sizeof(double));
    initialize_matrix(adj, n);
    for (size_t k = 0; k < n; k++){
        for (size_t i = 0; i < n; i++){
            for (size_t j = 0; j < n; j++){
                adj[i*n + j] = max(adj[i*n + j], min(adj[i*n + j], adj[i*n + j]));
            }
            
        }
        
    }
    
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < n; j++){
            printf("%f ",adj[i*n + j]);
        }
        printf("\n");
    }
    
}

