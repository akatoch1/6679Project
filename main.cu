#include <stdio.h>

int main(int argc, char**argv) {
    // m = number of constraints
    // n = number of variables (no slack vars)
    // b = highest value of constraint value 
    
    // Initialize host variables
    double* xB_h = (double*) malloc( sizeof(double)*m );
    for (unsigned int i=0; i < m; i++) { xB_h[i] = n+1+i; }

    double* xN_h = (double*) malloc( sizeof(double)*n );
    for (unsigned int i=0; i < n; i++) { xB_h[i] = i+1; }

    double* cB_h = (double*) malloc( sizeof(double)*m );
    for (unsigned int i=0; i < m; i++) { cB_h[i] = 0; }

    float* cN_h = (float*) malloc( sizeof(float)*n );
    for (unsigned int i=0; i < n; i++) { cN_h[i] = (rand()%10-5); }

    float* b_h = (float*) malloc( sizeof(float)*m );
    for (unsigned int i=0; i < m; i++) { b_h[i] = rand()%b; }

    float* B_h = (float*) malloc( sizeof(float)*m*m );
    for (unsigned int i=0; i < m; i++) { 
        for (unsigned int j=0; j < m; j++) { 
            if (i==j) {
                B_h[i][j] = 1
            } else {
                B_h[i][j] = 0
            }
        }
    }

    float* N_h = (float*) malloc( sizeof(float)*m*n );
    for (unsigned int i=0; i < m; i++) { 
        for (unsigned int j=0; j < n; j++) { 
            B_h[i][j] = (rand()%10-5);
        }
    }

    float* s_h = (float*) malloc( sizeof(float)*(m+1)*(n+1));
    // Assign vals for 0,0
    // Assign vals for m,0
    // Assign vals for 0,n
    // Assign vals for m,n

    // Allocate device variables
    float* s_d;
    cuda_ret = cudaMalloc((void**) &s_d, sizeof(float)*(m+1)*(n+1));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");
    
    float* xB_d;
    cuda_ret = cudaMalloc((void**) &xB_d, sizeof(float)*m);
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    // Copy tableau to Device
    cuda_ret = cudaMemcpy(s_d, s_h, sizeof(float)*(m+1)*(n+1), cudaMemcpyHostToDevice);
	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");

    cuda_ret = cudaMemcpy(xB_d, xB_h, sizeof(float)*m, cudaMemcpyHostToDevice);
	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");

    // Iterate until optimal is found

        // Copy first line of the Simplex tableau to host

        // Find index of entering variable k

        // Copy index k to device

        // Kernel 1 to process ratio column

        // Copy ratio column to host

        // Find the index of the leaving variable r

        // Copy index r to device

        // Kernel 2 to update the line r of the Simplex tableau

        // Kernel 3 to update Simplex tableau

        // Kernel 4 to Update column k of the Simplex Tableau

        // Check if should stop

    // Calculate optimal value and return it

}
