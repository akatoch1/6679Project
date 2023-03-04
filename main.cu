#include <stdio.h>
#include "kernel1.cu"
#include "kernel2.cu"
#include "kernel3.cu"
#include "kernel4.cu"

int main(int argc, char**argv) {
    unsigned int m; // number of constraints
    unsigned int n; // number of variables (no slack vars)
    unsigned int b; // b = highest value of constraint value 
   
    m = atoi(argv[0]);
    n = atoi(argv[1]);
    b = atoi(argv[2]);
    
    // Initialize host variables
    double* xB_h = (double*) malloc( sizeof(double)*m );
    for (unsigned int i=0; i < m; i++) { xB_h[i] = n+1+i; }

    double* xN_h = (double*) malloc( sizeof(double)*n );
    for (unsigned int i=0; i < n; i++) { xB_h[i] = i+1; }

    double* cB_h = (double*) malloc( sizeof(double)*m );
    for (unsigned int i=0; i < m; i++) { cB_h[i] = 0; }

    double* cN_h = (double*) malloc( sizeof(double)*n );
    for (unsigned int i=0; i < n; i++) { cN_h[i] = (rand()%10-5); }

    double* b_h = (double*) malloc( sizeof(double)*m );
    for (unsigned int i=0; i < m; i++) { b_h[i] = rand()%b; }

    double* B_h = (double*) malloc( sizeof(double)*m*m );
    for (unsigned int i=0; i < m; i++) { 
        for (unsigned int j=0; j < m; j++) { 
            if (i==j) {
                B_h[i][j] = 1;
            } else {
                B_h[i][j] = 0;
            }
        }
    }

    double* N_h = (double*) malloc( sizeof(double)*m*n );
    for (unsigned int i=0; i < m; i++) { 
        for (unsigned int j=0; j < n; j++) { 
            B_h[i][j] = (rand()%10-5);
        }
    }

    double* tab_h = (double*) malloc( sizeof(double)*(m+1)*(n+1));
    // Assign vals for 0,0
    // Assign vals for m,0
    // Assign vals for 0,n
    // Assign vals for m,n

    double* objLine_h = (double*) malloc( sizeof(double)*n);

    int* k_h = malloc(sizeof(int));
    
    int* r_h = malloc(sizeof(int));

    double* theta_h = (double*) malloc( sizeof(double)*(m+1));

    // Allocate device variables
    double* tab_d;
    cuda_ret = cudaMalloc((void**) &tab_d, sizeof(double)*(m+1)*(n+1));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");
    
    double* xB_d;
    cuda_ret = cudaMalloc((void**) &xB_d, sizeof(double)*m);
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    double* objLine_d;
    cuda_ret = cudaMalloc((void**) &objLine_d, sizeof(double)*n);
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    int* k_d;
    cuda_ret = cudaMalloc((void**) &k_d, sizeof(int));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    int* r_d;
    cuda_ret = cudaMalloc((void**) &r_d, sizeof(int));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    double* columnk_d;
    cuda_ret = cudaMalloc((void**) &columnk_d, sizeof(double)*(m+1));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    double* theta_d;
    cuda_ret = cudaMalloc((void**) &columnk_d, sizeof(double)*(m+1));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    // Copy tableau to Device
    cuda_ret = cudaMemcpy(s_d, s_h, sizeof(double)*(m+1)*(n+1), cudaMemcpyHostToDevice);
	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");

    cuda_ret = cudaMemcpy(xB_d, xB_h, sizeof(double)*m, cudaMemcpyHostToDevice);
	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");

    // Iterate until optimal is found

        // Copy first line of the Simplex tableau to host
        for (unsigned int i=0; i < m; i++) {objLine_d[i]= s_d[0][i+1] ;}

        cuda_ret = cudaMemcpy(objLine_h, objLine_d, sizeof(double)*n, cudaMemcpyDeviceToHost);
    	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to host");
        
        cudaDeviceSynchronize();

        // Find index of entering variable k
        double minValue = 1.0;
        int minIndex = 0;
        for (unsigned int i=0; i < m; i++) {
            if (objLine_h[i] < minValue) {
                minValue = objLine_h;
                minIndex = i;
            }
        }
        if (minValue>=0) {
            break;
        } else {
            k_h = minIndex + 1;
        }

        // Copy index k to device
        cuda_ret = cudaMemcpy(k_d, k_h, sizeof(int), cudaMemcpyHostToDevice);
	        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");

        cudaDeviceSynchronize();

        // Kernel 1 to process ratio column
        const unsigned int THREADS_PER_BLOCK = 512;
        const unsigned int numBlocks1 = m/THREADS_PER_BLOCK + 1;
        dim3 gridDim(numBlocks1, 1, 1), blockDim(THREADS_PER_BLOCK, 1, 1);
        kernel1<<<gridDim, blockDim>>>(tab_d, theta_d, columnk_d, k_d);

        cudaDeviceSynchronize();

        // Copy ratio column to host
        cuda_ret = cudaMemcpy(theta_h, theta_d, sizeof(double)*(m+1), cudaMemcpyDeviceToHost);
    	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to host");
        
        // Find the index of the leaving variable r

        // Copy index r to device

        // Kernel 2 to update the line r of the Simplex tableau
        const unsigned int numBlocks2 = n/THREADS_PER_BLOCK + 1;
        dim3 gridDim(numBlocks2, 1, 1), blockDim(THREADS_PER_BLOCK, 1, 1);
        kernel2<<<gridDim, blockDim>>>(tab_d, columnk_d, k_d, r_d);

        cudaDeviceSynchronize();

        // Kernel 3 to update Simplex tableau
        const unsigned int numBlocksX3 = m/THREADS_PER_BLOCK + 1;
        const unsigned int numBlocksY3 = n/THREADS_PER_BLOCK + 1;
        dim3 gridDim(numBlocksX3, numBlocksY3, 1), blockDim(THREADS_PER_BLOCK, THREADS_PER_BLOCK, 1);
        kernel3<<<gridDim, blockDim>>>(tab_d, columnk_d, k_d, r_d);

        cudaDeviceSynchronize();

        // Kernel 4 to Update column k of the Simplex Tableau
        const unsigned int numBlocks4 = n/THREADS_PER_BLOCK + 1;
        dim3 gridDim(numBlocks4, 1, 1), blockDim(THREADS_PER_BLOCK, 1, 1);
        kernel4<<<gridDim, blockDim>>>(tab_d, columnk_d, k_d, r_d);

        cudaDeviceSynchronize();

    // Calculate optimal value and return it

}
