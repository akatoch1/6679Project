#include <stdio.h>
#include "kernel1.cu"
#include "kernel2.cu"
#include "kernel3.cu"
#include "kernel4.cu"
#include "support.h"
void printMatrix(double *m, int r, int c) {
  for (int i = 0; i < r; i++)
    {
      for (int j = 0; j < c; j++)
      {
        printf("%f ", m[i*c + j]);
	}	    
      printf("\n");
    }
}
void printMatrix2(double **m, int r, int c) {
  for (int i = 0; i < r; i++)
    {
      for (int j = 0; j < c; j++)
      {
        printf("%f ", m[i][j]);
        }
      printf("\n");
    }
}

int main(int argc, char**argv) {


    unsigned int m; // number of constraints
    unsigned int n; // number of variables (no slack vars)
     // b = highest value of constraint value 
    unsigned int b;
    //unsigned int hb;

    //m = atoi(argv[1]);
    //n = atoi(argv[2]);
    //b = atoi(argv[3]);
    m = 2;
    n = 2;
    
    cudaError_t cuda_ret;
    // Initialize host variables
    double* xB_h = (double*) malloc( sizeof(double)*m ); // index's of the variables in the Basis.
    for (unsigned int i=0; i < m; i++) { xB_h[i] = n+i; }
    
    double* xN_h = (double*) malloc( sizeof(double)*n ); // index's of the variables not in the Basis.
    for (unsigned int i=0; i < n; i++) { xN_h[i] = i; }
    
    double* cB_h = (double*) malloc( sizeof(double)*m ); // C values for the basic variables
    for (unsigned int i=0; i < m; i++) { cB_h[i] = 0; }
    
    
    double* cN_h = (double*) malloc( sizeof(double)*n ); // C values for the non-basic variables
    for (unsigned int i=0; i < n; i++) { cN_h[i] = (rand()%10-5); }
    cN_h[0] = 40;
    cN_h[1] = 30;
    double* b_h = (double*) malloc( sizeof(double)*m ); // Right hand side values 
    for (unsigned int i=0; i < m; i++) { b_h[i] = rand()%b; }
    b_h[0] = 12;
    b_h[1] = 16;
    double* svec = (double*) malloc( sizeof(double)*(m+1));
    



    //double* B_h = (double*) malloc( sizeof(double)*m*m ); // Constraint coefficents of basic variables
    double* B_h[m];
      for (int i = 0; i < m; i++) {
        B_h[i] = (double*)malloc(m * sizeof(double));
      }  
    for (unsigned int i=0; i < m; i++) { 
        for (unsigned int j=0; j < m; j++) { 
            if (i==j) {
                B_h[i][j] = 1;
            } else {
                B_h[i][j] = 0;
            }
        }
    }


    double* cBB	= (double*) malloc( sizeof(double)*m);
    for	(int i = 0; i < m; i++) {
    	double z = 0;
    	for (int j = 0; j < m; j++) {
            z += cB_h[j] * B_h[j][i];
     	}
	cBB[i] = z;
    }
    double z = 0;
    for (int i = 0; i < m; i++) {
    	z += cBB[i] * b_h[i];
    }
    svec[0] = z;
    for (int i = 1; i < m+1; i++) {
    	double sum = 0;
	for (int j = 0; j < m; j++) {
	    sum += B_h[i-1][j] * b_h[j];
	}
	svec[i] = sum;
    }   


       
    //double* N_h = (double*) malloc( sizeof(double)*m*n ); // Constraint coefficents of non-basic variables
    double* N_h[m];
    for (int i = 0; i < m; i++) {
        N_h[i] = (double*)malloc(n * sizeof(double));
    }
   
    fflush(stdout);
    for (unsigned int i=0; i < m; i++) { 
        for (unsigned int j=0; j < n; j++) { 
            N_h[i][j] = (rand()%10-5);
        }
    }
    N_h[0][0] = 1;
    N_h[0][1] = 1;
    N_h[1][0] = 2;
    N_h[1][1] = 1;
    
    
    double* tab_h = (double *) malloc(((n+1) * (m+1)) * sizeof(double));
    double* columnk_h = (double *) malloc((m+1) * sizeof(double));
 
 
    // Assign vals for 0,0
    // Assign vals for m,0
    // Assign vals for 0,n
    // Assign vals for m,n
    double* smat[m+1];
    for (int i = 0; i < m+1; i++) {
        smat[i] = (double*) malloc(n * sizeof(double));
    }
    //update first row of smat
    for (int i = 0; i < n; i++) {
        double val = 0;
        for (int j = 0; j < m; j++) {
            val += cBB[j] * N_h[j][i];
        }
        smat[0][i] = val - cN_h[i];
    }
    for (int i = 1; i < m+1; i++) {

        for (int j = 0; j < n; j++) {
	    double val = 0;
	    for (int k = 0; k < m; k++) {
	    	
		val += B_h[i-1][k] * N_h[k][j];
		
	    	
	    }
	    smat[i][j] = val;
        }
    }

    //update tab
  //  for (int i = 0; i < m+1; i++) {
   // 	tab_h[i * (n+1)] = svec[i];
  //  }    
    //printMatrix2(smat, m+1, n);
    for (int i = 0; i < m+1; i++) {
    	for (int j = 0; j < n+1; j++) {
	    tab_h[i * (n+1) + (j+1)] = smat[i][j];
	}   
    }	
    for (int i = 0; i < m+1; i++) {
        tab_h[i * (n+1)] = svec[i];
    }
    double* objLine_h = (double*) malloc( sizeof(double)*(n+1)); // Objective line used to determine k
    
    int k_h;
    


    double* theta_h = (double*) malloc( sizeof(double)*(m+1)); // Ratio of right-hand side to k row

    // Allocate device variables
    double* tab_d;
    cuda_ret = cudaMalloc((void**) &tab_d, sizeof(double)*(m+1)*(n+1));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");
    

 
    double* xB_d;
    cuda_ret = cudaMalloc((void**) &xB_d, sizeof(double)*m);
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");
    

	
    int k_d;
    //cuda_ret = cudaMalloc((void**) &k_d, sizeof(int));
//	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    int r_d;
    cuda_ret = cudaMalloc((void**) &r_d, sizeof(int));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    double* columnk_d;
    cuda_ret = cudaMalloc((void**) &columnk_d, sizeof(double)*(m+1));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    double* theta_d;
    cuda_ret = cudaMalloc((void**) &theta_d, sizeof(double)*(m+1));
	if(cuda_ret != cudaSuccess) FATAL("Unable to allocate device memory");

    
    // Copy tableau to Device

    cuda_ret = cudaMemcpy(theta_d, theta_h, sizeof(double)*(m+1), cudaMemcpyHostToDevice);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");

    cuda_ret = cudaMemcpy(tab_d, tab_h, sizeof(double)*(m+1)*(n+1), cudaMemcpyHostToDevice);
	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");

    cuda_ret = cudaMemcpy(xB_d, xB_h, sizeof(double)*m, cudaMemcpyHostToDevice);

	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");
   
    //printMatrix(tab_h, m+1, n+1);	
    // Iterate until optimal or infeasible/unbounded is found
    bool continueVar = true;
    const int THREADS_PER_BLOCK = 16;    
    while (continueVar == true) {
        // Copy first line of the Simplex tableau to host
        for (unsigned int i=0; i < n+1; i++) {objLine_h[i]= tab_h[i];}


        
        cudaDeviceSynchronize();

        // Find index of entering variable k
        double minValue = 1.0;
        int minIndex = 0;
        for (unsigned int i=1; i < n+1; i++) {
            if (objLine_h[i] < minValue) {
                minValue = objLine_h[i];
                minIndex = i;
            }
        }
        if (minValue>=0) {
            printf("Optimal Value");
            continueVar = false;
            break;
        } else {
            k_h = minIndex;
        }





        cudaDeviceSynchronize();
	printMatrix(tab_h, m+1, n+1);
	printf("\n");
        // Kernel 1 to process ratio column

        const unsigned int numBlocks1 = m/THREADS_PER_BLOCK + 1;
        dim3 gridDim(numBlocks1, 1, 1), blockDim(THREADS_PER_BLOCK, 1, 1);
        kernel1<<<gridDim, blockDim>>>(tab_d, theta_d, columnk_d, k_h, n);

        cudaDeviceSynchronize();
	
        // Copy ratio column to host
        cuda_ret = cudaMemcpy(theta_h, theta_d, sizeof(double)*(m+1), cudaMemcpyDeviceToHost);
    	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to host");
	
	cuda_ret = cudaMemcpy(columnk_h, columnk_d, sizeof(double)*(m+1), cudaMemcpyDeviceToHost);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to host");        

	cuda_ret = cudaMemcpy(columnk_d, columnk_h, sizeof(double)*(m+1), cudaMemcpyHostToDevice);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to host");
	// Find the index of the leaving variable r
        minValue = 10000.0;
        minIndex = 0;
	int r_h;
        for (unsigned int i=1; i < m+1; i++) {
            if (theta_h[i] < minValue && theta_h[i] > 0) {
                minValue = theta_h[i];
                minIndex = i;
            }
        }
        if (minValue==10000.0) {
            printf("unbounded");
            continueVar = false;
            break;
        } else {
            r_h = minIndex;
        }
        // Copy index r to device

	//printf("%d\n", r_h);
        // Kernel 2 to update the line r of the Simplex tableau
        const unsigned int numBlocks2 = n/THREADS_PER_BLOCK + 1;
        dim3 gridDim2(numBlocks2, 1, 1);
	dim3 blockDim2(THREADS_PER_BLOCK, 1, 1);

        kernel2<<<gridDim2, blockDim2>>>(tab_d, columnk_d, k_h, r_h, n);
	
        cudaDeviceSynchronize();
	
	cuda_ret = cudaMemcpy(tab_h, tab_d, sizeof(double)*(m+1)*(n+1), cudaMemcpyDeviceToHost);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");
	
	printMatrix(tab_h, m+1,	 n+1);
	printf("\n");
	cuda_ret = cudaMemcpy(tab_d, tab_h, sizeof(double)*(m+1)*(n+1), cudaMemcpyHostToDevice);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");
	
	cuda_ret = cudaMemcpy(columnk_h, columnk_d, sizeof(double)*(m+1), cudaMemcpyDeviceToHost);
	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");	
	
	cuda_ret = cudaMemcpy(columnk_d, columnk_h, sizeof(double)*(m+1), cudaMemcpyHostToDevice);
	if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");	
	
        // Kernel 3 to update Simplex tableau
        const unsigned int numBlocksX3 = m/THREADS_PER_BLOCK + 1;
        const unsigned int numBlocksY3 = n/THREADS_PER_BLOCK + 1;
        dim3 gridDim3(numBlocksX3, numBlocksY3), blockDim3(THREADS_PER_BLOCK, THREADS_PER_BLOCK);
	
        kernel3<<<gridDim3, blockDim3>>>(tab_d, columnk_d, k_h, r_h, n, m);
	printf("\n");
        cudaDeviceSynchronize();
	cuda_ret = cudaMemcpy(tab_h, tab_d, sizeof(double)*(m+1)*(n+1), cudaMemcpyDeviceToHost);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");



	printMatrix(tab_h, m+1,n+1);
	printf("\n");
        // Kernel 4 to Update column k of the Simplex Tableau
	cuda_ret = cudaMemcpy(tab_d, tab_h, sizeof(double)*(m+1)*(n+1), cudaMemcpyHostToDevice);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");	

	cuda_ret = cudaMemcpy(columnk_h, columnk_d, sizeof(double)*(m+1), cudaMemcpyDeviceToHost);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");

//	printMatrix(columnk_h, 1, m+1);
//	printf("%d\n", k_h);
//	printf("%d\n", r_h);

        cuda_ret = cudaMemcpy(columnk_d, columnk_h, sizeof(double)*(m+1), cudaMemcpyHostToDevice);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");
	
        const unsigned int numBlocks4 = n/THREADS_PER_BLOCK + 1;
        dim3 gridDim4(numBlocks4, 1), blockDim4(THREADS_PER_BLOCK, 1);
        kernel4<<<gridDim4, blockDim4>>>(tab_d, columnk_d, k_h, r_h, n);
		
        cudaDeviceSynchronize();
    	cuda_ret = cudaMemcpy(tab_h, tab_d, sizeof(double)*(m+1)*(n+1), cudaMemcpyDeviceToHost);
        if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to device");
       	printMatrix(tab_h, m+1, n+1);
	printf("############################################\n");
    }
    // Calculate optimal value and return it
    //cuda_ret = cudaMemcpy(objLine_h, objLine_d, sizeof(double)*n, cudaMemcpyDeviceToHost);
    //if(cuda_ret != cudaSuccess) FATAL("Unable to copy memory to host");
    
}
