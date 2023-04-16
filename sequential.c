#include <stdio.h>
#include <stdbool.h>
#include <string.h>

int main(int argc, char**argv) {
    int m = 2;
    int n = 2;
    int c[2] = {40,30}; //
    int A[2][2]; //
    A[0][0] = 1; //
    A[0][1] = 1; //
    A[1][0] = 2; //
    A[1][1] = 1; //
    int b[2] = {12,16}; //

    double xB_h[m];
    for (unsigned int i=0; i < m; i++) { xB_h[i] = n+i;}
    
    double xN_h[n];
    for (unsigned int i=0; i < n; i++) { xN_h[i] = i;}
  
    double cB_h[m];
    for (unsigned int i=0; i < m; i++) { cB_h[i] = 0;}
   
    double cN_h[n];
    for (unsigned int i=0; i < n; i++) { cN_h[i] = c[i];}
 
    double b_h[m];
    for (unsigned int i=0; i < m; i++) { b_h[i] = b[i];}

    double B_h[m][m];
    for (unsigned int i=0; i < m; i++) { 
        for (unsigned int j=0; j < m; j++) { 
            if (i==j) {
                B_h[i][j] = 1;
            } else {
                B_h[i][j] = 0;
            }
        }
    }

    double N_h[m][n];
    for (unsigned int i=0; i < m; i++) { 
        for (unsigned int j=0; j < m; j++) { 
            N_h[i][j] = A[i][j];
        }
    }

    double svec[m+1];
    double cBB[m];
    double z;
    for	(int i = 0; i < m; i++) {
    	z = 0;
    	for (int j = 0; j < m; j++) {
            z += cB_h[j] * B_h[j][i];
     	}
	    cBB[i] = z;
    }
    z = 0;
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

    
    double tab[m+1][n+1];

    double smat[m+1][m+1];
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
    for (int i = 0; i < m+1; i++) {
    	tab[i][0] = svec[i];
    }    
    
    for (int i = 0; i < m+1; i++) {
    	for (int j = 1; j < n+1; j++) {
	    tab[i][j] = smat[i][j-1];
	}   
    }	

// //////////////////////////////////////////////////        
//     for (int i = 0; i < m+1; i++) {
//         for (int j = 0; j < n+1; j++) {
//             printf("%.1f ", tab[i][j]);
//         }
//         printf("\n");        
//     }
//     printf("\n");
// //////////////////////////////////////////////////
    int kIndex;
    bool continueVar = true;
    while (continueVar == true) {
        // Find the entering variable's index  
        double minValue = 1.0;
        int minIndex = 0;
        for (unsigned int i=1; i < n+1; i++) {
            if (tab[0][i] < minValue) {
                minValue = tab[0][i];
                minIndex = i;
            }
        }
        if (minValue>=0) {
            printf("Optimal Value\n");
            continueVar = false;
            break;
        } else {
            kIndex = minIndex;
        }

        // Record k column
        double colK[m];
        double test[m];
        for (unsigned int i=0; i < m+1; i++) {
            colK[i] = tab[i][kIndex];
        }


        // Calculate the theta values
        double theta[m];
        for (unsigned i = 1;  i < n+1; i++) {
            theta[i] =  tab[i][0]/tab[i][kIndex];
        }
        
        // Determine the leaving variables
        int rIndex;
        minValue = 10000.0;
        minIndex = 0;
        for (unsigned int i=1; i < m+1; i++) {
            if (theta[i] < minValue && theta[i] > 0) {
                minValue = theta[i];
                minIndex = i;
            }
        }
        if (minValue==10000.0) {
            printf("unbounded\n");
            continueVar = false;
            break;
        } else {
            rIndex = minIndex;
        }
        
        // Calculate new values for row r
        double divider = tab[rIndex][kIndex];
        for (unsigned int i=0; i < n+1; i++) {
            tab[rIndex][i] = tab[rIndex][i] / divider;
        }
        
        // Calculate new values rest of tableau
        double holderTab[m+1][n+1];
        for (unsigned int i=0; i < m+1; i++) {
            for (unsigned int j=0; j < n+1; j++) {
                if (i != rIndex) {
                    holderTab[i][j] = tab[i][j] - tab[rIndex][j]*tab[i][kIndex];
                } else {
                    holderTab[i][j] = tab[i][j];
                }
            }
        }

        for (unsigned int i=0; i < m+1; i++) {
            for (unsigned int j=0; j < n+1; j++) {
                tab[i][j] = holderTab[i][j];
            }
        }

        double w = colK[rIndex];
        for (unsigned int i=0; i < n+1; i++) {
            if (i!=rIndex) {
                tab[i][kIndex] = -1 * colK[i]/w;
            }
            if (i==rIndex) {
                tab[i][kIndex] = 1/w;
            }
        }
// ////////////////////////////////////////////////        
//       for (int i = 0; i < m+1; i++) {
//            for (int j = 0; j < n+1; j++) {
//                printf("%.1f ", tab[i][j]);
//            }
//            printf("\n");        
//        }
//        printf("\n");
//        printf("iteration end\n\n");
// ////////////////////////////////////////////////
    }
    ////////////////////////////////////////////////        
      for (int i = 0; i < m+1; i++) {
           for (int j = 0; j < n+1; j++) {
               printf("%.1f ", tab[i][j]);
           }
           printf("\n");        
       }
////////////////////////////////////////////////
}