#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(int argc, char *argv[])
{
    int n=argv[0];
    int m=argv[1];
    int b=argv[2];
    int c=argv[3];
    float vArray[m];
    srand(time(NULL));
    int throwAway = round((((float)rand())/RAND_MAX) * b);
    for (int i=0;i<m;i++) {
        vArray[i] = round((((float)rand())/RAND_MAX) * b);
    }

    float aMatrix[n][m];
    for (int i=0;i<n;i++) {
        for (int j=0;j<m;j++) {
            aMatrix[i][j] = round(((((float)rand())/RAND_MAX)-0.5) * b);
        }
    }

    float bArray[n];
    for (int i=0;i<n;i++) {
        bArray[i] = 0;
        for (int j=0;j<m;j++) {
            bArray[i] = bArray[i] +  aMatrix[i][j]*vArray[j];
        }
        bArray[i] = bArray[i] + round(((((float)rand())/RAND_MAX)-0.5) * 3);
    }

    float cArray[m];
    for (int i=0;i<m;i++) {
        cArray[i] = round(((((float)rand())/RAND_MAX)-0.5) * c);
    }


    //Parallelized
    float clockParallelizedStart = clock();
        //Call Function 
    float clockParallelizedEnd = clock();
    float parallelizedTime = (float)(clockParallelizedEnd - clockParallelizedStart)/CLOCKS_PER_SEC;
    printf("%f\n",parallelizedTime);
    //Sequential
    float clockSequentialStart = clock();
        //Call Function 
    float clockSequentialEnd = clock();
    float sequentialTime = (float)(clockSequentialEnd - clockSequentialStart)/CLOCKS_PER_SEC;
    printf("%f\n\n",sequentialTime);
}