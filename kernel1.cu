__global__ void vecAddKernel(double **tab, double *theta, double *colk, int k) 
{
int i = blockDim.x * blockIdx.x + threadIdx.x;
double w = tab[i][k];
colk[i] = w;
theta[i] = tab[i][1]/w;

}