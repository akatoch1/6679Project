__global__ void kernel1(double *tab, double *theta, double *colk, int k, int n) 
{
int i = blockDim.x * blockIdx.x + threadIdx.x;

double w = tab[i * (n+1) + k];

colk[i] = w;
theta[i] = tab[i * (n+1)]/w;

}