__global__ void kernel3(double *tab, double *colk, int k, int r, int n) 
{

int i = blockDim.x * blockIdx.x + threadIdx.x;
int j = blockIdx.y * blockIdx.y + threadIdx.y;
__shared__ double w[16];
if (threadIdx.y == 0 && threadIdx.x < 16)
{
w[threadIdx.x] = colk[blockIdx.y * blockDim.y + threadIdx.x];
}
__syncthreads();
printf("%f ", tab[i * (n+1) + j]);
if (j == r) return;
//tab[i * (n+1) + j] = tab[i * (n+1) + j] - w[threadIdx.y] * tab[r * (n+1) + i];
tab[i * (n+1) + j] = tab[i * (n+1) + j] - tab[i * (n+1) + k] * tab[r * (n+1) + j];

}