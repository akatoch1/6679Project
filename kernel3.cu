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
if (j == r) return;
tab[j * (n+1) + i] = tab[j * (n+1) + i] - w[threadIdx.y] * tab[r * (n+1) + i];
}