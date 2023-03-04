__global__ void kernel3(double **tab, double *colk, int k, int r) 
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
tab[j][i] = tab[j][i] - w[threadIdx.y] * tab[r][i];
}