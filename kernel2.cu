__global__ void kernel2(double *tab, double *colk, int k, int r, int n)
{

int i = blockDim.x * blockIdx.x + threadIdx.x;
if (i >= n+1) {
   return;
}
__shared__ double w;
if (threadIdx.x == 0) w = colk[r];

__syncthreads();


tab[r * (n+1) + i] = tab[r * (n+1) + i]/w;
}