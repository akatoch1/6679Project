__global__ void kernel4(double *tab, double *colk, int k, int r, int n) 
{

int j = blockDim.x * blockIdx.x + threadIdx.x;

__shared__ double w;
if (threadIdx.x == 0) w = colk[r];
__syncthreads();
if (j == r) {
   tab[j*(n+1) + k] = 1/w;
}
else {

     tab[j * (n+1) + k] = -colk[j]/w;

}
}