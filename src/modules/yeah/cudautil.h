/*

cuda notes




naming conventions

meaning				MACRO	var		dimensions		ulternative naming
-----------------------------------------------------------------------------------------------------------------
nGPU per node			NGPU	ngpu
nthreads per block		TperB	tperb		TperBx TperBy TperBz	BlockDim, BD, dim_block
nblock per grid			BperG	bperg		BperGx BperGy BperGz	GridDim,  GD, dim_grid
nthreads per grid		TperG	tperg	
thread index of block range		bidx
thread index of grid range		gidx






-------------

CE1 (cudaAPI (...));
CE2 (cudaAPI (...), "Error message");

-------------

kernel <<< ... >>> (...);
cutilCheckMsgAndSync ("kernel failed\n");

-------------

// call this function before the first kernal invocation
cudaGetLastError ();



*/




#ifndef _CUDAUTIL_
#define _CUDAUTIL_

#define CUDAVERSION 5










// force the load to go through the read-only data cache, for CC > 3.5
#if 1
#define CUDA_LDG_P(p) __ldg (p)
#define CUDA_LDG_D(d) __ldg (&d)
#endif
#if 0
#define CUDA_LDG_P(p) *p
#define CUDA_LDG_D(d) d
#endif




// cutil has been removed since CUDA 5
#include <nvidia_gpucomputingsdk_4.2.9_c_common_inc/cutil_inline.h>





// CUDA API Error-Checking Wrapper
inline void CE1(const cudaError_t rv)
{
  if (rv != cudaSuccess) {
    printf("CUDA error %d, %s.\n", rv, cudaGetErrorString (rv));
    exit (1);
  }
}


// CUDA API Error-Checking Wrapper
inline void CE2(const cudaError_t rv, char* pMsg)
{
  if (rv != cudaSuccess) {
    printf("CUDA error %d on %s, %s.\n", rv, pMsg, cudaGetErrorString (rv));
    cutilCheckMsg(pMsg);
    exit (1);
  }
}






#define CUDAMALLOC(ptr, sz, type) \
  cudaMalloc ((void **) &ptr, sz)

#define CUDAFREE(...) \
  cudaFree (__VA_ARGS__)

#define CUDAMEMCPY(dst, src, sz, direction) \
  CE2 (cudaMemcpy (dst, src, sz, direction), #dst)


#if CUDAVERSION == 4
#define CUDAMEMCPYTOSYMBOL(dst, src, type) \
  CE2 (cudaMemcpyToSymbol (#dst, src, sizeof (type), 0, cudaMemcpyHostToDevice), #dst)
#elif CUDAVERSION == 5
#define CUDAMEMCPYTOSYMBOL(dst, src, type) \
  CE2 (cudaMemcpyToSymbol (dst, src, sizeof (type), 0, cudaMemcpyHostToDevice), #dst)
#else
  printf ("cuda version not supportted by this tool set\n");
#endif



// call "cutilDeviceSynchronize" after each kernel
#define CUDAKERNELSYNC(funcname, dim_grid, dim_block, ...) \
  funcname <<< dim_grid, dim_block >>> (__VA_ARGS__); \
  cutilCheckMsgAndSync ("#funcname kernel failed\n")

#define CUDAKERNELASYNC(funcname, dim_grid, dim_block, ...) \
  funcname <<< dim_grid, dim_block >>> (__VA_ARGS__)

#define CUDAKERNELSTREAMSYNC(funcname, dim_grid, dim_block, n, stream, ...) \
  funcname <<< dim_grid, dim_block, n, stream >>> (__VA_ARGS__);	\
  cutilCheckMsgAndSync ("#funcname kernel failed\n")

#define CUDAKERNELSTREAM(funcname, dim_grid, dim_block, n, stream, ...) \
  funcname <<< dim_grid, dim_block, n, stream >>> (__VA_ARGS__);	\
  cutilCheckMsg ("#funcname kernel failed\n")


#endif


