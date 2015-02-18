/*
 * implement some CUDA APIs in C stdlib
 * to help write universal code for both CPU/GPU
 */


#ifndef _CDUAUTIL_CPU_
#define _CUDAUTIL_CPU_



#define cudaFuncSetCacheConfig(...) \
  ;

#define cudaSetDevice(...) \
  ;

#define cudaGetLastError() \
  ;



#define CUDAMALLOC(ptr, sz, type) \
  ptr = (type) malloc (sz)

#define CUDAFREE(...) \
  free (__VA_ARGS__);

#define CUDAMEMCPY(dst, src, sz, direction) \
  memcpy (dst, src, sz);

#define CUDAMEMCPYTOSYMBOL(dst, src, type)\
  dst = *src;

#define CUDAKERNELSYNC(funcname, dim_grid, dim_block, ...) \
  funcname (__VA_ARGS__);





#endif



