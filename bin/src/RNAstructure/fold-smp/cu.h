#ifndef CU_H
#define CU_H

#ifdef __CUDACC__

#define CU(x) if ((x) != cudaSuccess) {;die("%s:%d: CUDA error: %s", __FILE__, __LINE__,cudaGetErrorString(x)); }
#define DEV __device__
#define HOST __host__
#define GLOBAL __global__

#else

#define DEV
#define HOST
#define GLOBAL

#endif /* __CUDACC__ */

#endif /* CU_H */
