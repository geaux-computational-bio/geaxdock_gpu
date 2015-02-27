#ifndef  GPU_CUH
#define  GPU_CUH


__global__ void InitCurand_d ();

__global__ void ResetCounter_d (const int, const int);

__global__ void ExchangeReplicas_d (const int, const int);

__global__ void MonteCarlo_Init_d (const int, const int);

__global__ void MonteCarlo_d (const int, const int, const int, const int);





__device__ void CalcRmsd_d (const int, Ligand * __restrict__);

__device__ void Move_d (const int, Ligand * __restrict__, const float);

__device__ void CalcEnergy_d (const int, Ligand * __restrict__, const Protein * __restrict__);

__device__ void CombineEnergy_d (const int, Energy *);

__device__ void CalcMcc_d (const int, Ligand * __restrict__, const Protein * __restrict__);

__device__ void InitRefMatrix_d (const int, Ligand * __restrict__, const Protein * __restrict__);

__forceinline__ __device__ void Accept_d (const int, Ligand * __restrict__, const float, const int);




// util_d.cu

__device__ void InitAcs_d (const int);

//__device__ void InitLigRecord_d (const int, const int, const int);

//__forceinline__ __device__ void BackupLigCoord_d (const int, Ligand *);

__device__ void RecordLigand_d (const int, const int, const int, const int, const int, const Ligand *);



__forceinline__ __device__ float MyRand_d ();

//__forceinline__ __device__ int Minimal_int_d (const int, const int);

__forceinline__ __device__ void SumReduction1D_d (const int, float *);

__forceinline__ __device__ void SumReduction1D_5_d (const int, float *, float *, float *, float *, float *);

__forceinline__ __device__ void SumReduction_int_1D_4_d (const int, int *, int *, int *, int *);

__forceinline__ __device__ void SumReduction2D_d (float a[BDy][BDx]);

__forceinline__ __device__ void SumReduction2D_2_d (float a[BDy][BDx], float b[BDy][BDx]);



__forceinline__ __device__ float NormPdf(float x, float loc, float scale);

__forceinline__ __device__ float CauchyPdf(float x, float loc, float scale);

__forceinline__ __device__ float LogisticPdf(float x, float loc, float scale);

__forceinline__ __device__ float WaldPdf(float x, float loc, float scale);

__forceinline__ __device__ float LaplacePdf(float x, float loc, float scale);


#endif
