#ifndef  GPU_CUH
#define  GPU_CUH



void ExchangeReplicas_d (const int, const int);

void MonteCarlo_Init_d (const int, const int);

void MonteCarlo_d (const int, const int, const int, const int);





void Move_d (Ligand * __restrict__);

void CalcEnergy_d (Ligand * __restrict__, const Protein *);

void Accept_d (Ligand * __restrict__, const float, const int);





// util_d_cpu.C

void InitAcs_d ();

void InitLigRecord_d (const int, const int);

void RecordLigand_d (const int, const int, const int, const int, const Ligand *);

float MyRand_d ();



#endif
