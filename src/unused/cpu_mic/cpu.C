#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "dock.h"
#include "toggle.h"
//#include "cpu.h"
#include "size_gpu.cuh"


// _dc stands for gpu constant memory

// array pointers
Protein *prt_dc;
Psp *psp_dc;
Kde *kde_dc;
Mcs *mcs_dc;
EnePara *enepara_dc;
Temp *temp_dc;

Ligand *lig_dc;
Replica *replica_dc;
float *etotal_dc;
LigMoveVector *ligmovevector_dc;
LigRecord *ligrecord_dc;
TmpEnergy *tmpenergy_dc;
int *acs_dc;




// monte carlo parameters
int steps_per_exchange_dc;
int steps_per_dump_dc;
int steps_total_dc;
int is_random_dc;
float * move_scale_dc;


float enepara_lj0_dc;
float enepara_lj1_dc;
float enepara_el0_dc;
float enepara_el1_dc;
float enepara_a1_dc;
float enepara_b1_dc;
float enepara_kde2_dc;
float enepara_kde3_dc;



// residue numbers (of per replica)
int lna_dc;
int pnp_dc;
int pnk_dc;
int n_pos_dc;

// replica numbers
int n_lig_dc;
int n_prt_dc;
int n_tmp_dc;
int n_rep_dc;




//#include "exchangereplicas_d.cu"
#include "montecarlo_cpu.C"
#include "move_cpu.C"
#include "calcenergy_cpu.C"
#include "accept_cpu.C"
#include "util_d_cpu.C"

