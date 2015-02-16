
// _dc stands for gpu constant memory

// array pointers
__constant__ Protein *prt_dc;
__constant__ Psp *psp_dc;
__constant__ Kde *kde_dc;
__constant__ Mcs *mcs_dc;
__constant__ EnePara *enepara_dc;
__constant__ Temp *temp_dc;

__constant__ Ligand *lig_dc;
__constant__ Replica *replica_dc;
__constant__ float *etotal_dc;
__constant__ LigMoveVector *ligmovevector_dc;
__constant__ LigRecord *ligrecord_dc;
__constant__ int *acs_mc_dc;
__constant__ int *acs_temp_exchg_dc;
__constant__ ConfusionMatrix *ref_matrix_dc;


// PRNG seeds
__constant__ int seed_dc;
__constant__ curandState *curandstate_dc;



// monte carlo parameters
__constant__ int steps_per_exchange_dc;
__constant__ int steps_per_dump_dc;
__constant__ int steps_total_dc;
__constant__ float * move_scale_dc;

__constant__ float enepara_lj0_dc;
__constant__ float enepara_lj1_dc;
__constant__ float enepara_el0_dc;
__constant__ float enepara_el1_dc;
__constant__ float enepara_a1_dc;
__constant__ float enepara_b1_dc;
__constant__ float enepara_kde2_dc;
__constant__ float enepara_kde3_dc;


// replica numbers
__constant__ int n_lig_dc;
__constant__ int n_prt_dc;
__constant__ int n_tmp_dc;
__constant__ int n_rep_dc;

// residue numbers (of per replica)
__constant__ int lna_dc;
__constant__ int pnp_dc;
__constant__ int pnk_dc;
__constant__ int pos_dc;



#include "kernel_cuda_l1_exchangereplicas.cu"
#include "kernel_cuda_l1_initcurand.cu"
#include "kernel_cuda_l1_montecarlo.cu"
#include "kernel_cuda_l2_accept.cu"
#include "kernel_cuda_l2_calcenergy.cu"
#include "kernel_cuda_l2_calcmcc.cu"
#include "kernel_cuda_l2_calcrmsd.cu"
#include "kernel_cuda_l2_move.cu"
#include "kernel_cuda_l3_combineenergy.cu"
#include "kernel_cuda_l3_util.cu"

