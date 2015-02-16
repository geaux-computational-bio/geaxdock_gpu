
// mode 0: 0 swap 1 , 2 swap 3 , 4 swap 5 , ...
// mode 1: 0 , 1 swap 2 , 3 swap 4 , ...
// mode 3: random
// mode 4: not exchange

__global__ void
ExchangeReplicas_d (const int mode_l, const int mode_t)
{

  if (blockIdx.x == 0) {
    const int bidx = blockDim.x * threadIdx.y + threadIdx.x;      // within a TB

    // should implement a function to check if MAX_XXX is in proper size
    // eg: assert ((MAXREP > n_rep) && (not exaust the shared memory))

    
    //__shared__ int idx_tmp[MAXREP];
    //__shared__ int idx_lig[MAXREP];
    

    __shared__ int temps[MAXTMP];
    __shared__ int temp_orders[MAXTMP];
    __shared__ float energies[MAXTMP];


    if (bidx == 0) {

      for (int p = 0; p < n_prt_dc; ++p) {
	for (int l = 0; l < n_lig_dc; ++l) {

	  // copy index from the replica structure
	  for (int t = 0; t < n_tmp_dc; ++t) {
	     const int flatten_addr =
	       n_tmp_dc * n_lig_dc * p + n_lig_dc * t + l;
	     temps[t] = replica_dc[flatten_addr].idx_tmp;
	     energies[t] = etotal_dc[flatten_addr];
	     temp_orders[temps[t]] = t;
	  }


	  // exchange temperature
	  if (mode_t < 2) {
	    const int maxpair = n_tmp_dc / 2 - (n_tmp_dc % 2 == 0) * (mode_t == 1);
	    for (int pair = 0 ; pair < maxpair; ++pair) {
	      const int i1 = (pair << 1) + mode_t;
	      const int i2 = i1 + 1;
	      const int o1 = temp_orders[i1];
	      const int o2 = temp_orders[i2];
	      const int tt = temps[o1];
	      
	      float etot1 = energies[o1];
	      float etot2 = energies[o2];
	      float minus_beta1 = temp_dc[temps[o1]].minus_beta;
	      float minus_beta2 = temp_dc[temps[o2]].minus_beta;
	      
	      float delta = (minus_beta1 - minus_beta2) * (etot2 - etot1);
	      float exchange_prob = expf(delta);
	      // printf("exchange_prob: %f\n", exchange_prob);
	      // printf("etot1: %f, minus_beta1: %f\n", etot1, minus_beta1);
	      // printf("etot2: %f, minus_beta2: %f\n", etot2, minus_beta2);
	      if (exchange_prob > MyRand_d()) {
		temps[o1] = temps[o2];
		temps[o2] = tt;
		acs_temp_exchg_dc[o1] += 1;
		acs_temp_exchg_dc[o2] += 1;
	      }
	    }
	  }

 
	  // copy index to the replica structure
	  for (int t = 0; t < n_tmp_dc; ++t) {
	    const int flatten_addr =
	      n_tmp_dc * n_lig_dc * p + n_lig_dc * t + l;
	    replica_dc[flatten_addr].idx_tmp = temps[t];
	  }

        }
      }


    }





/*

    for (int r = bidx; r < n_rep_dc; r += TperB) {
      idx_tmp[r] = replica_dc[r].idx_tmp;
      idx_lig[r] = replica_dc[r].idx_lig;
    }
    
    __syncthreads ();

    // exchange temperature
    if (mode_t < 2) {
      int bidx_max = n_rep_dc / 2 + n_rep_dc % 2 - mode_t;
      int left_t = (bidx << 1) + mode_t;
      int right_t = left_t + 1;
      if (bidx < bidx_max) {
	int t = idx_tmp[left_t];
	idx_tmp[left_t] = idx_tmp[right_t];
	idx_tmp[right_t] = t;
      }
    }



    // exchange ligand
    if (mode_l < 2) {
      int bidx_max = n_rep_dc / 2 + n_rep_dc % 2 - mode_t;
      int left_l = bidx << 1 + mode_l;
      int right_l = left_l + 1;
      if (bidx < bidx_max) {
	int l = idx_lig[left_l];
	idx_lig[left_l] = idx_tmp[right_l];
	idx_lig[right_l] = l;
      }
    }
 

    __syncthreads ();

    for (int r = bidx; r < n_rep_dc; r += TperB) {
      replica_dc[r].idx_tmp = idx_tmp[r];
      replica_dc[r].idx_lig = idx_lig[r];
    }

*/

  }


}





