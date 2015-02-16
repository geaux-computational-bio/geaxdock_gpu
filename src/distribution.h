/*
==============================================================================================
     __________________ ____ ___              .___             __      _________   _____   
    /  _____/\______   \    |   \           __| _/____   ____ |  | __ /   _____/  /     \  
   /   \  ___ |     ___/    |   /  ______  / __ |/  _ \_/ ___\|  |/ / \_____  \  /  \ /  \ 
   \    \_\  \|    |   |    |  /  /_____/ / /_/ (  <_> )  \___|    <  /        \/    Y    \
    \______  /|____|   |______/           \____ |\____/ \___  >__|_ \/_______  /\____|__  /
           \/                                  \/           \/     \/        \/         \/ 

      GPU-accelerated hybrid-resolution ligand docking using Replica Exchange Monte Carlo

==============================================================================================
*/

#ifndef __DISTRIBUTION_H_
#define __DISTRIBUTION_H_

////////////////////////////////////////////////////////////////////////////////
// high decoys energy distribuion parameters

#define VDW_NORM_HIGH_LOC          0.05286683f
#define VDW_NORM_HIGH_SCALE        0.29916642f

#define ELE_CAUCHY_HIGH_LOC        -0.3653806f
#define ELE_CAUCHY_HIGH_SCALE      0.00681102f

#define PMF_LOGISTIC_HIGH_LOC      -0.0221777f 
#define PMF_LOGISTIC_HIGH_SCALE    0.11847532f

#define HPC_WALD_HIGH_LOC          -1.0041231f
#define HPC_WALD_HIGH_SCALE        0.1151713f

#define HDB_NORM_HIGH_LOC          0.60986808f
#define HDB_NORM_HIGH_SCALE        0.30974999f

#define DST_LOGISTIC_HIGH_LOC      -0.8246522f
#define DST_LOGISTIC_HIGH_SCALE    0.05218529f

#define PSP_LOGISTIC_HIGH_LOC      -0.2346999f
#define PSP_LOGISTIC_HIGH_SCALE    0.09871713f

#define KDE_WALD_HIGH_LOC          -1.0112055f
#define KDE_WALD_HIGH_SCALE        0.09882933f

#define LHM_LOGISTIC_HIGH_LOC      -0.3377100f
#define LHM_LOGISTIC_HIGH_SCALE    0.09282947f


////////////////////////////////////////////////////////////////////////////////
// low decoys energy distribuion parameters

#define VDW_NORM_LOW_LOC           0.15818486f
#define VDW_NORM_LOW_SCALE         0.29032471f

#define ELE_CAUCHY_LOW_LOC         -0.3647824f
#define ELE_CAUCHY_LOW_SCALE       0.00534225f

#define PMF_LOGISTIC_LOW_LOC       0.04257942f
#define PMF_LOGISTIC_LOW_SCALE     0.08801551f

#define HPC_WALD_LOW_LOC           -1.0019512f
#define HPC_WALD_LOW_SCALE         0.10312611f

#define HDB_LOGISTIC_LOW_LOC       0.76967708f
#define HDB_LOGISTIC_LOW_SCALE     0.11224029f

#define DST_LOGISTIC_LOW_LOC       -0.7509618f
#define DST_LOGISTIC_LOW_SCALE     0.08757299f

#define PSP_LAPLACE_LOW_LOC        -0.2022064f
#define PSP_LAPLACE_LOW_SCALE      0.07854212f

#define KDE_WALD_LOW_LOC           -1.004605f
#define KDE_WALD_LOW_SCALE         0.02967548f

#define LHM_LOGISTIC_LOW_LOC       -0.2387023f
#define LHM_LOGISTIC_LOW_SCALE     0.12134983f

#endif
