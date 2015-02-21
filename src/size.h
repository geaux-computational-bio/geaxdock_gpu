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


#ifndef __SIZE_H_
#define __SIZE_H_

#define     NAMEMAX 64
/* max lig aton name */

#define MAXPRO 10000
/* protein residues */

#define MAXLIG 100
/* ligand heavy atoms */

#define MAXEN1  20
/* protein confs */

#define MAXEN2 70
/* ligand confs */

#define MAXLIG_NUM 200
/* MAX ligand in one ligand .sdf file */

#define MAXTMP 40
/* number of temperature replicas */

#define MAXREP 2048
/* max replicas */

#define  MAXSWP 2000
/* max swapping pairs */

#define  MAXLIB 100
/* library cmps */

#define  MAXSDF  500
/* sdf length */

#define MAXTP1 30
/* point types */

#define MAXTP2 24
/* atom types */

#define MAXTP3 50
/* point types (for ele) */

#define MAXTP4 20
/* residue types */

#define MAXFP1 1024
/* smiles */

#define  MAXFP2 168
/* maccs */

#define  MAXWEI 10
/* energy terms */

#define MAXKDE 10000
/* kde points */

#define MAXPOS 1000
/* position restraints, number of mcs */

#define MAXMCS 256
/* mcs fields, number of field in a mcs */

#define MINLIGRMSD 0.1
/* minimum rmsdf value of ligand ensembles from native */



#define INITTEMP 10000.0f
/* temperature in the first replica */


#define PI 3.1415926535f

/* value used for return statement in data.C for missing atoms */
#define BADKDE 50

// maximum string length for file names
#define MAXSTRINGLENG 128

// if mcs equal to 728492, it is invalid
#define MCS_INVALID_COORD 728492

// if mcs equal to 728492, it is invalid
#define CMCC_INVALID_VAL 10000.0f

// boltzman constant
#define BOLTZMANN_CONST 1.0f

// monte carlo steps
#define STEPS_PER_DUMP 100



#endif

