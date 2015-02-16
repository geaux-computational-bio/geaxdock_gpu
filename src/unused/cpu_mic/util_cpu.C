#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <sys/stat.h>

#include "size.h"
#include "toggle.h"
#include "dock.h"
#include "util.h"
#include "load.h"

#include <yeah/quicksort.h>
#include <yeah/mkdir.h>

using namespace std;


void
Usage (char *bin)
{
  fprintf (stderr, "usage: %s [options]\n", bin);
  exit (1);
}


void
ParseArguments (int argc, char **argv, McPara * mcpara, ExchgPara * exchgpara, InputFiles * inputfiles)
{
  // default settings
#if IS_OUTPUT == 1
  char mydir[MAXSTRINGLENG] = "output_default";
#endif
  float floor_temp = 0.3f;   // lowest temperature in all replicas
  float ceiling_temp = 0.3f;   // highest temperature in all replicas
  int num_temp = 4;      // number of temperatures for the same ligand and protein conformations
  float t = 0.001f;  // translational scale
  float r = 3.1415f;  // rotational scale


 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-floor_temp")  && i < argc ) {
    floor_temp = atof(argv[i+1]);
  }
  if ( !strcmp(argv[i],"-ceiling_temp")  && i < argc ) {
    ceiling_temp = atof(argv[i+1]);
  }
  if ( !strcmp(argv[i],"-num_temp")  && i < argc ) {
    num_temp = atoi(argv[i+1]);
  }
  if ( !strcmp(argv[i],"-t")  && i < argc ) {
    t = atof(argv[i+1]);
  }
  if ( !strcmp(argv[i],"-r")  && i < argc ) {
    r = atof(argv[i+1]);
  }
  // if ( !strcmp(argv[i],"-l")  && i < argc ) { compounds_name = string(argv[i+1]); compounds_opt = true; }
  // if ( !strcmp(argv[i],"-s")  && i < argc ) { lhm_name       = string(argv[i+1]); lhm_opt       = true; }
  // if ( !strcmp(argv[i],"-o")  && i < argc ) { output_name    = string(argv[i+1]); output_opt    = true; }
  // if ( !strcmp(argv[i],"-i")  && i < argc ) { molid          = string(argv[i+1]);                       }
  // if ( !strcmp(argv[i],"-nr") && i < argc ) { remc_replicas  = atoi(argv[i+1]);                         }
  // if ( !strcmp(argv[i],"-ns") && i < argc ) { remc_steps     = atoi(argv[i+1]);                         }
  // if ( !strcmp(argv[i],"-nc") && i < argc ) { remc_cycles    = atoi(argv[i+1]);                         }
 }

#if 1
  inputfiles->lig_file.id = "1a07C1";
  inputfiles->lig_file.path = "../test/1a07C1.sdf";
  inputfiles->lig_file.molid = "MOLID";
  inputfiles->prt_file.path = "../test/1a07C.pdb";
  inputfiles->lhm_file.path = "../test/1a07C1.ff";
  inputfiles->lhm_file.ligand_id = "1a07C1";
  inputfiles->enepara_file.path = "../dat/gpudocksm.ff";
  inputfiles->weight_file.path = "../dat/08ff_opt";
  inputfiles->norpara_file.path_a = "../dat/nor_a";
  inputfiles->norpara_file.path_b = "../dat/nor_b";

  exchgpara->floor_temp = floor_temp;
  exchgpara->ceiling_temp = ceiling_temp;
  if (num_temp <= MAXTMP)
    exchgpara->num_temp = num_temp;
  else
    {
      cout << "setting number of temperatures exceeds MAXTMP" << endl;
      cout << "docking exiting ..." << endl;
      exit(1);
    }


  mcpara->steps_total = STEPS_TOTAL;
  mcpara->steps_per_dump = STEPS_PER_DUMP;
  mcpara->steps_per_exchange = STEPS_PER_EXCHANGE;

  // const float t = 1.0f;
  // const float r = 5.0f;

  mcpara->move_scale[0] = t;
  mcpara->move_scale[1] = t;
  mcpara->move_scale[2] = t;
  mcpara->move_scale[3] = r;
  mcpara->move_scale[4] = r;
  mcpara->move_scale[5] = r;

  // random walk toggle
  mcpara->is_random = IS_RANDOM;


#if IS_OUTPUT == 1
  // create output directory

  // obtain the time tag
  char mystime[MAXSTRINGLENG];
  time_t mytime = time (NULL);
  struct tm *mylocaltime;
  mylocaltime = localtime (&mytime);
  strftime (mystime, MAXSTRINGLENG, "%Y%m%d_%H%M%S", mylocaltime);

  // name the output directory using time tag
  if (strcmp (mydir, "output_default") == 0) {
    sprintf (mydir, "output_%s", mystime);
  }

  strcpy (mcpara->outputdir, mydir);
  MakeDir (mydir);

  // prefix of the file name
  const char h5file[MAXSTRINGLENG] = "a";
  strcpy (mcpara->outputfile, h5file);
#endif

#endif

}





void
OptimizeLigand (const Ligand0 * lig0, Ligand * lig, const ComplexSize complexsize)
{

  // data structure translation
  for (int i = 0; i < complexsize.n_lig; ++i) {
    const Ligand0 *src = &lig0[i];
    Ligand *dst = &lig[i];

    for (int residue = 0; residue < MAXLIG; ++residue) {
      dst->lig_point[residue].t = src->t[residue];
      dst->lig_point[residue].c = src->c[residue];
      dst->lig_point[residue].n = src->n[residue];
    }
    dst->lna = src->lna;



    // generate coord_orig
    dst->coord_orig = src->coord_orig;
    dst->coord_orig.center[0] = src->pocket_center[0];
    dst->coord_orig.center[1] = src->pocket_center[1];
    dst->coord_orig.center[2] = src->pocket_center[2];

  }

}

// move the ligand to its center
void InitLigCoord (Ligand * lig, const ComplexSize complexsize){
  for (int i = 0; i < complexsize.n_lig; ++i) {
    Ligand *mylig = &lig[i];
    
    ////////////////////////////////////////////////////////////////////////////////
    // debugging, move ligand far away from center
    mylig->coord_orig.center[0] += 1.5f;
    mylig->coord_orig.center[1] += 1.5f;
    mylig->coord_orig.center[2] += 1.5f;

    mylig->coord_new = mylig->coord_orig;
    for (int residue = 0; residue < mylig->lna; residue++) {
      mylig->coord_new.x[residue] += mylig->coord_new.center[0];
      mylig->coord_new.y[residue] += mylig->coord_new.center[1];
      mylig->coord_new.z[residue] += mylig->coord_new.center[2];
    }

    for (int i = 0; i < 6; ++i) {
      mylig->movematrix_old[i] = 0.0f;
    }
  }
}




void
CopyProteinResidue (const Protein0 * src, Protein * dst, const int residue_src,
		    const int residue_dst, const EnePara0 * enepara0)
{
  dst->prt_point[residue_dst].x = src->x[residue_src];
  dst->prt_point[residue_dst].y = src->y[residue_src];
  dst->prt_point[residue_dst].z = src->z[residue_src];

  dst->prt_point[residue_dst].t = src->t[residue_src];
  dst->prt_point[residue_dst].c = src->c[residue_src];

  int t = src->t[residue_src];
  int d = src->d[residue_src];
  int dt = t == 0 ? src->d[residue_src] + 30 : t;
  dst->prt_point[residue_dst].ele = enepara0->ele[dt];

  dst->prt_point[residue_dst].seq3r = src->seq3[src->r[residue_src]];

  dst->prt_point[residue_dst].c0_and_d12_or_c2 =
    (src->c[residue_src] == 0 && src->d[residue_src] == 12) || (src->c[residue_src] == 2);

  dst->prt_point[residue_dst].hpp = enepara0->hpp[d];

}

void
OptimizeProtein (const Protein0 * prt0, Protein * prt, const EnePara0 * enepara0,
		 const Ligand0 * lig0, const ComplexSize complexsize)
{
  // pocket center
  const float cx = lig0[0].pocket_center[0];
  const float cy = lig0[0].pocket_center[1];
  const float cz = lig0[0].pocket_center[2];

  for (int i = 0; i < complexsize.n_prt; ++i) {
    const Protein0 *src = &prt0[i];
    Protein *dst = &prt[i];
    const int pnp = src->pnp;

    // sort protein in increament order of t
    int *t = (int *) malloc (sizeof (int) * pnp);
    int *order = (int *) malloc (sizeof (int) * pnp);
    for (int j = 0; j < pnp; ++j) {
      t[j] = src->t[j];
      order[j] = j;
    }

#if 0
    QuickSort (t, order, 0, pnp - 1);
#endif

#if 0
    for (int j = 0; j < pnp; ++j) {
      printf ("%d ", t[j]);
    }
    putchar ('\n');

    for (int j = 0; j < pnp; ++j) {
      printf ("%d ", src->t[order[j]]);
    }
    putchar ('\n');
#endif

    dst->pnp = pnp;
    for (int j = 0; j < pnp; ++j)
      CopyProteinResidue (src, dst, order[j], j, enepara0);

    free (t);
    free (order);

    // assign the pocket center from the ligand structure to the protein sturcture
    dst->pocket_center[0] = cx;
    dst->pocket_center[1] = cy;
    dst->pocket_center[2] = cz;

  }

}





void
OptimizePsp (const Psp0 * psp0, Psp * psp, const Ligand * lig, const Protein * prt)
{
  for (int i = 0; i < MAXPRO; ++i) {
    for (int j = 0; j < MAXLIG; ++j) {
      psp->psp[j][i] = psp0->psp[i][j];
    }
  }

/*
  int lig_lna = lig[0].lna;
  int prt_pnp = prt[0].pnp;
  int count = 0;
*/

/*
  for (int i = 0; i < lig_lna; ++i) {
    int lig_t = lig[0].t[i];
    for (int j = 0; j < prt_pnp; ++j) {
      int pspidx2 = prt[0].seq3r[j];

      if (prt[0].c[j] == 2) {
	printf ("potentialy accessed psp \t%2d\t%2d\n", lig_t, pspidx2);
	count++;

      if (psp->psp[i][j] != 0)
        printf ("lig %3d \tprt %3d \t\tpsp\t%f\n", i, j, psp->psp[j][i]);

      }
    }
  }
*/

/*
  printf ("percentage %f / %f = %f", count, lig_lna * prt_pnp,
	  (float) count / (lig_lna * prt_pnp));
*/

}

void
OptimizeKde (const Kde0 * kde0, Kde * kde)
{
  for (int i = 0; i < MAXKDE; ++i) {
    kde->kde_point[i].x = kde0->x[i];
    kde->kde_point[i].y = kde0->y[i];
    kde->kde_point[i].z = kde0->z[i];
    kde->kde_point[i].t = kde0->t[i];
  }
  kde->pnk = kde0->pnk;
}

void
OptimizeMcs (const Mcs0 * mcs0, Mcs * mcs, const ComplexSize complexsize)
{
  // n_pos
  for (int i = 0; i < complexsize.n_pos; ++i) {
    mcs[i].tcc = mcs0[i].tcc;

    // n_mcs
    for (int j = 0; j < MAXMCS; ++j) {
      mcs[i].mcs_point[j].x = mcs0[i].x[j];
      mcs[i].mcs_point[j].y = mcs0[i].y[j];
      mcs[i].mcs_point[j].z = mcs0[i].z[j];
    }
  }

}

void
OptimizeEnepara (const EnePara0 * enepara0, EnePara * enepara)
{
  const float sqrt_2_pi = sqrtf (2.0f * PI);

  for (int i = 0; i < MAXTP2; ++i) {	// lig
    for (int j = 0; j < MAXTP1; ++j) {	// prt
      const float tmp = enepara0->vdw[j][i][0] * enepara0->lj[2];
      enepara->pa[i][j][0] = 2.0f * enepara0->vdw[j][i][1] * powf (tmp, 9.0f);
      enepara->pa[i][j][1] = 3.0f * enepara0->vdw[j][i][1] * powf (tmp, 6.0f);
    }
  }
  enepara->lj0 = enepara0->lj[0];
  enepara->lj1 = enepara0->lj[1];

  enepara->el0 = enepara0->el[0];
  enepara->a1 = 4.0f - 3.0f * enepara0->el[0];
  enepara->b1 = 2.0f * enepara0->el[0] - 3.0f;
  enepara->el1 = enepara0->el[1];

  for (int i = 0; i < MAXTP3; ++i)
    enepara->ele[i] = enepara0->ele[i];

  for (int i = 0; i < MAXTP2; ++i) {	// lig
    for (int j = 0; j < MAXTP1; ++j) {	// prt
      enepara->pmf[i][j][0] = enepara0->pmf[j][i][0];
      enepara->pmf[i][j][1] = enepara0->pmf[j][i][1];
      enepara->hdb[i][j][0] = enepara0->hdb[j][i][0];
      enepara->hdb[i][j][1] = 1.0f / enepara0->hdb[j][i][1];
      enepara->hdb[i][j][2] = -1.0f / (enepara0->hdb[j][i][1] * sqrt_2_pi);
    }
  }

  for (int i = 0; i < MAXTP4; ++i)
    enepara->hpp[i] = enepara0->hpp[i];
  for (int i = 0; i < MAXTP2; ++i) {
    enepara->hpl[i][0] = enepara0->hpl[i][0];
    enepara->hpl[i][1] = enepara0->hpl[i][1];
    enepara->hpl[i][2] = logf (1.0f / (enepara->hpl[i][1] * sqrt_2_pi));
  }

  enepara->kde2 = -0.5f / (enepara0->kde * enepara0->kde);
  enepara->kde3 = powf (enepara0->kde * sqrt_2_pi, 3.0f);

  for (int i = 0; i < MAXWEI; ++i) {
    enepara->w[i] = enepara0->w[i];
    // cout << enepara->w[i] << endl;
  }

  for (int i = 0; i < MAXWEI; ++i) {
    enepara->a_para[i] = enepara0->a_para[i];
    enepara->b_para[i] = enepara0->b_para[i];
    // cout << enepara->w[i] << endl;
  }
}




/*
void
SetWeight (EnePara * enepara)
{
#if 0
  for (int i = 0; i < MAXWEI; i++)
    enepara->w[i] = 1.0f;

  enepara->w[i] = 1.0f;
  enepara->w[i] = 1.0f;
  enepara->w[i] = 1.0f;
  enepara->w[i] = 1.0f;
  enepara->w[i] = 1.0f;
  enepara->w[i] = 1.0f;
  enepara->w[i] = 1.0f;
  enepara->w[i] = 1.0f;
  enepara->w[i] = 1.0f;
#endif

#if 1
  for (int i = 0; i < MAXWEI; i++)
    enepara->w[i] = 0.0f;
  enepara->w[0] = 1.0f;
#endif
}
*/

void
SetTemperature (Temp * temp, ExchgPara * exchgpara)
{
  int num_temp = exchgpara->num_temp;
  float floor_temp = exchgpara->floor_temp;
  float ceiling_temp = exchgpara->ceiling_temp;
  
  float beta_high = 1.0f / floor_temp;
  float beta_low = 1.0f / ceiling_temp;
  const float beta_ratio = exp (log (beta_high / beta_low) / (float) (num_temp - 1));
  // cout << "beta_ratio: " << beta_ratio << endl;

  float a = beta_low;
  for (int i = 0; i < num_temp; i++) {
    temp[i].order = i;
    temp[i].minus_beta = 0.0f - a;

    a *= beta_ratio;
  }

  // for (int i = 0; i < num_temp; i++) {
  //   temp[i].t = floor_temp;
  //   temp[i].minus_beta = -1.0f / temp[i].t;
  //   temp[i].order = i;
  // }
}

// replica[n_rep]
// replica[n_prt][n_tmp][n_lig]

void
SetReplica (Replica * replica, Ligand * lig, const ComplexSize complexsize)
{
  const int n_lig = complexsize.n_lig;
  const int n_prt = complexsize.n_prt;
  const int n_tmp = complexsize.n_tmp;

  for (int i = 0; i < n_prt; ++i) {
    for (int j = 0; j < n_tmp; ++j) {
      for (int k = 0; k < n_lig; ++k) {
	const int flatten_addr = n_tmp * n_lig * i + n_lig * j + k;
	replica[flatten_addr].idx_rep = flatten_addr;
	replica[flatten_addr].idx_prt = i;
	replica[flatten_addr].idx_tmp = j;
	replica[flatten_addr].idx_lig = k;

	lig[flatten_addr] = lig[k];	// duplicate ligand replicas
      }
    }
  }

}


void
SetMcLog (McLog * mclog)
{
  mclog->t0 = 0;
  mclog->t1 = 0;
  mclog->t2 = 0;
  mclog->ac = 0;

  for (int i = 0; i < MAXREP; ++i)
    mclog->acs[i] = 0;
}



// arg = 1      print title
// arg = 2      print content
// arg = 3      print all

void
PrintEnergy1 (const Energy * energy, const int step, const int arg)
{
  int a = arg & 0x1;
  int b = (arg >> 1) & 0x1;

  if (a == 1) {
    std::cout << "step" << ","
      << "etot" << "," << "elhm" << ","
      << "ekde" << "," << "edst" << ","
      << "eele" << "," << "epmf" << "," << "ehpc" << "," << "ehdb" << "," << "epsp" << "," << "evdw"
      << std::endl;
  }

  if (b == 1) {
    std::cout << step << ","
	      << energy->e[9] << ","    // total
	      << energy->e[7] << ","	// lhm
	      << energy->e[6] << ","	// kde
	      << energy->e[8] << ","	// dst
	      << energy->e[1] << ","	// ele
	      << energy->e[2] << ","	// pmf
	      << energy->e[5] << ","	// hpc
	      << energy->e[4] << ","	// hdb
	      << energy->e[3] << ","	// psp
	      << energy->e[0] << std::endl;	// vdw
  }

}

// arg = 1      print title
// arg = 2      print content
// arg = 3      print all


void
PrintCsv (const Energy * energy, const int idx_rep, const int step,
	      const int arg)
{

  char names[MAXWEI][8] = {
    "vdw",
    "ele",
    "pmf",
    "psp",
    "hdb",
    "hpc",
    "kde",
    "lhm",
    "dst",
    "total"
  };

  int a = arg & 0x1;
  int b = (arg >> 1) & 0x1;

  if (a == 1) {
    printf ("rep,step");
    for (int i = 0; i < MAXWEI; ++i)
      printf (",%s", names[i]);
    printf ("\n");
  }

  if (b == 1) {
    printf ("%4d,%5d", idx_rep, step);
    for (int i = 0; i < MAXWEI; ++i)
      printf (",%+14.10f", energy->e[i]);
    printf ("\n");
  }
}

void
PrintEnergy2 (const Energy * energy, const int idx_rep, const int step,
	      const int arg)
{

  char names[MAXWEI][8] = {
    "vdw",
    "ele",
    "pmf",
    "psp",
    "hdb",
    "hpc",
    "kde",
    "lhm",
    "dst",
    "total"
  };

  int a = arg & 0x1;
  int b = (arg >> 1) & 0x1;

  if (a == 1) {
    printf ("rep \tstep \t");
    for (int i = 0; i < MAXWEI; ++i)
      printf ("\t\t%s", names[i]);
    printf ("\n");
  }

  if (b == 1) {
    printf ("%4d \t%5d \t\t", idx_rep, step);
    for (int i = 0; i < MAXWEI; ++i)
      //printf ("\t%+.3e", energy->e[i]);
      printf ("\t%+14.10f", energy->e[i]);
    printf ("\n");
  }
}



void
PrintEnergy3 (const Energy * energy, const int idx_rep, const int step,
	      const int track, const int arg)
{

  char names[MAXWEI][8] = {
    "vdw",
    "ele",
    "pmf",
    "psp",
    "hdb",
    "hpc",
    "kde",
    "lhm",
    "dst",
    "total"
  };

  int a = arg & 0x1;
  int b = (arg >> 1) & 0x1;

  if (a == 1) {
    printf ("rep \tstep \ttrack ");
    for (int i = 0; i < MAXWEI; ++i)
      printf ("\t\t%s", names[i]);
    printf ("\n");
  }

  if (b == 1) {
    printf ("%4d \t%5d \t%d \t", idx_rep, step, track);
    for (int i = 0; i < MAXWEI; ++i)
      //printf ("\t%+.3e", energy->e[i]);
      printf ("\t%+10.6f", energy->e[i]);
    printf ("\n");
  }
}



void
PrintMoveVector (const float m[6], const int step)
{
  printf ("\t  %3d\t\t\t", step);
  for (int i = 0; i < 6; ++i) {
    printf (" %+f\t", m[i]);
  }
  printf("\n");
}



void
PrintMoveRecord (const LigRecord * ligrecord, const int steps_per_dump, const int replica,
		const int iter_begin, const int iter_end, const int arg)
{
  for (int s = iter_begin; s <= iter_end; ++s) {
    const LigRecordSingleStep *myrecord = &ligrecord[replica].step[s];
    PrintMoveVector (myrecord->movematrix, myrecord->step);
  }
  
}



// arg = 1      print title
// arg = 2      print content
// arg = 3      print all

void
PrintLigRecord (LigRecord * ligrecord,  int steps_per_dump,  int replica,
		int iter_begin,  int iter_end,  int arg)
{
  // print title
  PrintEnergy2 (NULL, NULL, NULL, 1);
  //PrintCsv (NULL, NULL, NULL, 1);

  for (int s = iter_begin; s <= iter_end; ++s) {
    const LigRecordSingleStep *myrecord = &ligrecord[replica].step[s];
    PrintEnergy2 (&myrecord->energy, replica, myrecord->step, arg);
    //PrintCsv (&myrecord->energy, replica, myrecord->step, arg);
  }

}





void
PrintRepRecord (const LigRecord * ligrecord, const int steps_per_dump, const int rep_begin,
		const int rep_end, const int iter_begin, const int iter_end, const int arg)
{
  printf ("\treplicas\n");

  printf ("step|\t");

  for (int r = rep_begin; r <= rep_end; ++r)
    printf ("%2d\t", r);
  putchar ('\n');

  printf ("----+");

  for (int r = rep_begin; r <= rep_end; ++r)
    printf ("--------");
  putchar ('\n');

  for (int s = iter_begin; s <= iter_end; ++s) {
    printf ("%3d |\t", s);

    for (int r = rep_begin; r <= rep_end; ++r) {
      const Replica *myrep = &ligrecord[r].step[s].replica;
      //printf ("%2d ", myrep->idx_prt);
      //printf ("%2d ", myrep->idx_tmp);
      //printf ("%2d ", myrep->idx_lig);

      printf ("%2d ", myrep->idx_rep);

      printf ("\t");
    }
    putchar ('\n');
  }

}





// print all temperature replicas of the same lig & prt
void
PrintRepRecord2 (LigRecord * ligrecord,  ComplexSize complexsize,
		 int steps_per_dump,  int idx_prt,  int idx_lig,
		 int iter_begin,  int iter_end,  int arg)
{
  printf ("temperature replicas with lig %d prt %d\n", idx_lig, idx_prt);

  printf ("MC step |\t");

  for (int t = 0; t < complexsize.n_tmp; ++t)
    printf ("%2d\t", t);
  putchar ('\n');

  printf ("--------+----");

  for (int t = 0; t < complexsize.n_tmp; ++t)
    printf ("--------");
  putchar ('\n');

  for (int s = iter_begin; s <= iter_end; ++s) {
    printf ("%5d   |\t", s);

    for (int t = 0; t < complexsize.n_tmp; ++t) {
      const int r =
	complexsize.n_tmp * complexsize.n_lig * idx_prt +
	complexsize.n_lig * t +
	idx_lig;
      const Replica *myrep = &ligrecord[r].step[s].replica;
      //printf ("%2d ", myrep->idx_prt);
      printf ("%2d ", myrep->idx_tmp);
      //printf ("%2d ", myrep->idx_lig);
      //printf ("%2d ", myrep->idx_rep);

      printf ("\t");
    }
    putchar ('\n');
  }
}






void
PrintLigCoord (const Ligand * lig, const int unused)
{
  const Ligand *mylig = &lig[0];
  const LigCoord *mycoord = &mylig->coord_new;

  //int range = mylig->lna;
  int range = 6;

  printf ("---------------------------------------------------------------------\n");

  for (int i = 0; i < range; ++i)
    printf ("%9d  ", i);
  putchar ('\n');
  for (int i = 0; i < range; ++i)
    printf ("%+9.4f  ", mycoord->x[i]);
  putchar ('\n');
  for (int i = 0; i < range; ++i)
    printf ("%+9.4f  ", mycoord->y[i]);
  putchar ('\n');
  for (int i = 0; i < range; ++i)
    printf ("%+9.4f  ", mycoord->z[i]);
  putchar ('\n');

  printf ("center xyz\t\t: %f  %4f  %4f\n", mycoord->center[0], mycoord->center[1],
	  mycoord->center[2]);

  printf ("---------------------------------------------------------------------\n");
}

void
PrintLigand (const Ligand * lig)
{

  // const LigCoord *mycoord = &lig->coord_new;
  const LigCoord *mycoord = &lig->coord_orig;
  printf ("center:\t\t%+10.6f\t%+10.6f\t%+10.6f\n", mycoord->center[0], mycoord->center[1], mycoord->center[2]);
  printf ("lna:\t\t%d\n", lig->lna);

  printf ("x \t\ty \t\tz \t\tc \t\t t \t n \tindex\n");
  printf ("-----------------------------------------------\n");
  const int lna = lig->lna;
  for (int i = 0; i < lna; ++i) {
    printf ("%+10.6f\t", mycoord->x[i]);
    printf ("%+10.6f\t", mycoord->y[i]);
    printf ("%+10.6f\t", mycoord->z[i]);
    printf ("%+10.6f\t", lig->lig_point[i].c);
    printf ("%2d\t", lig->lig_point[i].t);
    printf ("%2d\t", lig->lig_point[i].n);
    printf ("%3d\n", i);
  }

}

void
PrintProtein (const Protein * prt)
{
  printf ("pnp:\t\t%d\n", prt->pnp);

  printf ("x \t\ty \t\tz \t\t t \t c \t d \tindex\n");
  printf ("-----------------------------------------------\n");
  const int pnp = prt->pnp;
  for (int i = 0; i < pnp; ++i) {
    printf ("%+10.6f\t", prt->prt_point[i].x);
    printf ("%+10.6f\t", prt->prt_point[i].y);
    printf ("%+10.6f\t", prt->prt_point[i].z);
    printf ("%2d\t", prt->prt_point[i].t);
    printf ("%2d\t", prt->prt_point[i].c);
    printf ("%4d\n", i);
  }

}

void
PrintDataSize (const Ligand * lig,
	       const Protein * prt, const Psp * psp, const Kde * kde, const Mcs * mcs,
	       const EnePara * enepara)
{
  float lig_sz, prt_sz, psp_sz, kde_sz, mcs_sz, enepara_sz;

  lig_sz = sizeof (Ligand);
  prt_sz = sizeof (Protein);
  psp_sz = sizeof (Psp);
  kde_sz = sizeof (Kde);
  mcs_sz = sizeof (Mcs);
  enepara_sz = sizeof (EnePara);

  printf ("lig \t\tprt \t\tpsp \t\tkde \t\tmcs \t\tenepara\n");
  printf ("%f \t%f \t%f \t%f \t%f \t%f\t\t",
	  lig_sz / 1024, prt_sz / 1024, psp_sz / 1024, kde_sz / 1024, mcs_sz / 1024,
	  enepara_sz / 1024);
  printf ("full size (KB)\n");

  lig_sz = (3 * 2 + 4) * lig->lna * 4;
  prt_sz = (8) * prt->pnp * 4;
  psp_sz = 999;
  kde_sz = (4) * kde->pnk * 4;
  mcs_sz = 999;
  enepara_sz = 999;

  printf ("%f \t%f \t%f \t%f \t%f \t%f\t\t",
	  lig_sz / 1024, prt_sz / 1024, psp_sz / 1024, kde_sz / 1024, mcs_sz / 1024,
	  enepara_sz / 1024);
  printf ("effective size (KB)\n");
}

void
PrintSummary (const InputFiles * inputfiles, const McPara * mcpara, const Temp * temp,
	      const McLog * mclog, const ComplexSize * complexsize)
{
  putchar ('\n');

  printf ("ligand file\t\t\t");
  std::cout << inputfiles->lig_file.path << std::endl;
  printf ("protein file\t\t\t");
  std::cout << inputfiles->prt_file.path << std::endl;
  printf ("lhm file\t\t\t");
  std::cout << inputfiles->lhm_file.path << std::endl;
  printf ("enepara file\t\t\t");
  std::cout << inputfiles->enepara_file.path << std::endl;
  printf ("weight file\t\t\t");
  std::cout << inputfiles->weight_file.path << std::endl;

  printf ("output directory\t\t%s\n", mcpara->outputdir);
  printf ("out file (HDF5)\t\t\t%s/%s_XXXX.h5\n", mcpara->outputdir, mcpara->outputfile);

  printf ("steps_total\t\t\t%d\n", mcpara->steps_total);
  printf ("steps_per_dump\t\t\t%d\n", mcpara->steps_per_dump);
  printf ("steps_per_exchange\t\t%d\n", mcpara->steps_per_exchange);

  const size_t ligrecord_sz = sizeof (LigRecord) * complexsize->n_rep;
  printf ("per dump record size:\t\t%.3f MB\n", (float) ligrecord_sz / 1048576);

  printf ("movement scale (txyz, rxyz)\t");
  for (int i = 0; i < 6; ++i)
    printf("%.8f ", mcpara->move_scale[i]);
  printf ("\n");

  printf ("total replica\t\t\t%d\n", complexsize->n_rep);
  printf ("total prt conformations\t\t%d\n", complexsize->n_prt);
  printf ("total temperatures\t\t%d\n", complexsize->n_tmp);
  printf ("total ligand conformations\t%d\n", complexsize->n_lig);
  printf ("mcs\t\t\t\t%d\n", complexsize->n_pos);

  printf ("AR of total \t\t\t%d / %d \t%f\n",
	  mclog->ac,
	  mcpara->steps_total * complexsize->n_rep,
	  (float) mclog->ac / (mcpara->steps_total * complexsize->n_rep));

#if 0
  for (int t = 0; t < complexsize->n_tmp; ++t) {
    const int myreplica = complexsize->n_lig * t;
    printf ("AR of temperature[%d]=%f \t %d / %d \t%f\n",
	    t,
	    temp[t].t,
	    mclog->acs[myreplica],
	    mcpara->steps_total, (float) mclog->acs[myreplica] / mcpara->steps_total);
  }
#endif

#if 0

  for (int i = 0; i < complexsize->n_prt; ++i) {
    for (int j = 0; j < complexsize->n_tmp; ++j) {
      for (int k = 0; k < complexsize->n_lig; ++k) {
	const int flatten_addr =
	  complexsize->n_tmp * complexsize->n_lig * i + complexsize->n_lig * j + k;
	printf ("AR of %4d temperature[%d]=%f \t %d / %d \t%f\n", flatten_addr, j, temp[j].t,
		mclog->acs[flatten_addr], mcpara->steps_total,
		(float) mclog->acs[flatten_addr] / mcpara->steps_total);
      }
    }
  }
#endif


  const float mcpersec0 = mcpara->steps_total * complexsize->n_rep / mclog->t0;
  printf ("compute time\t\t\t%.3f seconds\n", mclog->t0);
  printf ("time per MC sweep per replica\t%.3f * 1e-6 seconds\n", 1e6 / mcpersec0);
  printf ("MC sweeps per second\t\t%.3f\n", mcpersec0);
  printf ("speedup over 843.75\t\t%.3f X\n", mcpersec0 / 843.75);

  const float mcpersec1 = mcpara->steps_total * complexsize->n_rep / mclog->t1;
  printf ("wall time\t\t\t%.3f seconds\n", mclog->t1);
  printf ("time per MC sweep per replica\t%.3f * 1e-6 seconds \n", 1e6 / mcpersec1);
  printf ("MC sweeps per second\t\t%.3f\n", mcpersec1);
  printf ("speedup over 843.75\t\t%.3f X\n", mcpersec1 / 843.75);


}

/*
float MyRand()
{
	//float randdd = (float) rand () / RAND_MAX;
	float randdd = 0.002f;
	return randdd;
}
*/


int
minimal_int (const int a, const int b)
{
  return a < b ? a : b;
}

