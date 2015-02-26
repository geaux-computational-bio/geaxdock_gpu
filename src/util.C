#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <list>
#include <vector>
#include <string>
#include <stdio.h>
#include <sstream>
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
Banner ()
{
  cout << "------------------------------------------------------------" << endl
       << "                         GeauxDock                          " << endl
       << "                        version 0.1                         " << endl << endl
       << "   GPU-accelerated mixed-resolution ligand docking using    " << endl
       << "                Replica Exchange Monte Carlo                " << endl
       << "------------------------------------------------------------" << endl << endl;
}


void
TraceBanner ()
{
  cout << "------------------------------------------------------------" << endl
       << "                         GeauxDock                          " << endl
       << "                        version 0.1                         " << endl << endl
       << "         Generate ligand conformation in sdf format     " << endl
       << "                    given the move vector" << endl
       << "------------------------------------------------------------" << endl << endl;
}




void
ParseArguments (int argc, char **argv, McPara * mcpara, ExchgPara * exchgpara,
		InputFiles * inputfiles)
{

  ////////////////////////////////////////////////////////////////////////////////
  // default settings
  float ts = 0.001f;		// translational scale
  float rs = 3.1415f;		// rotational scale

  exchgpara->floor_temp = 0.3f;
  exchgpara->ceiling_temp = 0.3f;

  mcpara->steps_total = 3000;
  mcpara->steps_per_dump = STEPS_PER_DUMP;
  mcpara->steps_per_exchange = 5;


#if IS_BAYE == 1
  inputfiles->norpara_file.path_a = "baye_nor_a";
  inputfiles->norpara_file.path_b = "baye_nor_b";
#elif IS_BAYE == 0
  inputfiles->norpara_file.path_a = "08_nor_a";
  inputfiles->norpara_file.path_b = "08_nor_b";
  // inputfiles->norpara_file.path_a = "08_nor_a";
  // inputfiles->norpara_file.path_b = "08_nor_b";
#endif

  inputfiles->enepara_file.path = "gpudocksm.ff";
  inputfiles->lig_file.molid = "MOLID";


  bool protein_on = false;
  bool compounds_on = false;
  bool lhm_one = false;
  bool hdf_on = false;

  for (int i = 0; i < argc; i++) {

    // complex files

    if (!strcmp (argv[i], "-id") && i < argc) {
      inputfiles->lhm_file.ligand_id = argv[i + 1];
    }
	// *.pdb
    if (!strcmp (argv[i], "-p") && i < argc) {
      inputfiles->prt_file.path = argv[i + 1];
      protein_on = true;
    }
	// *.sdf
    if (!strcmp (argv[i], "-l") && i < argc) {
      inputfiles->lig_file.path = argv[i + 1];
      compounds_on = true;
    }
	// *.ff
    if (!strcmp (argv[i], "-s") && i < argc) {
      inputfiles->lhm_file.path = argv[i + 1];
      lhm_one = true;
    }



	// parameter files

    if (!strcmp (argv[i], "-opt") && i < argc) {
      inputfiles->weight_file.path = argv[i + 1];
    }
    if (!strcmp (argv[i], "-na") && i < argc) {
      inputfiles->norpara_file.path_a = argv[i + 1];
    }
    if (!strcmp (argv[i], "-nb") && i < argc) {
      inputfiles->norpara_file.path_b = argv[i + 1];
    }
    if (!strcmp (argv[i], "-para") && i < argc) {
      inputfiles->enepara_file.path = argv[i + 1];
    }



    // MC steps

    if (!strcmp (argv[i], "-ns") && i < argc) {
      mcpara->steps_total = atoi (argv[i + 1]);
    }
    if (!strcmp (argv[i], "-nc") && i < argc) {
      mcpara->steps_per_exchange = atoi (argv[i + 1]);
    }




    // temperatures

    if (!strcmp (argv[i], "-floor_temp") && i < argc) {
      exchgpara->floor_temp = atof (argv[i + 1]);
    }
    if (!strcmp (argv[i], "-ceiling_temp") && i < argc) {
      exchgpara->ceiling_temp = atof (argv[i + 1]);
    }
    if (!strcmp (argv[i], "-nt") && i < argc) {
      int num_temp = atoi (argv[i + 1]);
      if (num_temp == 0) {
	cout << "number of temperature cannot set to be zero" << endl;
	cout << "docking exiting ..." << endl;
	exit (1);
      }
      else if (num_temp == 1) {
	exchgpara->num_temp = num_temp;
      }
      else {
	if ((num_temp <= MAXTMP) && (num_temp > 1))
	  exchgpara->num_temp = num_temp;
	else {
	  cout << "setting number of temperatures exceeds MAXTMP" << endl;
	  cout << "try modifying MAXTMP in size.h and compile again" << endl;
	  cout << "docking exiting ..." << endl;
	  exit (1);
	}
      }
    }

    // trace file
    if (!strcmp (argv[i], "-d") && i < argc) {
      inputfiles->hdf5_out_path = argv[i + 1];
      hdf_on = true;
    }
    
    // trace file
    if (!strcmp (argv[i], "-tr") && i < argc) {
      inputfiles->trace_file.path = argv[i + 1];
    }

    // ligand conformation sdf file corresponding to the trace file
    if (!strcmp (argv[i], "-lc") && i < argc) {
      inputfiles->lig_file.conf_path = argv[i + 1];
    }

    // move scale

    if (!strcmp (argv[i], "-t") && i < argc) {
      ts = atof (argv[i + 1]);
    }
    if (!strcmp (argv[i], "-r") && i < argc) {
      rs = atof (argv[i + 1]);
    }
  }

  if (!protein_on) {
    cout << "Provide target protein structure" << endl;
    exit (EXIT_FAILURE);
  }

  if (!compounds_on) {
    cout << "Provide compound library in SD format" << endl;
    exit (EXIT_FAILURE);
  }

  if (!lhm_one) {
    cout << "Provide LHM potentials" << endl;
    exit (EXIT_FAILURE);
  }

  mcpara->move_scale[0] = ts;
  mcpara->move_scale[1] = ts;
  mcpara->move_scale[2] = ts;
  mcpara->move_scale[3] = rs;
  mcpara->move_scale[4] = rs;
  mcpara->move_scale[5] = rs;


#if IS_OUTPUT == 1
  char mydir[MAXSTRINGLENG];

  char mystime[MAXSTRINGLENG];
  time_t mytime = time (NULL);
  struct tm *mylocaltime;
  mylocaltime = localtime (&mytime);
  strftime (mystime, MAXSTRINGLENG, "%Y%m%d_%H%M%S", mylocaltime);

  sprintf (mydir, "output_%s", mystime);
  strcpy (mcpara->outputdir, mydir);
  MakeDir (mydir);

  // prefix of the file name
  const char h5file[MAXSTRINGLENG] = "a";
  if (hdf_on == false)
    strcpy (mcpara->outputfile, h5file);
  else
    strcpy (mcpara->outputfile, inputfiles->hdf5_out_path.c_str());
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
      dst->t[residue] = src->t[residue];
      dst->c[residue] = src->c[residue];
      dst->n[residue] = src->n[residue];
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
void
InitLigCoord (Ligand * lig, const ComplexSize complexsize)
{
  for (int i = 0; i < complexsize.n_lig; ++i) {
    Ligand *mylig = &lig[i];

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
PlaceLigand (Ligand* mylig, float* movematrix_new)
{
  float rot[3][3];
  
  const float s1 = sinf (movematrix_new[3]);
  const float c1 = cosf (movematrix_new[3]);
  const float s2 = sinf (movematrix_new[4]);
  const float c2 = cosf (movematrix_new[4]);
  const float s3 = sinf (movematrix_new[5]);
  const float c3 = cosf (movematrix_new[5]);

  rot[0][0] = c1 * c2;
  rot[0][1] = c1 * s2 * s3 - c3 * s1;
  rot[0][2] = s1 * s3 + c1 * c3 * s2;
  rot[1][0] = c2 * s1;
  rot[1][1] = c1 * c3 + s1 * s2 * s3;
  rot[1][2] = c3 * s1 * s2 - c1 * s3;
  rot[2][0] = -1 * s2;
  rot[2][1] = c2 * s3;
  rot[2][2] = c2 * c3;

  LigCoord *coord_new = &mylig->coord_new;
  LigCoord *coord_orig = &mylig->coord_orig;


  const float cx = coord_orig->center[0];
  const float cy = coord_orig->center[1];
  const float cz = coord_orig->center[2];
  
  // iterate through all ligand residues
  // rotation and translation, and apply coordinate system transformation
  int lna = mylig->lna;
  for (int l = 0; l < lna; l += 1) {
    float x = coord_orig->x[l];
    float y = coord_orig->y[l];
    float z = coord_orig->z[l];
    coord_new->x[l] = rot[0][0] * x + rot[0][1] * y + rot[0][2] * z + movematrix_new[0] + cx;
    coord_new->y[l] = rot[1][0] * x + rot[1][1] * y + rot[1][2] * z + movematrix_new[1] + cy;
    coord_new->z[l] = rot[2][0] * x + rot[2][1] * y + rot[2][2] * z + movematrix_new[2] + cz;
  }
  
  
  for (int i = 0; i < 3; ++i) { 
    coord_new->center[i] = coord_orig->center[i] + movematrix_new[i];
  }
}

list < string > 
replaceLigandCoords(LigandFile * lig_file, Ligand * lig){
  list < string > l1_sdf;
  list < string >::iterator i1_sdf;
  string line1;
  ifstream compounds_file(lig_file->path.c_str());
  if (!compounds_file.is_open()) {
    cout << "cannot open ligand file" << endl;
    cout << "Cannot open " << lig_file->path << endl;
    exit(EXIT_FAILURE);
  }

  while (getline(compounds_file, line1))
    l1_sdf.push_back(line1);
  
  compounds_file.close();
  
  list < string > new_sdf;
  int line_num = 0;
  int lna = lig->lna;
  int atom_num = 0;
  char xyz[100];
  const LigCoord *mycoord = &lig->coord_new;
  
  for (i1_sdf = l1_sdf.begin(); i1_sdf != l1_sdf.end(); line_num++, i1_sdf++){
    string old_line = *i1_sdf;
    if (line_num > 3 && line_num < (4 + lna)) {
      // coordinates lines, replace
      float x = mycoord->x[atom_num];
      float y = mycoord->y[atom_num];
      float z = mycoord->z[atom_num];
      atom_num += 1;
      sprintf(xyz, "%10.4f%10.4f%10.4f", x, y, z);

      string rest = old_line.substr(30);
      string coords = string(xyz);
      string new_line = coords + rest;

      new_sdf.push_back(new_line);
    }
    else {
      // use the original lines
      new_sdf.push_back(old_line);
    }
  }
  
  return new_sdf;
}

void
PrintLigCoord2File(const Ligand * lig, string ofn)
{
  const LigCoord *mycoord = &lig->coord_new;
  
  ofstream myfile;
  myfile.open (ofn.c_str());

  myfile << "lna:\t" << lig->lna << endl;
  myfile << "@<BEGIN>ATOM\n";

  const int lna = lig->lna;
  myfile.precision(4);
  myfile << fixed;
  for (int i = 0; i < lna; ++i) {
    myfile << mycoord->x[i] << "  ";
    myfile << mycoord->y[i] << "  ";
    myfile << mycoord->z[i] << "  " << endl;
  }

  myfile << "@<END>ATOM" << endl;

  myfile.close();
}



void
CopyProteinResidue (const Protein0 * src, Protein * dst, const int residue_src,
		    const int residue_dst, const EnePara0 * enepara0)
{
  dst->x[residue_dst] = src->x[residue_src];
  dst->y[residue_dst] = src->y[residue_src];
  dst->z[residue_dst] = src->z[residue_src];

  dst->t[residue_dst] = src->t[residue_src];
  dst->c[residue_dst] = src->c[residue_src];

  int t = src->t[residue_src];
  int d = src->d[residue_src];
  int dt = t == 0 ? src->d[residue_src] + 30 : t;
  dst->ele[residue_dst] = enepara0->ele[dt];

  dst->seq3r[residue_dst] = src->seq3[src->r[residue_src]];

  dst->c0_and_d12_or_c2[residue_dst] =
    (src->c[residue_src] == 0 && src->d[residue_src] == 12) || (src->c[residue_src] == 2);

  dst->hpp[residue_dst] = enepara0->hpp[d];

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
    kde->x[i] = kde0->x[i];
    kde->y[i] = kde0->y[i];
    kde->z[i] = kde0->z[i];
    kde->t[i] = kde0->t[i];
  }
  kde->pnk = kde0->pnk;
}

void
OptimizeMcs (const Mcs0 * mcs0, Mcs * mcs, const ComplexSize complexsize)
{
  // pos
  for (int i = 0; i < complexsize.pos; ++i) {
    mcs[i].tcc = mcs0[i].tcc;

    // n_mcs
    for (int j = 0; j < MAXMCS; ++j) {
      mcs[i].x[j] = mcs0[i].x[j];
      mcs[i].y[j] = mcs0[i].y[j];
      mcs[i].z[j] = mcs0[i].z[j];
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
      enepara->p1a[i][j] = 2.0f * enepara0->vdw[j][i][1] * powf (tmp, 9.0f);
      enepara->p2a[i][j] = 3.0f * enepara0->vdw[j][i][1] * powf (tmp, 6.0f);
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
      enepara->pmf0[i][j] = enepara0->pmf[j][i][0];
      enepara->pmf1[i][j] = enepara0->pmf[j][i][1];
      enepara->hdb0[i][j] = enepara0->hdb[j][i][0];
      enepara->hdb1[i][j] = 1.0f / enepara0->hdb[j][i][1];
    }
  }

  for (int i = 0; i < MAXTP4; ++i)
    enepara->hpp[i] = enepara0->hpp[i];
  for (int i = 0; i < MAXTP2; ++i) {
    enepara->hpl0[i] = enepara0->hpl[i][0];
    enepara->hpl1[i] = enepara0->hpl[i][1];
    enepara->hpl2[i] = logf (1.0f / (enepara->hpl1[i] * sqrt_2_pi));
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
  printf ("Setting up temperature replicas\n");
  printf ("================================================================================\n");

  int num_temp = exchgpara->num_temp;
  float floor_temp = exchgpara->floor_temp;
  float ceiling_temp = exchgpara->ceiling_temp;

  if (num_temp == 1) {
    float my_temp = floor_temp;
    float beta = 1.0f / (BOLTZMANN_CONST * my_temp);

    for (int i = 0; i < num_temp; i++) {
      printf ("temp # %d\t\t\t%f\n", i, my_temp);

      temp[i].order = i;
      temp[i].minus_beta = 0.0f - beta;
    }
  }
  else {
    const float temp_ratio = exp (log (ceiling_temp / floor_temp) / (float) (num_temp - 1));
    float a = floor_temp;
    for (int i = 0; i < num_temp; i++) {
      float my_temp = a;
      float my_beta = 1.0f / (BOLTZMANN_CONST * my_temp);
      printf ("temp # %d\t\t\t%f\n", i, my_temp);

      temp[i].order = i;
      temp[i].minus_beta = 0.0f - my_beta;

      a *= temp_ratio;
    }
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
  mclog->ac_mc = 0;

  for (int i = 0; i < MAXREP; ++i)
    mclog->acs_mc[i] = 0;
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
    std::cout << step << "," << energy->e[9] << ","	// total
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
PrintStepTrack (const LigRecordSingleStep * step_record, const int step, const int arg)
{

  char names[11][30] = {
    "step",
    "prt_conf",
    "lig_conf",
    "temp_idx",
    "total",
    "x",
    "y",
    "z",
    "theta",
    "gama",
    "fi"
  };

  int a = arg & 0x1;
  int b = (arg >> 1) & 0x1;

  if (a == 1) {
    for (int i = 0; i < 11; ++i)
      printf (",%s", names[i]);
    printf ("\n");
  }

  if (b == 1) {
    printf (",%d", step_record->step);
    printf (",%d", step_record->replica.idx_prt);
    printf (",%d", step_record->replica.idx_lig);
    printf (",%d", step_record->replica.idx_tmp);
    printf (",%f", step_record->energy.e[MAXWEI - 1]);
    for (int i = 0; i < 6; i++)
      printf (",%f", step_record->movematrix[i]);
    printf ("\n");
  }
}

void
PrintCsv (const Energy * energy, const int idx_rep, const int step, const int arg)
{

  char names[MAXWEI][8] = {
    "vdw", // 0
    "ele", // 1
    "pmf", // 2
    "psp", // 3
    "hdb", // 4
    "hpc", // 5
    "kde", // 6
    "lhm", // 7
    "dst", // 8
    "total" // 9
  };

  int a = arg & 0x1;
  int b = (arg >> 1) & 0x1;

  if (a == 1) {
    printf ("rep step");
    for (int i = 0; i < MAXWEI; ++i)
      printf (" %s", names[i]);
    printf ("\n");
  }

  if (b == 1) {
    printf ("%d %d", idx_rep,  step);
    for (int i = 0; i < MAXWEI; ++i)
      // printf (" %+14.10f"  energy->e[i]);
      printf (" %.4f",  energy->e[i]);
    printf ("\n");
  }
}

void
PrintEnergy2 (const Energy * energy, const int idx_rep, const int step, const int arg)
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
  printf ("\n");
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


void
PrintTrack (LigRecord * ligrecord, int steps_per_dump, int replica,
	    int iter_begin, int iter_end, int arg)
{
  // print title
  // PrintEnergy2 (NULL, NULL, NULL, 1);
  PrintStepTrack (NULL, NULL, 1);

  for (int s = iter_begin; s <= iter_end; ++s) {
    const LigRecordSingleStep *myrecord = &ligrecord[replica].step[s];
    // PrintEnergy2 (&myrecord->energy, replica, myrecord->step, arg);
    PrintStepTrack (myrecord, myrecord->step, arg);
  }

}

// arg = 1      print title
// arg = 2      print content
// arg = 3      print all

void
PrintLigRecord (LigRecord * ligrecord, int steps_per_dump, int replica,
		int iter_begin, int iter_end, int arg)
{
  // print title
  // PrintEnergy2 (NULL, NULL, NULL, 1);
  PrintCsv (NULL, NULL, NULL, 1);

  for (int s = iter_begin; s <= iter_end; ++s) {
    const LigRecordSingleStep *myrecord = &ligrecord[replica].step[s];
    // PrintEnergy2 (&myrecord->energy, replica, myrecord->step, arg);
    PrintCsv (&myrecord->energy, replica, myrecord->step, arg);
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
PrintRepRecord2 (LigRecord * ligrecord, ComplexSize complexsize,
		 int steps_per_dump, int idx_prt, int idx_lig,
		 int iter_begin, int iter_end, int arg)
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
	complexsize.n_tmp * complexsize.n_lig * idx_prt + complexsize.n_lig * t + idx_lig;
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
  printf ("center:\t\t%+10.6f\t%+10.6f\t%+10.6f\n", mycoord->center[0], mycoord->center[1],
	  mycoord->center[2]);
  printf ("lna:\t\t%d\n", lig->lna);

  printf ("x \t\ty \t\tz \t\tc \t\t t \t n \tindex\n");
  printf ("-----------------------------------------------\n");
  const int lna = lig->lna;
  for (int i = 0; i < lna; ++i) {
    printf ("%+10.6f\t", mycoord->x[i]);
    printf ("%+10.6f\t", mycoord->y[i]);
    printf ("%+10.6f\t", mycoord->z[i]);
    printf ("%+10.6f\t", lig->c[i]);
    printf ("%2d\t", lig->t[i]);
    printf ("%2d\t", lig->n[i]);
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
    printf ("%+10.6f\t", prt->x[i]);
    printf ("%+10.6f\t", prt->y[i]);
    printf ("%+10.6f\t", prt->z[i]);
    printf ("%2d\t", prt->t[i]);
    printf ("%2d\t", prt->c[i]);
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

  // inputs and outputs
  printf ("================================================================================\n");
  printf ("Inputs and Outputs\n");
  printf ("================================================================================\n");
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

  printf ("steps_per_dump\t\t\t%d\n", mcpara->steps_per_dump);

  const size_t ligrecord_sz = sizeof (LigRecord) * complexsize->n_rep;
  printf ("per dump record size:\t\t%.3f MB\n", (float) ligrecord_sz / 1024 / 1024);

  printf ("================================================================================\n");

  // Replica Exchange Monte carlo parameters
  printf ("Replica Exchange Monte Carlo parameters\n");
  printf ("================================================================================\n");
  printf ("steps_total\t\t\t%d\n", mcpara->steps_total);
  printf ("steps_per_dump\t\t\t%d\n", mcpara->steps_per_dump);
  printf ("steps_per_exchange\t\t%d\n", mcpara->steps_per_exchange);

  printf ("translational scale\t\t");
  for (int i = 0; i < 3; ++i)
    printf ("%.8f ", mcpara->move_scale[i]);
  printf ("\n");

  printf ("rotational scale\t\t");
  for (int i = 3; i < 6; ++i)
    printf ("%.8f ", mcpara->move_scale[i]);
  printf ("\n");

  printf ("ligand conformations\t\t%d\n", complexsize->n_lig);
  printf ("prt conformations\t\t%d\n", complexsize->n_prt);
  printf ("temperatures\t\t\t%d\n", complexsize->n_tmp);
  printf ("replica ensembles\t\t%d\n", complexsize->n_rep);

  printf ("size_lig\t\t\t%d\n", complexsize->lna);
  printf ("size_prt\t\t\t%d\n", complexsize->pnp);
  printf ("size_pnk\t\t\t%d\n", complexsize->pnk);
  printf ("size_mcs\t\t\t%d\n", complexsize->pos);


  printf ("AR of MC \t\t\t%d / %d \t%f\n",
	  mclog->ac_mc,
	  mcpara->steps_total * complexsize->n_rep,
	  (float) mclog->ac_mc / (mcpara->steps_total * complexsize->n_rep));


#if 0
  for (int t = 0; t < complexsize->n_tmp; ++t) {
    const int myreplica = complexsize->n_lig * t;
    printf ("AR of temperature[%d]=%f \t %d / %d \t%f\n",
	    t,
	    temp[t].t,
	    mclog->acs_mc[myreplica],
	    mcpara->steps_total, (float) mclog->acs_mc[myreplica] / mcpara->steps_total);
  }
#endif

#if 0

  for (int i = 0; i < complexsize->n_prt; ++i) {
    for (int j = 0; j < complexsize->n_tmp; ++j) {
      for (int k = 0; k < complexsize->n_lig; ++k) {
	const int flatten_addr =
	  complexsize->n_tmp * complexsize->n_lig * i + complexsize->n_lig * j + k;
	printf ("AR of %4d temperature[%d]=%f \t %d / %d \t%f\n", flatten_addr, j, temp[j].t,
		mclog->acs_mc[flatten_addr], mcpara->steps_total,
		(float) mclog->acs_mc[flatten_addr] / mcpara->steps_total);
      }
    }
  }
#endif

  printf ("AR of temp exchange \t\t%d / %d \t%f\n",
	  mclog->ac_temp_exchg,
	  mcpara->steps_total / mcpara->steps_per_exchange * complexsize->n_rep,
	  (float) mclog->ac_temp_exchg / (mcpara->steps_total / mcpara->steps_per_exchange *
					  complexsize->n_rep));

  printf ("================================================================================\n");

  // performance
  printf ("Performance\n");
  printf ("================================================================================\n");

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
  printf ("================================================================================\n");


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


vector < string >
splitByWhiteSpace(string s) {
  vector < string > tokens;
  istringstream ss (s);
  while (!ss.eof())
    {
      string x;
      getline(ss, x, ' ');
      tokens.push_back(x);
    }

  return tokens;
}


// return 1 if two movevector are the same else return 0
int
sameVector(float *v1, float *v2)
{
  float threshold = 0.01;
  for (int i = 0; i < 6; i++) {
    float diff = fabs(v1[i] - v2[i]);
    if (diff > threshold)
      return 0;
  }
  return 1;
}

int
checkRedundancy(vector < LigRecordSingleStep > &records,
                int idx_rep,
                LigRecord * ligrecord)
{
  // rare array of float used to be compared
  // expext no initial move-vectors be the same as this one
  float current_matrix[6] = {3.10, 1.6, 4.8, 1.2, 0.03, 0.08};
  
  for (int i = 0; i < STEPS_PER_DUMP; i++) {
    LigRecordSingleStep *myrecord = &ligrecord[idx_rep].step[i];
    float *movematrix = myrecord->movematrix;

    if (!sameVector(current_matrix, movematrix))
      {
        LigRecordSingleStep rec;
        memcpy(&rec, myrecord, sizeof(LigRecordSingleStep));
        records.push_back(rec);

        // copied for the next comparison
        memcpy(current_matrix, movematrix, sizeof(current_matrix));
      }
  }

  return records.size();
}
