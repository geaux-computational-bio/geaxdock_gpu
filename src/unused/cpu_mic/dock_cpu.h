#ifndef  DOCK_CPU_H
#define  DOCK_CPU_H




struct LigPoint
{
  int t;		// atom type                            used
  float c;		// atom charge                          used
  int n;		// atom number                          used
};



struct Ligand
{
  // coord_center is always under lab system
  // coord_xyz for "orig" is under ligand_ceter system
  // coord_xyz for "new" is under lab system
  LigCoord coord_orig;
  LigCoord coord_new;


  LigPoint lig_point[MAXLIG];



  // translation x y z, rotation x y z
  float movematrix_old[6];       // old matrix
  float movematrix_new[6];       // trail matrix

  Energy energy_old;		//                                      used
  Energy energy_new;		//                                      used
  int lna;			// number of ligand atoms               used
};




struct PrtPoint
{
  float x;		// residue x coord
  float y;		// residue y coord
  float z;		// residue z coord

  int t;		// effective point type
  int c;		// effective point class


  float ele;            // dt = prt->t[i] == 0 ? prt->d[i] + 30 : prt->t[i];
                                // enepara->ele[dt]

  int seq3r;            // prt->seq3r[i] == prt->seq3[prt->r[i]];
  int c0_and_d12_or_c2; // (prt_c == 0 && prt_d == 12)) || (prt_c == 2)
  float hpp;            // enepara->hpp[prt->d[i]]
};




struct Protein
{
  PrtPoint prt_point[MAXPRO];
  int pnp;			// number of protein effective points
  float pocket_center[3];
};




struct Psp
{
  float psp[MAXLIG][MAXPRO];                                           // replaced
  //float sparse_psp[MAXLIG][MAXPRO];                                    // replaced

};


struct KdePoint
{
  float x;		// KDE x coord                          used
  float y;		// KDE y coord                          used
  float z;		// KDE z coord                          used
  int t;		// KDE atom type                        used
};


struct Kde
{
  KdePoint kde_point[MAXKDE];
  int pnk;			// number of kde points                 used
};



struct McsPoint
{
  float x;              //                         used, mcs->y[lig_n]
  float y;              //                         used
  float z;              //                         used
};


struct Mcs
{
  McsPoint mcs_point[MAXMCS];
  float tcc;                    //                         used
};



struct EnePara
{
  // L-J
  float pa[MAXTP2][MAXTP1][2];
  float lj0, lj1;

  // electrostatic
  float el0;
  float el1;
  float ele[MAXTP3];
  float a1; // 4.0f - 3.0f * el0;
  float b1; // 2.0f * el0 - 3.0f;

  // contact
  float pmf[MAXTP2][MAXTP1][2];  // interchange to [lig][prt]
  float hdb[MAXTP2][MAXTP1][3];  // interchange to [lig][prt]

  // hydrophobic
  float hpp[MAXTP4];
  float hpl[MAXTP2][3];

  // kde
  float kde2; // -0.5f / (kde * kde)
  float kde3; // powf (kde * sqrtf (2.0f * PI), 3.0f)

  // weights for energy terms

  float a_para[MAXWEI];         // the a parameter in normalization
  float b_para[MAXWEI];         // the b parameter in normalization
  float w[MAXWEI];

};



#endif // DOCK_CPU_H

