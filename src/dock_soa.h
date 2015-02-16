#ifndef  DOCK_SOA_H
#define  DOCK_SOA_H

// struct of array



struct Ligand
{
  // coord_center is always under lab system
  // coord_xyz for "orig" is under ligand_ceter system
  // coord_xyz for "new" is under lab system
  LigCoord coord_orig;
  LigCoord coord_new;

  // translation x y z, rotation x y z
  float movematrix_old[6];       // old matrix
  float movematrix_new[6];       // trail matrix

  Energy energy_old;		//                                      used
  Energy energy_new;		//                                      used

  int t[MAXLIG];		// atom type                            used
  float c[MAXLIG];		// atom charge                          used
  int n[MAXLIG];		// atom number                          used
                                // n == index???

  int lna;			// number of ligand atoms               used


  // confusion matrix
  int native_confusion_matx[MAXLIG][MAXPRO];
  int decoy_confusion_matx[MAXLIG][MAXPRO];
};





struct Protein
{
  float x[MAXPRO];		// residue x coord
  float y[MAXPRO];		// residue y coord
  float z[MAXPRO];		// residue z coord

  int t[MAXPRO];		// effective point type
  int c[MAXPRO];		// effective point class


  float ele[MAXPRO];            // dt = prt->t[i] == 0 ? prt->d[i] + 30 : prt->t[i];
                                // enepara->ele[dt]

  int seq3r[MAXPRO];            // prt->seq3r[i] == prt->seq3[prt->r[i]];
  int c0_and_d12_or_c2[MAXPRO]; // (prt_c == 0 && prt_d == 12)) || (prt_c == 2)
  float hpp[MAXPRO];            // enepara->hpp[prt->d[i]]

  int pnp;			// number of protein effective points

  float pocket_center[3];
};




struct Psp
{
  float psp[MAXLIG][MAXPRO];                                           // replaced
  //float sparse_psp[MAXLIG][MAXPRO];                                    // replaced

};



struct Kde
{
  float x[MAXKDE];		// KDE x coord                          used
  float y[MAXKDE];		// KDE y coord                          used
  float z[MAXKDE];		// KDE z coord                          used
  int t[MAXKDE];		// KDE atom type                        used

  int pnk;			// number of kde points                 used
};




struct Mcs
{
  float x[MAXMCS];              //                         used, mcs->y[lig_n]
  float y[MAXMCS];              //                         used
  float z[MAXMCS];              //                         used

  float tcc;                    //                         used
};



struct EnePara
{
  // L-J
  float p1a[MAXTP2][MAXTP1];
  float p2a[MAXTP2][MAXTP1];
  float lj0, lj1;

  // electrostatic
  float el0;
  float el1;
  float ele[MAXTP3];
  float a1; // 4.0f - 3.0f * el0;
  float b1; // 2.0f * el0 - 3.0f;

  // contact
  float pmf0[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float pmf1[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float hdb0[MAXTP2][MAXTP1];  // interchange to [lig][prt]
  float hdb1[MAXTP2][MAXTP1];  // interchange to [lig][prt]

  // hydrophobic
  float hpp[MAXTP4];
  float hpl0[MAXTP2];
  float hpl1[MAXTP2];
  float hpl2[MAXTP2];

  // kde
  float kde2; // -0.5f / (kde * kde)
  float kde3; // powf (kde * sqrtf (2.0f * PI), 3.0f)

  // weights for energy terms
  float w[MAXWEI];
  float a_para[MAXWEI];         // the a parameter in normalization
  float b_para[MAXWEI];         // the b parameter in normalization
};



#endif // DOCK_SOA_H

