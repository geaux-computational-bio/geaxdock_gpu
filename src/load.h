#ifndef  LOAD_H
#define  LOAD_H

#include <vector>
#include <string>

#include "dock.h"

using namespace std;

void loadTrace (TraceFile *, float *);


void trimLigand (InputFiles * inputfiles, Ligand0 * lig);

// move the ligand to the protein pocket center
void moveLigand2PocketCenter(float *, Ligand0 *);

// move the ligand to its center frame
void moveLigand2ItsCenterFrame (Ligand0 *);

// write the property of the ligand, such as SMILES, does not wirte coordinates
void writeLigProperty( vector < string > , Ligand0 *);

// write properties of the atoms
void writeLigAtomProperty( vector < string > , Ligand0 *);

// write the coords to the default ligand
void writeDefaultLigAtomCoord (vector < string >, Ligand0 *);

// write the coords to the ensemble ligand
void writeEnsAtomCoord (string coord_line, Ligand0 * mylig);

// load ligand and return the total number of conformers
int loadOneLigand (vector < string >, Ligand0 *);

int getLigEnsembleTotal ( vector < string > );

vector < vector < string > > readLigandSections (string);
vector < string > getLigEnsembleCoords (vector < string > );
vector < float > getLigEnsembleRmsd (vector < string > );


void loadLigConf (LigandFile *);
void loadLigand (InputFiles *, Ligand0 *);
void loadLigand_bk (LigandFile *, Ligand0 *);

void loadPrtConf (ProteinFile *, Protein0 *);
void loadProtein (ProteinFile *, Protein0 *);

void loadPocketCenter (string, float *);
void loadLHM (LhmFile *, Psp0 *, Kde0 *, Mcs0 *);
void loadEnePara (EneParaFile *, EnePara0 *);

void loadWeight(WeightFile *, EnePara0 *);
void loadNorPara(NorParaFile *, EnePara0 *);  // load the normalization parameter a and b

vector < vector < float > > read2D(TraceFile *);


#endif // LOAD_H

