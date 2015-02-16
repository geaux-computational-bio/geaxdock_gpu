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


#include "size.h"
#include "data.h"

using namespace std;

std::string three2oneS (std::string resnam1)
{
  if (resnam1 == "ALA") {
    return "A";
  }
  else if (resnam1 == "CYS") {
    return "C";
  }
  else if (resnam1 == "ASP") {
    return "D";
  }
  else if (resnam1 == "GLU") {
    return "E";
  }
  else if (resnam1 == "PHE") {
    return "F";
  }
  else if (resnam1 == "GLY") {
    return "G";
  }
  else if (resnam1 == "HIS") {
    return "H";
  }
  else if (resnam1 == "ILE") {
    return "I";
  }
  else if (resnam1 == "LYS") {
    return "K";
  }
  else if (resnam1 == "LEU") {
    return "L";
  }
  else if (resnam1 == "MET") {
    return "M";
  }
  else if (resnam1 == "ASN") {
    return "N";
  }
  else if (resnam1 == "PRO") {
    return "P";
  }
  else if (resnam1 == "GLN") {
    return "Q";
  }
  else if (resnam1 == "ARG") {
    return "R";
  }
  else if (resnam1 == "SER") {
    return "S";
  }
  else if (resnam1 == "THR") {
    return "T";
  }
  else if (resnam1 == "VAL") {
    return "V";
  }
  else if (resnam1 == "TRP") {
    return "W";
  }
  else if (resnam1 == "TYR") {
    return "Y";
  }

  else {
    cout << "Unknown residue passed to three2one: " << resnam1 << endl <<
      endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}

char
three2oneC (std::string resnam1)
{
  if (resnam1 == "ALA") {
    return 'A';
  }
  else if (resnam1 == "CYS") {
    return 'C';
  }
  else if (resnam1 == "ASP") {
    return 'D';
  }
  else if (resnam1 == "GLU") {
    return 'E';
  }
  else if (resnam1 == "PHE") {
    return 'F';
  }
  else if (resnam1 == "GLY") {
    return 'G';
  }
  else if (resnam1 == "HIS") {
    return 'H';
  }
  else if (resnam1 == "ILE") {
    return 'I';
  }
  else if (resnam1 == "LYS") {
    return 'K';
  }
  else if (resnam1 == "LEU") {
    return 'L';
  }
  else if (resnam1 == "MET") {
    return 'M';
  }
  else if (resnam1 == "ASN") {
    return 'N';
  }
  else if (resnam1 == "PRO") {
    return 'P';
  }
  else if (resnam1 == "GLN") {
    return 'Q';
  }
  else if (resnam1 == "ARG") {
    return 'R';
  }
  else if (resnam1 == "SER") {
    return 'S';
  }
  else if (resnam1 == "THR") {
    return 'T';
  }
  else if (resnam1 == "VAL") {
    return 'V';
  }
  else if (resnam1 == "TRP") {
    return 'W';
  }
  else if (resnam1 == "TYR") {
    return 'Y';
  }

  else {
    cout << "Unknown residue passed to three2one: " << resnam1 << endl <<
      endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}

std::string one2three (std::string resnam1)
{
  if (resnam1 == "A") {
    return "ALA";
  }
  else if (resnam1 == "C") {
    return "CYS";
  }
  else if (resnam1 == "D") {
    return "ASP";
  }
  else if (resnam1 == "E") {
    return "GLU";
  }
  else if (resnam1 == "F") {
    return "PHE";
  }
  else if (resnam1 == "G") {
    return "GLY";
  }
  else if (resnam1 == "H") {
    return "HIS";
  }
  else if (resnam1 == "I") {
    return "ILE";
  }
  else if (resnam1 == "K") {
    return "LYS";
  }
  else if (resnam1 == "L") {
    return "LEU";
  }
  else if (resnam1 == "M") {
    return "MET";
  }
  else if (resnam1 == "N") {
    return "ASN";
  }
  else if (resnam1 == "P") {
    return "PRO";
  }
  else if (resnam1 == "Q") {
    return "GLN";
  }
  else if (resnam1 == "R") {
    return "ARG";
  }
  else if (resnam1 == "S") {
    return "SER";
  }
  else if (resnam1 == "T") {
    return "THR";
  }
  else if (resnam1 == "V") {
    return "VAL";
  }
  else if (resnam1 == "W") {
    return "TRP";
  }
  else if (resnam1 == "Y") {
    return "TYR";
  }

  else {
    cout << "Unknown residue passed to one2three: " << resnam1 << endl <<
      endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}

int
getResCode (std::string r_name)
{
  if (r_name == "ALA")
    return 0;
  else if (r_name == "ARG")
    return 1;
  else if (r_name == "CYS")
    return 2;
  else if (r_name == "PRO")
    return 3;
  else if (r_name == "LEU")
    return 4;
  else if (r_name == "ILE")
    return 5;
  else if (r_name == "GLU")
    return 6;
  else if (r_name == "HIS")
    return 7;
  else if (r_name == "TRP")
    return 8;
  else if (r_name == "THR")
    return 9;
  else if (r_name == "MET")
    return 10;
  else if (r_name == "ASP")
    return 11;
  else if (r_name == "GLY")
    return 12;
  else if (r_name == "PHE")
    return 13;
  else if (r_name == "LYS")
    return 14;
  else if (r_name == "GLN")
    return 15;
  else if (r_name == "SER")
    return 16;
  else if (r_name == "ASN")
    return 17;
  else if (r_name == "TYR")
    return 18;
  else if (r_name == "VAL")
    return 19;
  else {
    cout << "Unknown residue: " << r_name << endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}

int
getResCodeOne (std::string r_name)
{
  if (r_name == "A")
    return 0;
  else if (r_name == "R")
    return 1;
  else if (r_name == "C")
    return 2;
  else if (r_name == "P")
    return 3;
  else if (r_name == "L")
    return 4;
  else if (r_name == "I")
    return 5;
  else if (r_name == "E")
    return 6;
  else if (r_name == "H")
    return 7;
  else if (r_name == "W")
    return 8;
  else if (r_name == "T")
    return 9;
  else if (r_name == "M")
    return 10;
  else if (r_name == "D")
    return 11;
  else if (r_name == "G")
    return 12;
  else if (r_name == "F")
    return 13;
  else if (r_name == "K")
    return 14;
  else if (r_name == "Q")
    return 15;
  else if (r_name == "S")
    return 16;
  else if (r_name == "N")
    return 17;
  else if (r_name == "Y")
    return 18;
  else if (r_name == "V")
    return 19;
  else {
    cout << "Unknown residue: " << r_name << endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}

std::string getResName (int r_code)
{
  if (r_code == 0)
    return "ALA";
  else if (r_code == 1)
    return "ARG";
  else if (r_code == 2)
    return "CYS";
  else if (r_code == 3)
    return "PRO";
  else if (r_code == 4)
    return "LEU";
  else if (r_code == 5)
    return "ILE";
  else if (r_code == 6)
    return "GLU";
  else if (r_code == 7)
    return "HIS";
  else if (r_code == 8)
    return "TRP";
  else if (r_code == 9)
    return "THR";
  else if (r_code == 10)
    return "MET";
  else if (r_code == 11)
    return "ASP";
  else if (r_code == 12)
    return "GLY";
  else if (r_code == 13)
    return "PHE";
  else if (r_code == 14)
    return "LYS";
  else if (r_code == 15)
    return "GLN";
  else if (r_code == 16)
    return "SER";
  else if (r_code == 17)
    return "ASN";
  else if (r_code == 18)
    return "TYR";
  else if (r_code == 19)
    return "VAL";
  else {
    cout << "Unknown residue code: " << r_code << endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}

int
getPntCode (std::string r_name)
{
  if (r_name == "CA")
    return 0;
  else if (r_name == "PP")
    return 1;
  else if (r_name == "A1")
    return 2;
  else if (r_name == "C1")
    return 3;
  else if (r_name == "D1")
    return 4;
  else if (r_name == "E1")
    return 5;
  else if (r_name == "E2")
    return 6;
  else if (r_name == "F1")
    return 7;
  else if (r_name == "F2")
    return 8;
  else if (r_name == "H1")
    return 9;
  else if (r_name == "H2")
    return 10;
  else if (r_name == "I1")
    return 11;
  else if (r_name == "K1")
    return 12;
  else if (r_name == "K2")
    return 13;
  else if (r_name == "L1")
    return 14;
  else if (r_name == "M1")
    return 15;
  else if (r_name == "M2")
    return 16;
  else if (r_name == "N1")
    return 17;
  else if (r_name == "P1")
    return 18;
  else if (r_name == "Q1")
    return 19;
  else if (r_name == "Q2")
    return 20;
  else if (r_name == "R1")
    return 21;
  else if (r_name == "R2")
    return 22;
  else if (r_name == "S1")
    return 23;
  else if (r_name == "T1")
    return 24;
  else if (r_name == "V1")
    return 25;
  else if (r_name == "W1")
    return 26;
  else if (r_name == "W2")
    return 27;
  else if (r_name == "Y1")
    return 28;
  else if (r_name == "Y2")
    return 29;
  else {
    cout << "Unknown effective point: " << r_name << endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}

std::string getPntName (int r_code)
{
  if (r_code == 0)
    return "CA";
  else if (r_code == 1)
    return "PP";
  else if (r_code == 2)
    return "A1";
  else if (r_code == 3)
    return "C1";
  else if (r_code == 4)
    return "D1";
  else if (r_code == 5)
    return "E1";
  else if (r_code == 6)
    return "E2";
  else if (r_code == 7)
    return "F1";
  else if (r_code == 8)
    return "F2";
  else if (r_code == 9)
    return "H1";
  else if (r_code == 10)
    return "H2";
  else if (r_code == 11)
    return "I1";
  else if (r_code == 12)
    return "K1";
  else if (r_code == 13)
    return "K2";
  else if (r_code == 14)
    return "L1";
  else if (r_code == 15)
    return "M1";
  else if (r_code == 16)
    return "M2";
  else if (r_code == 17)
    return "N1";
  else if (r_code == 18)
    return "P1";
  else if (r_code == 19)
    return "Q1";
  else if (r_code == 20)
    return "Q2";
  else if (r_code == 21)
    return "R1";
  else if (r_code == 22)
    return "R2";
  else if (r_code == 23)
    return "S1";
  else if (r_code == 24)
    return "T1";
  else if (r_code == 25)
    return "V1";
  else if (r_code == 26)
    return "W1";
  else if (r_code == 27)
    return "W2";
  else if (r_code == 28)
    return "Y1";
  else if (r_code == 29)
    return "Y2";
  else {
    cout << "Unknown effective point code: " << r_code << endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}

int
getLigCode (std::string r_name)
{
  if (r_name == "Br")
    return 0;
  else if (r_name == "C.1")
    return 1;
  else if (r_name == "C.2")
    return 2;
  else if (r_name == "C.3")
    return 3;
  else if (r_name == "C.ar")
    return 4;
  else if (r_name == "C.cat")
    return 5;
  else if (r_name == "Cl")
    return 6;
  else if (r_name == "F")
    return 7;
  else if (r_name == "I")
    return 8;
  else if (r_name == "N.1")
    return 9;
  else if (r_name == "N.2")
    return 10;
  else if (r_name == "N.3")
    return 11;
  else if (r_name == "N.4")
    return 12;
  else if (r_name == "N.am")
    return 13;
  else if (r_name == "N.ar")
    return 14;
  else if (r_name == "N.pl3")
    return 15;
  else if (r_name == "O.2")
    return 16;
  else if (r_name == "O.3")
    return 17;
  else if (r_name == "O.co2")
    return 18;
  else if (r_name == "P.3")
    return 19;
  else if (r_name == "S.2")
    return 20;
  else if (r_name == "S.3")
    return 21;
  else if (r_name == "S.O")
    return 22;
  else if (r_name == "S.O2")
    return 23;
  else {
    return BADKDE;
  }

  exit (EXIT_FAILURE);
}

std::string getLigName (int r_code)
{
  if (r_code == 0)
    return "Br";
  else if (r_code == 1)
    return "C.1";
  else if (r_code == 2)
    return "C.2";
  else if (r_code == 3)
    return "C.3";
  else if (r_code == 4)
    return "C.ar";
  else if (r_code == 5)
    return "C.cat";
  else if (r_code == 6)
    return "Cl";
  else if (r_code == 7)
    return "F";
  else if (r_code == 8)
    return "I";
  else if (r_code == 9)
    return "N.1";
  else if (r_code == 10)
    return "N.2";
  else if (r_code == 11)
    return "N.3";
  else if (r_code == 12)
    return "N.4";
  else if (r_code == 13)
    return "N.am";
  else if (r_code == 14)
    return "N.ar";
  else if (r_code == 15)
    return "N.pl3";
  else if (r_code == 16)
    return "O.2";
  else if (r_code == 17)
    return "O.3";
  else if (r_code == 18)
    return "O.co2";
  else if (r_code == 19)
    return "P.3";
  else if (r_code == 20)
    return "S.2";
  else if (r_code == 21)
    return "S.3";
  else if (r_code == 22)
    return "S.O";
  else if (r_code == 23)
    return "S.O2";
  else {
    cout << "Unknown ligand atom code: " << r_code << endl;
    exit (EXIT_FAILURE);
  }

  exit (EXIT_FAILURE);
}
