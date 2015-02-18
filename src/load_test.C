#include <stdio.h>

#include "size.h"
#include "dock.h"
#include "load.h"

#include "gtest/gtest.h"

TEST (Load_Ligand, 1 a07C1_sdf)
{
  Ligand0 *lig0 = new Ligand0[MAXEN2];
  printf ("testing ligand loading ...\n");

  delete[]lig0;
}

