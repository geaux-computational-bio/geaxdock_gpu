//#ifndef  _MKDIR_H_
//#define  _MKDIR_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>


void
MakeDir (const char *mydir)
{
  struct stat st;

  // new dir
  if (stat (mydir, &st) != 0) {
    mkdir (mydir, 0777);
  }

  // the dir exists
  else {
    //printf ("program terminated: output dir conflict\n");
    //exit (1);

    printf ("the output directory %s exists, force overwriting? y/n ", mydir);
    char input;
    scanf ("%c", &input);
    if (input == 'n') {
      printf ("program terminated\n");
      exit (1);
    }
  }

}


//#endif


