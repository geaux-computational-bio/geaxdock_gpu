// simulate the 1D array of "ligand record" as a 3D array "ligrecord[rep][steps]"
#define LIGRECORD(r, s) (ligrecord[steps_per_dump * (r) + (s)])
#define LIGRECORD_DC(r, s) (ligrecord_dc[steps_per_dump_dc * (r) + (s)])

