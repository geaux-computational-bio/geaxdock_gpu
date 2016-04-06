#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <string>
#include <cstdio>

#include "size.h"
#include "dock.h"
#include "run.h"
#include "util.h"
#include "toggle.h"
#include "load.h"
#include "stats.h"
#include "post_mc.h"
#include "boost/program_options.hpp"


namespace {
  const size_t ERROR_IN_COMMAND_LINE = 1;
  const size_t SUCCESS = 0;
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;
}

int main(int argc, char **argv) {
  // Banner ();

  try {
    std::string pdb_path, sdf_path, ff_path, id, para;

    McPara mcpara = McPara();
    ExchgPara exchgpara = ExchgPara();
    InputFiles inputfiles = InputFiles();

    // default values
    float ts = 0.02f, rs = 0.08f;
    exchgpara.floor_temp = 0.3f;
    exchgpara.ceiling_temp = 0.3f;
    exchgpara.num_temp = 1;
    mcpara.steps_total = STEPS_PER_DUMP;
    mcpara.steps_per_dump = STEPS_PER_DUMP;
    mcpara.steps_per_exchange = 10;
    inputfiles.enepara_file.path = "gpudocksm.ff";
    inputfiles.lig_file.molid = "MOLID";

    namespace po = boost::program_options;

    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "Print help messages")

      ("pdb,p", po::value<std::string>(&inputfiles.prt_file.path)->required(),"protein path (PDB)")
      ("sdf,l", po::value<std::string>(&inputfiles.lig_file.path)->required(),"ligand path (SDF)")
      ("ff,s", po::value<std::string>(&inputfiles.lhm_file.path)->required(), "force field path")
      ("para", po::value<std::string>(&inputfiles.enepara_file.path)->required(), "parameter file path")
      ("id,i", po::value<std::string>(&inputfiles.lhm_file.ligand_id)->required(), "complex id")
      ("csv,o", po::value<std::string>(&inputfiles.trace_file.path)->required(), "trajectories")

      ("nc", po::value<int>(&mcpara.steps_per_exchange), "")
      ("nt", po::value<int>(&exchgpara.num_temp), "number of temperatures")
      ("ceiling_temp", po::value<float>(&exchgpara.ceiling_temp), "ceiling temperature")
      ("floor_temp", po::value<float>(&exchgpara.floor_temp), "floor temperature")
      (",t", po::value<float>(&ts), "translational scale")
      (",r", po::value<float>(&rs), "rotational scale")
      ;

    mcpara.move_scale[0] = ts;
    mcpara.move_scale[1] = ts;
    mcpara.move_scale[2] = ts;
    mcpara.move_scale[3] = rs;
    mcpara.move_scale[4] = rs;
    mcpara.move_scale[5] = rs;

    po::variables_map vm;
    try {
      po::store(po::parse_command_line(argc, argv, desc), vm);
      if (vm.count("help")) {
        std::cout << "GeauxDock usage:" << std::endl << desc << std::endl;
        return SUCCESS;
      }

      po::notify(vm);
    }
    catch (po::error & e) {
      std::cerr << "Command line parse error: " << e.what() << std::endl
                << "GeauxDock will now exit" << std::endl;
      return ERROR_IN_COMMAND_LINE;
    }


    // run application
    srand (time (0));

    McLog *mclog = new McLog;

    // load into preliminary data structures
    Ligand0 *lig0 = new Ligand0[MAXEN2];
    Protein0 *prt0 = new Protein0[MAXEN1];
    Psp0 *psp0 = new Psp0;
    Kde0 *kde0 = new Kde0;
    Mcs0 *mcs0 = new Mcs0[MAXPOS];
    EnePara0 *enepara0 = new EnePara0;

    loadLigand (&inputfiles, lig0);
    loadProtein (&inputfiles.prt_file, prt0);
    loadLHM (&inputfiles.lhm_file, psp0, kde0, mcs0);
    loadEnePara (&inputfiles.enepara_file, enepara0);

    // sizes
    ComplexSize complexsize;
    complexsize.n_prt = inputfiles.prt_file.conf_total;	// number of protein conf
    complexsize.n_tmp = exchgpara.num_temp;	// number of temperature
    complexsize.n_lig = inputfiles.lig_file.conf_total;	// number of ligand conf
    complexsize.n_rep = complexsize.n_lig * complexsize.n_prt * complexsize.n_tmp;
    complexsize.lna = inputfiles.lig_file.lna;
    complexsize.pnp = inputfiles.prt_file.pnp;
    complexsize.pnk = kde0->pnk;
    complexsize.pos = inputfiles.lhm_file.pos;	// number of MCS positions

    // data structure optimizations
    Ligand *lig = new Ligand[complexsize.n_rep];
    Protein *prt = new Protein[complexsize.n_prt];
    Psp *psp = new Psp;
    Kde *kde = new Kde;
    Mcs *mcs = new Mcs[complexsize.pos];
    EnePara *enepara = new EnePara;
    Temp *temp = new Temp[complexsize.n_tmp];
    Replica *replica = new Replica[complexsize.n_rep];

    OptimizeLigand (lig0, lig, complexsize);
    OptimizeProtein (prt0, prt, enepara0, lig0, complexsize);
    OptimizePsp (psp0, psp, lig, prt);
    OptimizeKde (kde0, kde);
    OptimizeMcs (mcs0, mcs, complexsize);
    OptimizeEnepara (enepara0, enepara);

    delete[]lig0;
    delete[]prt0;
    delete[]psp0;
    delete[]kde0;
    delete[]mcs0;
    delete[]enepara0;

    // initialize system
    InitLigCoord (lig, complexsize);
    SetTemperature (temp, &exchgpara);
    // SetTemperature (temp, mcpara, complexsize);
    SetReplica (replica, lig, complexsize);
    SetMcLog (mclog);

    // debug
    //PrintDataSize (lig, prt, psp, kde, mcs, enepara);
    //PrintLigand (lig);
    //PrintProtein (prt);

    std::map<int, std::vector<LigRecordSingleStep> > multi_reps_records;

    printf ("Start docking\n");
    Run (lig, prt, psp, kde, mcs, enepara, temp, replica, &mcpara, mclog,
         multi_reps_records, complexsize);

    // post_mc(multi_reps_records, lig, prt, enepara, mcpara);
#ifdef IS_OPT == 1
    opt_ff(multi_reps_records, lig, complexsize.n_lig, prt, enepara, &mcpara);
#endif // IS_OPT == 1
    // printStates(multi_reps_records[0], inputfiles.trace_file.path);

    PrintSummary (&inputfiles, &mcpara, temp, mclog, &complexsize);

    // clean up

    delete[]mclog;
    delete[]lig;
    delete[]prt;
    delete[]psp;
    delete[]kde;
    delete[]mcs;
    delete[]enepara;
    delete[]temp;
    delete[]replica;

    return 0;
  }
  catch (std::exception & e) {
    std::cerr << "Unhandled Exception reached the top of main: " << e.what()
              << "GeauxDock will now exit" << std::endl;
    return ERROR_UNHANDLED_EXCEPTION;
  }

}
