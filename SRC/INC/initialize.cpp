#include                    "simulation_tools.h"
#include                    "initialize.h"
#include                    "system_implementation.h"
#include                    "system.h"
#include                    "library.h"
#include                    <cstdlib>
#include                    <cstdio>
#include                    <cstring>
#include                    <cmath>
#include                    <iostream>
#include                    <chrono> // c++11
#include                    <string>
#include                    <algorithm>

double delta_gamma;
double max_gamma;
unsigned int SAMPLES;

bool prepare(Jamming::Spheres &system, int argc, char **argv) {
  if (argc != 10) {
    std::cout << " $ " << TERMCOL2 << "./MAIN.bin N dim phi dg gmax OUT/ seed VERB SAMPLES" << TERMEND << std::endl << std::endl;
    std::cout << TERMCOL3 << " N" << TERMEND << " \t Number of particles, e.g. 1000 (int>0)" << std::endl;
    std::cout << TERMCOL3 << " dim" << TERMEND << " \t Space dimension, e.g. 3 (int>0)" << std::endl;
    std::cout << TERMCOL3 << " phi" << TERMEND << " \t Packing fraction (double):" << std::endl;
    std::cout << " \t\t -1 to find jamming" << std::endl;
    std::cout << " \t\t >0 otherwise" << std::endl;

    std::cout << TERMCOL3 << " dg" << TERMEND << " \t Step in shear strain, e.g. 0.0001" << std::endl;
    std::cout << TERMCOL3 << " gmax" << TERMEND << " \t Maximum accumulated shear strain, e.g. 0.1" << std::endl;

    std::cout << TERMCOL3 << " OUT/" << TERMEND << " \t Folder where to output data; it must end with a slash (dir/)" << std::endl;
    std::cout << TERMCOL3 << " seed" << TERMEND << " \t Random number generator's seed (int)" << std::endl;
    std::cout << TERMCOL3 << " VERB" << TERMEND << " \t Verbose mode (bool)" << std::endl;
    std::cout << TERMCOL3 << " SAMPLES" << TERMEND << " \t Number of runs; outputs in OUT/<n>, n=0..SAMPLES-1" << std::endl;
    return false;
  }

  if (system.__VERBOSE) std::cout << TERMCOL2 << "Setting parameters" << TERMEND << std::endl;

  system.N = atoi(argv[1]);
  if (system.N < 1) {
    std::cout << TERMCOL1 << "N >= 1!" << TERMEND << std::endl;
    return false;
  }

  system.dim = atoi(argv[2]);
  if (system.dim <= 1) {
    std::cout << TERMCOL1 << "DIM > 1!" << TERMEND << std::endl;
    return false;
  }

  system.phi = atof(argv[3]);
  delta_gamma = atof(argv[4]);
  max_gamma = atof(argv[5]);

  OUTPUT_FOLDER = std::string(argv[6]);

  SEED = atoi(argv[7]);

  system.__VERBOSE = atoi(argv[8]);
  SAMPLES = atoi(argv[9]);

  if (SEED > 0) {
    set_seed(std::to_string(SEED));
  } else {
    SEED = std::chrono::system_clock::now().time_since_epoch().count();
    set_seed(std::to_string(SEED));
  }

  return true;
}

bool run(Jamming::Spheres &system, double thresh_ratio, double MRC) {
/*phi < phi_j   -->  noshearsubcells if not shearing, Nc >= Nc_min
                -->  with shear, I should implement shearsubcells (possibly buggy)
  phi >= phi_j  -->  neighborlists, small thresh, large MRC; I can shear the system*/

  system.load();
  system.METRICS = std::unique_ptr<Jamming::Metrics::metrics>(new Jamming::Metrics::leesedwards(0)); // GAMMA = 0 at the beginning
  system.MINIMIZER = std::unique_ptr<Jamming::Minimizers::minimizer>(new Jamming::Minimizers::md_fire()); // parameters?

  if (system.phi < 0) jam(system);

  system.set_radius();

  system.NLISTS = std::unique_ptr<Jamming::Neighborhoods::neighborhood>(new Jamming::Neighborhoods::neighborlists(system.XS, thresh_ratio*system.R, MRC, *system.METRICS));
  system.INTERACTION = std::unique_ptr<Jamming::Interactions::interaction>(new Jamming::Interactions::harmonic(1, system.R)); // coupling, cutoff

  if (system.__VERBOSE >= Jamming::VERBOSE_ALL) {
    std::cout << " > " << TERMCOL3 << "Metrics: " << TERMEND << system.METRICS->__NAME << std::endl; 
    std::cout << " > " << TERMCOL3 << "Neighborhood: " << TERMEND << system.NLISTS->__NAME << std::endl;
    std::cout << " > " << TERMCOL3 << "Interaction: " << TERMEND << system.INTERACTION->__NAME << std::endl;
    std::cout << " > " << TERMCOL3 << "Minimizer: " << TERMEND << system.MINIMIZER->__NAME << std::endl;
  }

  system.reload();
  system.minimize();

  if (system.__VERBOSE) std::cout << TERMCOL2 << "The system is ready" << TERMEND << std::endl;

  return true;
}

void print_info(Jamming::Spheres &system) {
  std::cout << TERMCOL1 << "\n+-------------------------------------------------------------------------------+\n" << TERMEND;
  std::cout << TERMCOL1 << "|  > SYSTEM INFO\033[0m (print_info)\t\t\t\t\t\t\t" << TERMCOL1 << "|\n" << TERMEND;
  std::cout << TERMCOL1 << "+" << TERMCOL2 << "-------------------------------------------------------------------------------" << TERMCOL1 << "+\n"  << TERMEND;

  std::cout << TERMCOL1 << "| N" << TERMEND << "\t\t\t= " << system.N << TERMCOL1 << "\t\t\t\t\t\t\t|\n" << TERMEND;
  std::cout << TERMCOL1 << "| dim" << TERMEND << "\t\t\t= " << system.dim << TERMCOL1 << "\t\t\t\t\t\t\t|\n" << TERMEND;
  std::cout << TERMCOL1 << "| phi" << TERMEND << "\t\t\t= " << system.phi << TERMCOL1 << "\t\t\t\t\t\t|\n" << TERMEND;
  std::cout << TERMCOL1 << "| R" << TERMEND << "\t\t\t= " << system.R << TERMCOL1 << "\t\t\t\t\t\t|\n" << TERMEND;

  std::cout << TERMCOL1 << "|" << TERMEND << "  > " << TERMCOL3 << "Metrics: \t\t  " << TERMEND << system.METRICS->__NAME << TERMCOL1 << "\t\t\t\t\t\t|\n" << TERMEND;
  std::cout << TERMCOL1 << "|" << TERMEND << "  > " << TERMCOL3 << "Neighborhood: \t  " << TERMEND << system.NLISTS->__NAME << TERMCOL1 << "\t\t\t\t\t\t|\n" << TERMEND;
  std::cout << TERMCOL1 << "|" << TERMEND << "  > " << TERMCOL3 << "Interaction: \t  " << TERMEND << system.INTERACTION->__NAME << TERMCOL1 << "\t\t\t\t\t\t|\n" << TERMEND;
  std::cout << TERMCOL1 << "|" << TERMEND << "  > " << TERMCOL3 << "Minimizer: \t  " << TERMEND << system.MINIMIZER->__NAME << TERMCOL1 << "\t\t\t\t\t\t|\n" << TERMEND;

  if (system.INTERACTION->__NAME == "harmonic")
    std::cout << TERMCOL1 << "| R (interaction)" << TERMEND << "\t= " << dynamic_cast<Jamming::Interactions::harmonic *>(system.INTERACTION.get())->Rc << TERMCOL1 << "\t\t\t\t\t\t|\n" << TERMEND;
  if (system.NLISTS->__NAME == "neighborlists")
    std::cout << TERMCOL1 << "| R (neighborlists)" << TERMEND << "\t= " << dynamic_cast<Jamming::Neighborhoods::neighborlists *>(system.NLISTS.get())->Rth << TERMCOL1 << "\t\t\t\t\t\t|\n" << TERMEND;
  else if (system.NLISTS->__NAME == "noshearsubcells")
    std::cout << TERMCOL1 << "| 1/Nc (noshearsubcells)" << TERMEND << "\t=" << 1/dynamic_cast<Jamming::Neighborhoods::noshearsubcells *>(system.NLISTS.get())->Nc << TERMCOL1 << "\t\t\t\t\t\t\t|\n" << TERMEND;

  std::cout << TERMCOL1 << "+-------------------------------------------------------------------------------+\n\n" << TERMEND;

  return;
}
