#include                    <cmath>
#include                    <cstdlib>
#include                    <iostream>
#include                    <fstream>
#include                    <sstream>
#include                    <string>
#include                    "INC/library.h"
#include                    "INC/simulation_tools.h"
#include                    "INC/initialize.h"
#include                    "INC/system.h"
#include                    "INC/system_implementation.h"

using namespace             std;


/*####  PROGRAM NOTES  #############################################################################
##                                                                                                ##
##  TODO:                                                                                         ##
## -- I AM ACTUALLY ASSUMING THAT |n| <= 1                                                        ##
## -- I want to improve the output code:                                                          ##
##    color codes? extension of cout?                                                             ##
## -- REWRITE observables()                                                                       ##
##                                                                                                ##
######  PHYSICS NOTES  #############################################################################
##                                                                                                ##
##  OUTPUT:                                                                                       ##
##      gamma (accumulated), energy (=>total), z, overlap (t-1,t: bc non additive), stress xy     ##
##          energy and extensive are intensive                                                    ##
##                                                                                                ##
##                                                                                                ##
##   I SHOULD COMPUTE:                                                                            ##
##     - eigenvector corresponding to smaller eigenvalue before avalanche                         ##
##       => participation ratio vs avalanche size (also in overlap)                               ##
##     - all this quantities all the way up to yielding; change in exponent; abrupt?              ##
##       => is there yielding? or shear-jamming and divergence of S?                              ##
##     - compare results for monodisperse LJ particles (I don't expect any scale free events)     ##
##     - overlap(0,t); as a function of delta_g or g_max: does it go to 0 or to a finite value?   ##
##                                                                                                ##
##   => do all this for N=64 128 256 512 1024 2048 4096 -> system size dependence (exponents)     ##
##   => do all this for dg=1e-3 1e-4 1e-5, perhaps some intermediate values -> scalings           ##
##                                                                                                ##
##################################################################################################*/

/*
  CONSTANT PRESSURE:
    the best way, imho, is to gradually increase the volume after each minimization following a shear step,
    until a target pressure is reached. Can I minimize in a bisection-like fashion? I am not so sure...

//curr_press = system.pressure();
//system.R *= exp(-(curr_press - 1e-5)/4/system.dim);
//system.NLISTS = std::unique_ptr<Jamming::Neighborhoods::neighborhood>(new Jamming::Neighborhoods::neighborlists(system.XS, 2.5*system.R, 100, *system.METRICS));
//system.INTERACTION = std::unique_ptr<Jamming::Interactions::interaction>(new Jamming::Interactions::harmonic(1, system.R)); // coupling, cutoff
//system.reload();

  NOT WORKING PROPERLY...
*/

extern double max_gamma;
extern double delta_gamma;

/*##########  MAIN  ##############################################################################*/

int main(int argc, char **argv) {
  std::cout << TERMCLEAR << TERMCOL1 << "  ==========  SIMULATION OF HARMONIC SPHERES UNDER SHEAR-STRAIN  ===============  " << TERMEND << std::endl;

  Jamming::Spheres system;
  if (!prepare(system, argc, argv)) return -1;

/*==========  START SAMPLES  =====================================================================*/

  for (unsigned int sample = 0; sample < SAMPLES; sample++) {
    std::cout << TERMCOL1 << " > SIMULATING SAMPLE #" << sample << " / " << SAMPLES << TERMEND << std::endl;
    run(system, 2.5, 10);  // set up and minimization; args; thresh_ratio -> 2.5*system.R, MCR (update count)

    auto n = OUTPUT_FOLDER + to_string(sample);
    ofstream OUT;
    OUT.open(n.c_str(), ios::out);
    print_info(system);  // this prints some information

/*----------  START  -----------------------------------------------------------------------------*/

    double curr_ene, curr_z, curr_s_xy;//, curr_press;  // note: energy and stress are extensive; number of contacts is per particle
    system.observables(curr_ene, curr_z, curr_s_xy);
    OUT << "# N = " << system.N << endl;
    OUT << "# dim = " << system.dim << endl;
    OUT << "# phi = " << system.phi << endl;
    OUT << "# delta_gamma = " << delta_gamma << endl;
    OUT << "# max_gamma = " << max_gamma << endl;
    OUT << "# R = " << system.R << endl;
    OUT << "# SEED = " << SEED << endl;
    OUT << "# GAMMA DGAMMA ENERGY STRESS DENERGY DSTRESS" << endl;
//    OUT << GAMMA << " " << delta_gamma << " " << curr_ene << " " << curr_s_xy << " " << DELTA_ene << " " << DELTA_s_xy << /*curr_z << " " << curr_press << " " <<*/ endl;

/*----------  SHEAR  -----------------------------------------------------------------------------*/

    for (double acc_gamma = delta_gamma; acc_gamma <= max_gamma; acc_gamma += delta_gamma) {
      double prev_ene = curr_ene, prev_s_xy = curr_s_xy;
      shear(system, delta_gamma);
      system.minimize();
      system.observables(curr_ene, curr_z, curr_s_xy);

      OUT << acc_gamma << " " << delta_gamma << " " << curr_ene << " " << curr_s_xy << " " << curr_ene - prev_ene << " " << curr_s_xy - prev_s_xy << endl;
    }

/*----------  END  -------------------------------------------------------------------------------*/

    OUT.close();

  }

/*==========  END SAMPLES  =======================================================================*/

  cout << TERMCOL1 << "\n  ==========  SIMULATION ENDS  =================================================  \n" << TERMEND << endl;

  return 0;
}
