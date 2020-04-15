#include                    <iostream>
#include                    <chrono>
#include                    "simulation_tools.h"
#include                    "system_implementation.h"
#include                    "system.h"

void shear(Jamming::Spheres &system, double g) { // g is a DELTA_GAMMA
  if (system.__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mShearing the system\033[0m (::shear, " << g << ")" << std::endl;

  //reload metrics -> g:
  system.METRICS->GAMMA += g;

  for (unsigned int i = 0; i < system.N; ++i) {
//    for (unsigned int d = 0; d < system.dim; ++d) <COMPUTE X_PREV>  // <-- needed to remove the affine transformation, or to study overlaps
    system.XS(i,0) += g*system.XS(i,1);  // AFFINE transformation
    system.METRICS->regularize(system.XS, 0);
  }

  //reload nlists? neighborlists: depends ong g..., noshearsubcells WRONG, leesedwardssubcells yes -> new neighborhood::update_on_shear(g) ???
  system.NLISTS->update_on_shear(system.XS, *system.METRICS);

  return;
}

void jam(Jamming::Spheres &system) {
double PREC_E = 1e-6;
double PREC_PHI = 1e-4;

  if (system.__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mFinding the jammed configuration\033[0m (jam)" << std::endl;
  double PHI_0 = 0, PHI_1 = 1, e = 0;

  unsigned int verb_bak = system.__VERBOSE;
  if (system.__VERBOSE == Jamming::VERBOSE_MINIMAL) system.__VERBOSE = Jamming::VERBOSE_NONE;

  std::chrono::time_point<std::chrono::system_clock> __start_time, __stop_time;

  if (system.__VERBOSE >= Jamming::VERBOSE_ALL) __start_time = std::chrono::system_clock::now();

  do {
    system.phi = (PHI_0 + PHI_1)/2;
    system.set_radius();
    system.NLISTS.release();
    system.NLISTS = std::unique_ptr<Jamming::Neighborhoods::neighborhood>(new Jamming::Neighborhoods::noshearsubcells(system.XS, (unsigned int)(std::max(ceil(1./2/system.R - 1), 1.0)), *system.METRICS));

    if (system.__VERBOSE >= Jamming::VERBOSE_ALL)
      std::cout << " > \033[1;37mNeighborhood: \033[0m" << system.NLISTS->__NAME << std::endl;
    system.INTERACTION = std::unique_ptr<Jamming::Interactions::interaction>(new Jamming::Interactions::harmonic(1, system.R));
    system.reload();
    system.minimize();

    e = system.H();
    if (system.__VERBOSE >= Jamming::VERBOSE_ALL)
      std::cout << "\t\033[1;36mphi=\033[0m" << system.phi << ", \033[1;36mR=\033[0m" << system.R << ": \033[1;36mE=\033[0m" << e << std::endl;
    if (e > PREC_E) PHI_1 = system.phi;
    else PHI_0 = system.phi;
  } while( (fabs(PHI_1 - PHI_0) > PREC_PHI) || (e > PREC_E) );

  if (system.__VERBOSE >= Jamming::VERBOSE_ALL) {
    __stop_time = std::chrono::system_clock::now();
    std::cout << "\t\033[1;31mTime to find jamming:\033[0m " << ((std::chrono::duration<double>)(__stop_time - __start_time)).count() << "s\n";
  }

//  system.phi += PREC_PHI;
  system.__VERBOSE = verb_bak;
  return;
}
