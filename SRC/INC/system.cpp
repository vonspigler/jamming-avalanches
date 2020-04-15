#include                    "system.h"
#include                    <cstdlib>
#include                    <iostream>
#include                    <fstream>
#include                    <sstream>
#include                    <chrono> // c++11
#include                    "library.h"

extern std::string OUTPUT_FOLDER;

namespace Jamming {

namespace Metrics { metrics::~metrics() {} }
namespace Neighborhoods { neighborhood::~neighborhood() {} }
namespace Interactions { interaction::~interaction() {} }
namespace Minimizers { minimizer::~minimizer() {} }

Spheres::Spheres() {
  __VERBOSE = 0;
  ITER_MAX = 10000;
}

bool Spheres::minimize() {
  std::chrono::time_point<std::chrono::system_clock> __start_time, __stop_time;

  if (__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mMinimizing\033[0m (::minimize)" << std::endl;
  if (__VERBOSE >= Jamming::VERBOSE_ALL)
    __start_time = std::chrono::system_clock::now();

  MINIMIZER->reset();
  unsigned int counter = 0;

  do {
//    if (__VERBOSE) std::cout << "[#" << counter << "] " << H() << " " << XS(0,0) << std::endl;
    MINIMIZER->minimize_step(XS, VS, AS, *METRICS, *NLISTS, *INTERACTION);
    ++counter;
//    if (counter % NLIST_EVERY == 0) NLISTS.reload(NLIST_RF_SQ*R*R);  -----> in minimizer
  } while ((counter < ITER_MAX*N) && !MINIMIZER->check_convergence());

  if (__VERBOSE >= Jamming::VERBOSE_ALL) {
    __stop_time = std::chrono::system_clock::now();
    std::cout << "\t\033[1;31mTime:\033[0m " << ((std::chrono::duration<double>)(__stop_time - __start_time)).count() << "s\n";
  }

  if (counter == ITER_MAX*N) return false;
  return true;
}

void Spheres::load() {
  if (__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mLoading the particles\033[0m (::load)" << std::endl;

  XS = vectors<double>(N, dim);
  VS = vectors<double>(N, dim);
  AS = vectors<double>(N, dim);

  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int d = 0; d < dim; ++d) {
      XS(i,d) = rand_double(1);
      VS(i,d) = 0;
      AS(i,d) = 0;
    }
  }

  return;
}

void Spheres::set_radius() {
  if (__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mSetting R = R(phi, N, dim)\033[0m (::set_radius)" << std::endl;

  __dim_param = tgamma(dim/2.0 + 1)/pow(M_PI, dim/2.0);
  if (phi > 0) R = pow(phi/N*__dim_param, 1.0/dim);
  else std::cout << "\033[1;31mERROR: trying to set R = R(phi) with phi < 0!\033[0m" << std::endl;

  return;
}

void Spheres::reload() { // computes new R;
  if (__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mReloading the particles\033[0m (::reload)" << std::endl;
  NLISTS->reload(XS, *METRICS);
  return;
}

void Spheres::perturb(double r) {
  if (__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mPerturbing the system\033[0m (::perturb, " << r << ")" << std::endl;

  for (unsigned int i = 0; i < N; ++i)
    for (unsigned int d = 0; d < dim; ++d)
      XS(i,d) += R*rand_double(r);
}

double Spheres::H() {
  return INTERACTION->H(XS, *METRICS, *NLISTS);
}

void Spheres::observables(double &e, double &z, double &sxy) {
  if (__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mComputing observables\033[0m (::observables)" << std::endl;

  double t_z = 0, t_sxy = 0, dist_sq;

  for (unsigned int i = 0; i < N; ++i) {
    Jamming::Neighborhoods::uset neighbors_set;
    NLISTS->neighbors_of(i, neighbors_set);

    for (const auto j : neighbors_set) {
      dist_sq = METRICS->dist_sq(i, j, XS, XS);


      if (dist_sq < 4.0*R*R) ++t_z;
      t_sxy -= INTERACTION->stress_component(i, j, 0, 1, XS, *METRICS, *NLISTS);
    }
  }

  e = INTERACTION->H(XS, *METRICS, *NLISTS);  // extensive
  z = t_z/N;                                  // intensive
  sxy = t_sxy/2;                              // extensive

  return;
}

double Spheres::pressure() {
  if (__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mComputing pressure\033[0m (::observables)" << std::endl;

  double p = 0;

  for (unsigned int i = 0; i < N; ++i) {
    Jamming::Neighborhoods::uset neighbors_set;
    NLISTS->neighbors_of(i, neighbors_set);

    for (auto j : neighbors_set)
      for (unsigned int d = 0; d < dim; ++d) p -= INTERACTION->stress_component(i, j, d, d, XS, *METRICS, *NLISTS);
  }

  return -p/6;  // factor 2 bc sum_i<j = 1/2 sum_i!=j
}

void Spheres::print_conf(std::string name) const {
  if (__VERBOSE >= Jamming::VERBOSE_MINIMAL)
    std::cout << "\033[1;36mPrinting configuration\033[0m (::print_conf, `" << OUTPUT_FOLDER << "conf_" << name.c_str() << "')" << std::endl;

  std::ofstream OUT_conf;
  std::string n(OUTPUT_FOLDER);
  n += name;
  n += ".conf";
  OUT_conf.open(n.c_str(), std::ios::out);

  for (unsigned int i = 0; i < N; ++i) {
    OUT_conf << i << " " << R << " ";
    for (unsigned int d = 0; d < dim; ++d) OUT_conf << XS(i,d) << " ";
    OUT_conf << std::endl;
  }

  OUT_conf.close();

  return;
}

}  // namespace Jamming
