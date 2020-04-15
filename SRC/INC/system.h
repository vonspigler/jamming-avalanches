#ifndef                     HEADERS_SYSTEM
#define                     HEADERS_SYSTEM

#include                    "library.h"
#include                    <cstdlib>
#include                    <string>
#include                    <memory>
#include                    <unordered_set> // c++11

extern std::string OUTPUT_FOLDER;

namespace Jamming {

namespace Metrics {

class metrics { // this includes the boundary conditions;
public:
  virtual ~metrics() = 0;
  virtual double dist_sq(unsigned int, unsigned int, const vectors<double> &, const vectors<double> &) const = 0;
  virtual double dist_1d(unsigned int, unsigned int, unsigned int, const vectors<double> &, const vectors<double> &) const = 0;
  virtual void regularize(vectors<double> &, int) const = 0;    // this function checks whether some particles are outside the unit box and puts them back
  std::string __NAME;                                           // second parameter: regularize only along dimension <d>
  double GAMMA;  // total shear-strain
};

}  // namespace Metrics

namespace Neighborhoods {

using uset = std::unordered_set<unsigned int>;
using uptr_uset = std::unique_ptr<uset []>;
using uptr_uint = std::unique_ptr<unsigned int []>;

class neighborhood {
public:
  virtual ~neighborhood() = 0;
  virtual void neighbors_of(unsigned int, uset &) const = 0;
  // reload is a function that recomputes the neighborhoods;
  // update can be called at every minimization step (eg md_fire), and decides when/what to reload, either by calling reload() or doing something else
  // update_on_shear is called when the system is sheared
  virtual void update(const vectors<double> &, const Metrics::metrics &) = 0;  // keep internal counter, in case
  virtual void reload(const vectors<double> &, const Metrics::metrics &) = 0;  // if needed, change threshold first; it should also call, or do also the same as, update()
  virtual void update_on_shear(const vectors<double> &, const Metrics::metrics &) = 0;
  std::string __NAME;
  unsigned int N;
  unsigned int dim;
};

}  // namespace Neighborhoods

namespace Interactions {

class interaction {
public:
  virtual ~interaction() = 0;
  virtual double H(const vectors<double> &, const Metrics::metrics&, const Neighborhoods::neighborhood &) const = 0;
  virtual void F(const vectors<double> &, vectors<double> &, const Metrics::metrics &, const Neighborhoods::neighborhood &) const = 0;  // XS, AS, METRICS, NEIGHBORHOOD
  virtual double stress_component(unsigned int, unsigned int, unsigned int, unsigned int, const vectors<double> &, const Metrics::metrics &, const Neighborhoods::neighborhood &) const = 0;
  std::string __NAME;
};

}  // namespace Interactions

namespace Minimizers {

class minimizer {
public:
  virtual ~minimizer() = 0;
  virtual void reset() = 0;
  virtual void minimize_step(vectors<double> &, vectors<double> &, vectors<double> &, const Metrics::metrics &, Neighborhoods::neighborhood &, const Interactions::interaction &) = 0; // XS VS AS; neighborhood is not CONST
  virtual bool check_convergence() const = 0; // it returns 1 if converged
  std::string __NAME;
};

}  // namespace Minimizers

const unsigned int VERBOSE_NONE = 0;
const unsigned int VERBOSE_MINIMAL = 1;
const unsigned int VERBOSE_ALL = 10;

class Spheres {
public:
  Spheres();
  void load();
  void set_radius();  // according to current phi
  void reload();
  double H();
  bool minimize();
  void perturb(double);
  void observables(double &, double &, double &); // energy, stress, z, for now
  double pressure();                              // pressure; temporary
  void print_conf(std::string) const;

  unsigned int __VERBOSE;  // degree of verbosity in terminal; 0 - minimal; 1 - function calls; 10 - object names, parameters and time
  unsigned int ITER_MAX;   // maximum number of cycles per particle for minimize()
  double R;
  unsigned int N;
  unsigned int dim;
  double phi;
  vectors<double> XS;
  vectors<double> VS;
  vectors<double> AS;

  // I use unique_ptr so that the objects are automatically destroyed when I change the specificities of the system
  std::unique_ptr<Metrics::metrics> METRICS;
  std::unique_ptr<Neighborhoods::neighborhood> NLISTS;
  std::unique_ptr<Interactions::interaction> INTERACTION;
  std::unique_ptr<Minimizers::minimizer> MINIMIZER;

private:
  double __dim_param;
};

}  // namespace Jamming

#endif
