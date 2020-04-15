#ifndef                     HEADERS_SYSTEM_IMPLEMENTATION
#define                     HEADERS_SYSTEM_IMPLEMENTATION

#include                    "system.h"
#include                    "library.h"

namespace Jamming {

/*##############################################################################  METRICS  ######*/

namespace Metrics {

class leesedwards : public metrics {
public:
  leesedwards(double);
  double dist_sq(unsigned int, unsigned int, const vectors<double> &, const vectors<double> &) const;
  double dist_1d(unsigned int, unsigned int, unsigned int, const vectors<double> &, const vectors<double> &) const;
  void regularize(vectors<double> &, int) const;
};

}  // namespace Metrics

/*########################################################################  NEIGHBORHOODS  ######*/

namespace Neighborhoods {

class neighborlists : public neighborhood {
public:
  neighborlists(const vectors<double> &, double, double, const Metrics::metrics &);  // XS, threshold, max reload counter, metrics
  void reload(const vectors<double> &, const Metrics::metrics &);
  void update(const vectors<double> &, const Metrics::metrics &);                    // keep internal counter
  void update_on_shear(const vectors<double> &, const Metrics::metrics &);
  void neighbors_of(unsigned int i, uset &) const;

  double Rth; // threshold
  unsigned int reload_counter;
  unsigned int MAX_RELOAD_COUNTER;

private:
  uptr_uset _neighbors_of;
};

/*-----------------------------------------------------------------------------------------------*/

class noshearsubcells : public neighborhood {
public:
  noshearsubcells(const vectors<double> &, double, const Metrics::metrics &);  // XS, num of subcells, metrics
  void reload(const vectors<double> &, const Metrics::metrics &);
  void update(const vectors<double> &, const Metrics::metrics &);
  void neighbors_of(unsigned int i, uset &) const;                              // it gives the list of particles neighboring with i (i not included)
  void update_on_shear(const vectors<double> &, const Metrics::metrics &);
  inline unsigned int cell_of(unsigned int i) const;                           // the cell of i
  void move(unsigned int, unsigned int, unsigned int);                         // particle, old cell, new cell

  unsigned int Nc;   // num of subcells
  uptr_uint POW_Nc;  // this array contains the powers Nc^d, for d in {1, ..., dim}, so that I do not have to call pow() every time...
  uptr_uint POW_3;   // """

private:
  uptr_uset _particles_in_cell;   // _particles_in_cell[c] is the list of particles in subcell c
  uptr_uint _cells_of_particles;  // this is the inverse, _cells_of_particles[i] is the subcell of particle i

  uptr_uset _neighbors_of;        // internal representation of the subcells' topology: _neighbors[c] is the list of neighboring subcells
  uptr_uint _num_of_particles;    // [c] -> number of particles in cell c
  bool __loaded_once;             // I use this to slightly optimize the reloading if the particles were already loaded into the subcells
};

}  // namespace Neighborhoods

/*#########################################################################  INTERACTIONS  ######*/

namespace Interactions {

class harmonic : public interaction {
public:
  harmonic(double, double);
  // XS, ...
  double H(const vectors<double> &, const Metrics::metrics &, const Neighborhoods::neighborhood &) const;                   // global energy and forces
  // XS, AS, ...
  void F(const vectors<double> &, vectors<double> &, const Metrics::metrics &, const Neighborhoods::neighborhood &) const;

  // stress term r_ij^d F_ij^d; i, j, d, XS, ...
  double stress_component(unsigned int, unsigned int, unsigned int, unsigned int, const vectors<double> &, const Metrics::metrics &, const Neighborhoods::neighborhood &) const;

  double eps;
  double Rc;
};

}  // namespace Interactions

/*###########################################################################  MINIMIZERS  ######*/

namespace Minimizers {

class md_fire : public minimizer {
public:
  md_fire();
  void reset();
  bool check_convergence() const;
  void minimize_step(vectors<double> &XS, vectors<double> &VS, vectors<double> &AS, const Metrics::metrics &, Neighborhoods::neighborhood &, const Interactions::interaction &);
  double MD_V_norm;
  double MD_A_norm;

private:
  double MD_PREC;
  double MD_alpha_start;
  double MD_delta_t_max;
  double MD_f_alpha;
  double MD_f_inc;
  double MD_f_dec;
  unsigned int MD_N_min;
  double MD_alpha;
  double MD_delta_t;
  unsigned int MD_counter;
};

}  // namespace Minimizers

}  // namespace Jamming

#endif
