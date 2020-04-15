#include                    <iostream>
#include                    "system.h"
#include                    "system_implementation.h"

namespace Jamming {

/*##############################################################################  METRICS  ######*/

namespace Metrics {

leesedwards::leesedwards(double g) {
  __NAME = "leesedwards";
  GAMMA = g;
}

double leesedwards::dist_sq(unsigned int i, unsigned int j, const vectors<double> &XS, const vectors<double> &YS) const {
  double dist, r = 0, del;

  short n;
  del = YS(j,1) - XS(i,1);

  if (del > 0.5) n = 1;
  else if (del < -0.5) n = -1;
  else n = 0;

  del = YS(j,0) - XS(i,0) - GAMMA*n;

  r += del;
  if (del > 0.5) --r;
  else if (del < -0.5) ++r;

  dist = r*r;

  for (unsigned int d = 1; d < XS.DD; ++d) {  // only d > 0
    r = 0;

    del = YS(j,d) - XS(i,d);

    r += del;
    if (del > 0.5) --r;
    else if (del < -0.5) ++r;

    dist += r*r;
  }

  return dist;
}

double leesedwards::dist_1d(unsigned int j, unsigned int i, unsigned int d, const vectors<double> &XS, const vectors<double> &YS) const {
  double del;

  if (d != 0) {
    del = YS(j,d) - XS(i,d);
    if (del > 0.5) return del - 1;
    else if (del < -0.5) return del + 1;
    else return del;

  } else {
    short n;
    del = YS(j,1) - XS(i,1);

    // I AM ACTUALLY ASSUMING THAT |n| <= 1 --> ok since the two particles will be belonging to neighboring cells
    if (del > 0.5) n = 1;
    else if (del < -0.5) n = -1;
    else n = 0;

    del = YS(j,0) - XS(i,0) - GAMMA*n;

    if (del > 0.5) return del - 1;
    else if (del < -0.5) return del + 1;
    else return del;
  }
}

void leesedwards::regularize(vectors<double> &XS, int d) const {
  int n_disp;

  if (d < 0) {
    for (unsigned int i = 0; i < XS.NN; ++i) {
      n_disp = floor(XS(i,1));
      if (n_disp != 0) XS(i,0) -= GAMMA*n_disp; // no VS, AQS
      for (unsigned int dd = 0; dd < XS.DD; ++dd) if ((XS(i,dd) >= 1) || (XS(i,dd) < 0)) XS(i,dd) -= floor(XS(i,dd));
    }

  } else {
    // assuming d < dim
    for (unsigned int i = 0; i < XS.NN; ++i) {
      n_disp = floor(XS(i,1));
      if (n_disp != 0) XS(i,0) -= GAMMA*n_disp; // no VS, AQS
      if ((XS(i,d) >= 1) || (XS(i,d) < 0)) XS(i,d) -= floor(XS(i,d));
    }
  }

  return;
}

}  // namespace Metrics

/*########################################################################  NEIGHBORHOODS  ######*/

namespace Neighborhoods {

neighborlists::neighborlists(const vectors<double> &XS, double th, double MRC, const Metrics::metrics &DIST) {
  __NAME = "neighborlists";
  N = XS.NN;
  dim = XS.DD;
  Rth = th;
  MAX_RELOAD_COUNTER = MRC;
  reload_counter = 0;
  _neighbors_of = uptr_uset(new uset [N]);

  reload(XS, DIST);
}

void neighborlists::reload(const vectors<double> &XS, const Metrics::metrics &DIST) {
  for (unsigned int i = 0; i < N; ++i) _neighbors_of[i].clear();

  double th_sq = Rth*Rth;

  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      if (j <= i) continue; // it's symmetric: i in N(j) iff j in N(i)

      if (DIST.dist_sq(i, j, XS, XS) < th_sq) {
        _neighbors_of[i].insert(j);
        _neighbors_of[j].insert(i);
      }
    }
  }

  reload_counter = 0;

  return;
}

void neighborlists::update(const vectors<double> &XS, const Metrics::metrics &DIST) {  // it must keep an internal counter
  if (++reload_counter >= MAX_RELOAD_COUNTER) reload(XS, DIST);
  return;
}

void neighborlists::update_on_shear(const vectors<double> &XS, const Metrics::metrics &DIST) {
// in principle, for small g I can not reload...
  reload(XS, DIST);
  return;
}

void neighborlists::neighbors_of(unsigned int i, uset &neighbors_set) const {
  neighbors_set.clear();
  for (auto j : _neighbors_of[i]) neighbors_set.insert(j);
  return;
}

/*-----------------------------------------------------------------------------------------------*/

noshearsubcells::noshearsubcells(const vectors<double> &XS, double NSC, const Metrics::metrics &DIST) {  // XS, num of subcells, metrics
  __NAME = "noshearsubcells";
  N = XS.NN;
  dim = XS.DD;
  Nc = NSC;

  POW_3 = uptr_uint(new unsigned int [dim + 2]);
  POW_Nc = uptr_uint(new unsigned int [dim + 2]);
  POW_3[0] = 1;
  POW_Nc[0] = 1;

  for (unsigned int d = 1; d < dim + 2; d++) {
    POW_3[d] = 3*POW_3[d - 1];
    POW_Nc[d] = Nc*POW_Nc[d - 1];
  }

  _neighbors_of = uptr_uset(new uset [POW_Nc[dim]]);  // IMPORTANT: TO CHECK : I think that _neighbors_of contains the cell itself!
  _cells_of_particles = uptr_uint(new unsigned int [N]);
  _particles_in_cell = uptr_uset(new uset [POW_Nc[dim]]);
  _num_of_particles = uptr_uint(new unsigned int [POW_Nc[dim]]);

  __loaded_once = false;

  for (unsigned int i = 0; i < N; ++i) _cells_of_particles[i] = 0;

  // this loop defines the hypercubic topology, ie _neighbors[]
  for (unsigned int c = 0; c < POW_Nc[dim]; ++c) {
    _num_of_particles[c] = 0;

    for (unsigned int n = 0; n < POW_3[dim]; ++n) {
      unsigned int neighbor = 0;

      for (unsigned int d = 0; d < dim; ++d) {
        unsigned int cell_digit = (unsigned int) mod(floor(1.0*c/POW_Nc[d]), Nc);
        int delta_digit = (int) (mod(floor(1.0*n/POW_3[d]), 3) - 1);
        int neighbor_digit = mod(int(cell_digit) + delta_digit, Nc);
        neighbor += neighbor_digit*POW_Nc[d];
      }

      _neighbors_of[c].insert(neighbor);
    }
  }

  reload(XS, DIST);
}

void noshearsubcells::reload(const vectors<double> &XS, const Metrics::metrics &DIST) {
  if (__loaded_once == false) {
    for (unsigned int c = 0; c < POW_Nc[dim]; ++c) _num_of_particles[c] = 0;

    for (unsigned int i = 0; i < N; ++i) {
      unsigned int c_new = 0;
      for (unsigned int d = 1; d < dim; d++) c_new += ((unsigned int) mod(floor(XS(i,d)*Nc), Nc))*POW_Nc[d];
      _cells_of_particles[i] = c_new;
      ++_num_of_particles[c_new];
      _particles_in_cell[c_new].insert(i);
    }

    __loaded_once = true;

  } else {

    for (unsigned int i = 0; i < N; ++i) {
      unsigned int c = cell_of(i), c_new = 0;
      for (unsigned int d = 1; d < dim; d++) c_new += ((unsigned int) mod(floor(XS(i,d)*Nc), Nc))*POW_Nc[d];

      if (c != c_new) move(i, c, c_new);
    }
  }

  return;
}

void noshearsubcells::update(const vectors<double> &XS, const Metrics::metrics &DIST) {  // keep internal counter
  unsigned c_new, c;

  for (unsigned int i = 0; i < XS.NN; ++i) {
    c_new = 0;
    for (unsigned int d = 0; d < XS.DD; ++d) c_new += ((unsigned int) mod(floor(XS(i,d)*Nc), Nc))*POW_Nc[d];

    c = _cells_of_particles[i];
    if (c != c_new) move(i, c, c_new);
  }

  return;
}

void noshearsubcells::update_on_shear(const vectors<double> &XS, const Metrics::metrics &DIST) {
// in principle, for small g I can not reload...
  reload(XS, DIST);
  return;
}

void noshearsubcells::move(unsigned int i, unsigned int c, unsigned int c_new) {
  _particles_in_cell[c].erase(i);
  --_num_of_particles[c];

  _cells_of_particles[i] = c_new;
  ++_num_of_particles[c_new];
  _particles_in_cell[c_new].insert(i);

  return;
}

void noshearsubcells::neighbors_of(unsigned int i, uset &neighbors_set) const {
  neighbors_set.clear();
  unsigned int c = _cells_of_particles[i];

// I have to erase(i), because _neighbors_of[c(i)] includes  c(i) itself
  for (const auto n : _neighbors_of[c])
    for (const auto j : _particles_in_cell[n])
      neighbors_set.insert(j);
  neighbors_set.erase(i);

  return;
}

inline unsigned int noshearsubcells::cell_of(unsigned int i) const {
  return _cells_of_particles[i];
}

}  // namespace Neighborhoods

/*#########################################################################  INTERACTIONS  ######*/

namespace Interactions {

harmonic::harmonic(double e, double c) {
  __NAME = "harmonic";
  eps = e;
  Rc = c;
}

double harmonic::H(const vectors<double> &XS, const Metrics::metrics &DIST, const Neighborhoods::neighborhood &NLISTS) const {
  double r = 0, __dist;

  for (unsigned int i = 0; i < NLISTS.N; ++i) {
    Jamming::Neighborhoods::uset neighbors_set;
    NLISTS.neighbors_of(i, neighbors_set);

    for (const auto j : neighbors_set) {
      double __dist_sq = DIST.dist_sq(i, j, XS, XS); // change name `DD'....

      if (__dist_sq < 4.0*Rc*Rc) {
        __dist = sqrt(__dist_sq);

        r += eps*(1 - __dist/2/Rc)*(1 - __dist/2/Rc);
      }
    }
  }

  return r/2;
}

void harmonic::F(const vectors<double> &XS, vectors<double> &AS, const Metrics::metrics &DIST, const Neighborhoods::neighborhood &NLISTS) const {
  double r, __dist;

  for (unsigned int i = 0; i < XS.NN; ++i) {
    Jamming::Neighborhoods::uset neighbors_set;
    NLISTS.neighbors_of(i, neighbors_set);

    for (unsigned int d = 0; d < XS.DD; ++d) {
      r = 0;

      for (const auto j : neighbors_set) {
        __dist = sqrt(DIST.dist_sq(i, j, XS, XS));
        if (__dist < 2*Rc) r += eps*(1 - __dist/2/Rc)*DIST.dist_1d(i, j, d, XS, XS)/Rc/__dist;
      }
      AS(i,d) = r;
    }
  }

  return;
}

// i, j, m, n... -> r_ij^d F_ij^d
double harmonic::stress_component(unsigned int i, unsigned int j, unsigned int m, unsigned int n, const vectors<double> &XS, const Metrics::metrics &DIST, const Neighborhoods::neighborhood &NLISTS) const {
  double dist_sq = DIST.dist_sq(i, j, XS, XS);
  if (dist_sq < 4*Rc*Rc) {
    double dist = sqrt(dist_sq);
    return eps*(1 - dist/2/Rc)*DIST.dist_1d(i, j, m, XS, XS)*DIST.dist_1d(i, j, n, XS, XS)/Rc/dist;
  } else {
    return 0;
  }
}

}  // namespace Interactions

/*###########################################################################  MINIMIZERS  ######*/

namespace Minimizers {

md_fire::md_fire() {
  __NAME = "md_fire";

  MD_PREC = 1e-4; // IMPROVE
  MD_alpha_start = 0.1;
  MD_delta_t_max = 0.01; // check!
  MD_f_alpha = 0.99;
  MD_f_inc = 1.1;
  MD_f_dec = 0.5;
  MD_N_min = 5;
}

void md_fire::reset() {
  MD_alpha = MD_alpha_start;
  MD_delta_t = MD_delta_t_max/10;
  MD_counter = 0;
  return;
}

bool md_fire::check_convergence() const {
  return bool(MD_A_norm < MD_PREC); // 1 if converged
}

void md_fire::minimize_step(vectors<double> &XS, vectors<double> &VS, vectors<double> &AS, const Metrics::metrics &DIST, Neighborhoods::neighborhood &NLISTS, const Interactions::interaction &INTERACTION) {
  INTERACTION.F(XS, AS, DIST, NLISTS);

  // Velocity Verlet, unit mass
  for (unsigned int i = 0; i < XS.NN; ++i) {
    for (unsigned int d = 0; d < XS.DD; ++d) {
      XS(i,d) += VS(i,d)*MD_delta_t + AS(i,d)*MD_delta_t*MD_delta_t/2;
      VS(i,d) += AS(i,d)*MD_delta_t/2; // updating V's and X's
    }
  }

  // -1 option: regularizes all dimensions
  DIST.regularize(XS, -1);
  // this function updates the neighborhoods:
  //   - in the case of neighborlists, the class will keep track of an internal counter, reset also via reset()
  //   - for noshearsubcells, it will be update at every step
  NLISTS.update(XS, DIST);

  INTERACTION.F(XS, AS, DIST, NLISTS);

  MD_A_norm = 0;
  MD_V_norm = 0;

  for (unsigned int i = 0; i < XS.NN; ++i)
    for (unsigned int d = 0; d < XS.DD; ++d) {
      VS(i,d) += AS(i,d)*MD_delta_t/2;
      MD_V_norm += VS(i,d)*VS(i,d);
      MD_A_norm += AS(i,d)*AS(i,d);
  }

  MD_V_norm = sqrt(MD_V_norm);
  MD_A_norm = sqrt(MD_A_norm);

  double A_norm, P = 0;

  // FIRE
  for (unsigned int i = 0; i < XS.NN; ++i)
    for (unsigned int d = 0; d < XS.DD; ++d) {
      A_norm = (MD_A_norm >= 1e-8 ? AS(i,d)/MD_A_norm : 0);  // I need to normalize A, but if the norm is already small I introduce numerical errors
      VS(i,d) = (1 - MD_alpha)*VS(i,d) + MD_alpha*MD_V_norm*A_norm;
      P += VS(i,d)*AS(i,d);
  }

  if (P <= 0) {
    for (unsigned int i = 0; i < XS.NN; ++i)
      for (unsigned int d = 0; d < XS.DD; ++d)
        VS(i,d) = 0;

    MD_alpha = MD_alpha_start;
    MD_counter = 0;
    MD_delta_t *= MD_f_dec;

  } else {
    if (MD_counter > MD_N_min) {
      MD_delta_t *= MD_f_inc;
      if (MD_delta_t > MD_delta_t_max) MD_delta_t = MD_delta_t_max;

      MD_alpha *= MD_f_alpha;
    }

    MD_counter++;
  }

  return;
}

}  // namespace Minimizers

}  // namespace Jamming
