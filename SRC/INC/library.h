#ifndef                     HEADERS_LIBRARY
#define                     HEADERS_LIBRARY

#include                    <memory> // c++11
#include                    <string>
#include                    <random>
#include                    <vector>

#define                     TERMCLEAR "\033[2J\033[H\n"
#define                     TERMEND "\033[0m"
#define                     TERMCOL1 "\033[1;31m"
#define                     TERMCOL2 "\033[1;36m"
#define                     TERMCOL3 "\033[1;37m"

extern int SEED;
extern unsigned int SAMPLES;

extern const double PREC;
extern const double PREC_PHI;
extern const double PREC_Z;
extern const unsigned int ITER_MAX;
extern double NLIST_RF_SQ;
extern unsigned int NLIST_EVERY;

extern std::default_random_engine __RAND_GENERATOR;
extern std::uniform_int_distribution<unsigned int> __RAND_INT;
extern std::uniform_real_distribution<double> __RAND_DOUBLE;
extern std::normal_distribution<double> __RAND_GAUSS;

inline unsigned int rand_int(unsigned int n)
  { return (unsigned int)(1.0*n*__RAND_INT(__RAND_GENERATOR)/__RAND_INT.max()); }
inline double rand_double(double n)
  { return n*__RAND_DOUBLE(__RAND_GENERATOR); }
inline double rand_gauss(double m, double s)
  { return __RAND_GAUSS(__RAND_GENERATOR)*s + m; }
inline double abs_v(double x) { return (x < 0 ? -x : x); }

std::string get_seed();
void set_seed(std::string);

inline bool draw_prob(double p) { return (1.0*rand()/RAND_MAX < p); }
inline double mod(double x, unsigned int n) { return x - n*floor(x/n); }

template <class T> class vectors {
private:
  T *data;

public:
  unsigned int NN, DD;
  vectors();
  vectors(unsigned int, unsigned int);
  vectors(const vectors &);
  ~vectors();

  vectors<T>& operator=(const vectors<T> &);

  inline T& operator()(unsigned int i, unsigned int d) {
//    if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
//    if (d >= DD) std::cout << "ERROR, d (" << d << ") must be < DIM! (" << DD << ")" << std::endl;
    return data[i*DD + d];
  };
  inline const T& operator()(unsigned int i, unsigned int d) const {
//    if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
//    if (d >= DD) std::cout << "ERROR, d (" << d << ") must be < DIM! (" << DD << ")" << std::endl;
    return data[i*DD + d];
  };
  inline const T* operator()(unsigned int i) const {
//    if (i >= NN) std::cout << "ERROR, i (" << i << ") must be < N! (" << NN << ")" << std::endl;
    return &data[i*DD];
  }; // READ-ONLY to avoid mess; call as const double* r = x(i);
};

extern template class vectors<double>;
extern template class vectors<unsigned int>;

#endif
