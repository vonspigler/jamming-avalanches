#ifndef                     HEADERS_INITIALIZE
#define                     HEADERS_INITIALIZE

#include                    "system.h"
#include                    "system_implementation.h"
#include                    "library.h"
#include                    <cstdlib>
#include                    <cmath>
#include                    <iostream>
#include                    <chrono> // c++11
#include                    <string>

extern double max_gamma;
extern double delta_gamma;

bool prepare(Jamming::Spheres &, int, char**);
bool run(Jamming::Spheres &, double, double);
void print_info(Jamming::Spheres &);

#endif
