#include <iostream>
#include <vector>
#include <math.h>
#include <set>
#include <map>
#include <cassert>
#include <random>

/*
  ex = exp(-al*(b-a0))
  dex = -al*ex, ddex = al*al*ex

  E = D*( ex*ex -2.0*ex + 1.0)
  F = 2.0*D*al*ex*(ex-1.0)
  C = = -dF = 2.0*D*al*al*(2.0*ex-1)*ex
  C = 0 => ex = 0.5
  ex=1 => C = 2*D*al*al
*/
#define LOGTWO 0.6931471805599453
#define PI 3.141592654
#define RTTWO 1.4142135623730951

double denergy(double d, double al, double a0, double b) {
  /*
  d,al,a0 = Morse parameters
  b = bond magnitude
  */
  return 2.0*d*al*exp(-al*(b-a0))*(1.0-exp(-al*(b-a0))); // -> 2.0*D*AL*AL*(b-a0)
};

double energy(double d, double al, double a0, double b) {
  return d*(1.0-exp(-al*(b-a0)))*(1.0-exp(-al*(b-a0))); // -> D*AL*AL*(b-a0)^2
};


/*
d,al,a0 = Morse parameters
a = lattice vector for bond
btm = transverse magnitude
blm = longitudinal magnitude
*/

double energy(double d, double al, double rt, double a0, double blm, double btm) {
  // longitudinal energy
  double el = d * ( 1.0 - exp(-al*(blm-a0)) ) * ( 1.0 - exp(-al*(blm-a0)) );
  // transverse energy
  double et = rt*d*al*al*btm*btm;
  return el + et;
};

double energy_long(double d, double al, double rt, double a0, double blm) {
  // longitudinal energy
  return d * ( 1.0 - exp(-al*(blm-a0)) ) * ( 1.0 - exp(-al*(blm-a0)) );
};

double energy_tran(double d, double al, double rt, double a0, double btm) {
  // transverse energy
  return rt*d*al*al*btm*btm;
};


double force(double d, double al, double rt, double a0, double blm, double btm, double *f) {
  // force longitudinal
  f[0] = 2.0 * d * al * exp(-al*(blm-a0)) * (exp(-al*(blm-a0))-1.0); // -> -2.0*d*al*al*db
  // force traditional
  f[1] = -2.0 * rt * d * al * al;

  return sqrt(f[0]*f[0]+f[1]*f[1]*btm*btm);
};
