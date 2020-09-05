#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>

#ifndef STATISTICS_H
#define STATISTICS_H


double logPoisson(int x, double lambda);

double ilogCumulativePoisson(int x, double lambda);


float gammln(float xx);
double gammlnD(double xx);

double betacfD(double a, double b, double x);

float factln(int n);

void gcf(float a, float x, float &gammcf, float& gln);

void gser(float a, float x, float &gamser, float &gln);

#define HYPERGEO_CACHE_SIZE 1000000

#endif
