#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <stdint.h>
#include "struc.h"
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>

}

using namespace std;

double myrand();
///
double myrandRange( double Min, double Max );
///
void smyrand( long seed );
///
double gengam(double aa,double bb);
///
double genbet(double aa,double bb);
///
double gennor(double av,double sd);
//
double gengau(double sd );
//
double gennorLtail( double av, double sd, double a );
//
double gennorRtail( double av, double sd, double a );
///
int genbinomial( int n, double p );
//
unsigned int genpoi( double );
//
uint32_t SampleFromDiscrete_new( std::vector<double> &cdf );
///
uint32_t SampleFromRange( uint32_t, uint32_t );
//
void genMultinomial(size_t K, unsigned int N, double *pbty, unsigned int *n );

double MyMulti(data_integer *Vector,int total);
double exponen(double x);
int SampleFromDiscrete_non_cum(vector<double> &pbty);
