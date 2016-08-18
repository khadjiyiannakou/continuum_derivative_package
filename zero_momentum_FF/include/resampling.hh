// AUTHOR: Konstantin Ottnad
// 
// DATE: 20140707
//
// Implements basic functions for resampling (bootstrap, jackknife).
//
// TODO:

#ifndef resampling_hh
#define resampling_hh

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <limits>
#include <cmath>
#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>

using namespace std;

int bootstrap(gsl_vector *IN, gsl_vector *&OUT, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed);
int bootstrap(gsl_vector *IN, gsl_vector *&OUT, double &ERR, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed);
int bootstrap(gsl_vector_complex *IN, gsl_vector_complex *&OUT, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed);
int bootstrap(gsl_vector_complex *IN, gsl_vector_complex *&OUT, gsl_complex &ERR, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed);
int bootstrap(gsl_matrix *IN, gsl_matrix *&OUT, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed);
int bootstrap(gsl_matrix *IN, gsl_matrix *&OUT, gsl_vector *&ERR, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed);
int bootstrap(gsl_matrix_complex *IN, gsl_matrix_complex *&OUT, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed);
int bootstrap(gsl_matrix_complex *IN, gsl_matrix_complex *&OUT, gsl_vector_complex *&ERR, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed);
// general function for bootstrapping (including blocking). GSL matrix index 1 index for IN is expected to run over gauge index and OUT is allocated with the appropriate number of samples 'no_samples' in the first index. The second index contains all other internal indices.


int jackknife(gsl_vector *IN, gsl_vector *&OUT, const unsigned int blocksize);
int jackknife(gsl_vector *IN, gsl_vector *&OUT, double &ERR, const unsigned int blocksize);
int jackknife(gsl_vector_complex *IN, gsl_vector_complex *&OUT, const unsigned int blocksize);
int jackknife(gsl_vector_complex *IN, gsl_vector_complex *&OUT, gsl_complex &ERR, const unsigned int blocksize);
int jackknife(gsl_matrix *IN, gsl_matrix *&OUT, const unsigned int blocksize);
int jackknife(gsl_matrix *IN, gsl_matrix *&OUT, gsl_vector *&ERR, const unsigned int blocksize);
int jackknife(gsl_matrix_complex *IN, gsl_matrix_complex *&OUT, const unsigned int blocksize);
int jackknife(gsl_matrix_complex *IN, gsl_matrix_complex *&OUT, gsl_vector_complex *&ERR, const unsigned int blocksize);
// general function for jacknifing (including blocking). GSL matrix index 1 index for IN is expected to run over gauge index and OUT is allocated with the appropriate number of samples = Ncount/blocksize+1 in the first index. The second index contains all other internal indices.

void polar_method(double *z); // returns a pair of independent, normal distributed random numbers
void resampling(const double mean, const double err, const unsigned int no_samples, double *results); // performs resampling (normal distributed)
#endif
