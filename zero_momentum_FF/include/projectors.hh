// AUTHOR: Konstantin Ottnad
//
// DATE: 2140217
//
// This files defines the projectors that are required in the analysis of nucleon n-pt functions
// 
// Using standard (Eucledian) chiral rep for Dirac matrices:
//
//        | 0         e_mu | 
// g_mu = |                |,  g_4 = 1,   g_5 = g_0 * g_1 * g_2 * g_3
//        | e_mu^dag  0    |
//
// where e_0 = -1,  e_k = -i tau_k, tau_k standard Pauli matrices

#ifndef projectors_hh
#define projectors_hh

#include <iomanip>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

using namespace std;

const unsigned int no_projectors = 31; // 25 + 6 gamma matrices

const double dirac_matrices[6][4][4][2] = 
 {{{{+0.,+0.},{+0.,+0.},{-1.,+0.},{+0.,+0.}},
   {{+0.,+0.},{+0.,+0.},{+0.,+0.},{-1.,+0.}},
   {{-1.,+0.},{+0.,+0.},{+0.,+0.},{+0.,+0.}},
   {{+0.,+0.},{-1.,+0.},{+0.,+0.},{+0.,+0.}}},

  {{{+0.,+0.},{+0.,+0.},{+0.,+0.},{+0.,-1.}},
   {{+0.,+0.},{+0.,+0.},{+0.,-1.},{+0.,+0.}},
   {{+0.,+0.},{+0.,+1.},{+0.,+0.},{+0.,+0.}},
   {{+0.,+1.},{+0.,+0.},{+0.,+0.},{+0.,+0.}}},

  {{{+0.,+0.},{+0.,+0.},{+0.,+0.},{-1.,+0.}},       
   {{+0.,+0.},{+0.,+0.},{+1.,+0.},{+0.,+0.}},          
   {{+0.,+0.},{+1.,+0.},{+0.,+0.},{+0.,+0.}},             
   {{-1.,+0.},{+0.,+0.},{+0.,+0.},{+0.,+0.}}},

  {{{+0.,+0.},{+0.,+0.},{+0.,-1.},{+0.,+0.}},       
   {{+0.,+0.},{+0.,+0.},{+0.,+0.},{+0.,+1.}},          
   {{+0.,+1.},{+0.,+0.},{+0.,+0.},{+0.,+0.}},             
   {{+0.,+0.},{+0.,-1.},{+0.,+0.},{+0.,+0.}}},

  {{{+1.,+0.},{+0.,+0.},{+0.,+0.},{+0.,+0.}},       
   {{+0.,+0.},{+1.,+0.},{+0.,+0.},{+0.,+0.}},          
   {{+0.,+0.},{+0.,+0.},{+1.,+0.},{+0.,+0.}},             
   {{+0.,+0.},{+0.,+0.},{+0.,+0.},{+1.,+0.}}},

  {{{+1.,+0.},{+0.,+0.},{+0.,+0.},{+0.,+0.}},       
   {{+0.,+0.},{+1.,+0.},{+0.,+0.},{+0.,+0.}},          
   {{+0.,+0.},{+0.,+0.},{-1.,+0.},{+0.,+0.}},             
   {{+0.,+0.},{+0.,+0.},{+0.,+0.},{-1.,+0.}}}};

extern gsl_matrix_complex *PROJECTOR[no_projectors];

void init_projectors();           // initialize the projectors in complex gsl matrices
void destroy_projectors();        // destroy the gsl matrices
void print_projector(const unsigned int index);  // prints a projector given by 'index' to stdout

#endif
