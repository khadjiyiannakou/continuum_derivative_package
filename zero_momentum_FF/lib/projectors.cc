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

#include "../include/projectors.hh"

using namespace std;

gsl_matrix_complex *PROJECTOR[no_projectors] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

gsl_matrix_complex *GAMMA[6] = {NULL, NULL, NULL, NULL, NULL, NULL};


// Creates the gsl matrix projectors
// STATUS: to be tested, check if aux is needed at all, otherwise remove it
void init_projectors()
{
  destroy_projectors(); // clean up, just in case ...
  gsl_complex I;        // complex i
  GSL_SET_COMPLEX(&I, 0., 1.);
  gsl_complex zero;
  GSL_SET_COMPLEX(&zero, 0., 0.);
  gsl_complex one;
  GSL_SET_COMPLEX(&one, 1., 0.);
  // initialize gsl version of gamma matrices
  for (unsigned int i=0; i<6; i++)
  { 
    GAMMA[i] = gsl_matrix_complex_alloc(4,4);
    gsl_matrix_complex_const_view tmp = gsl_matrix_complex_const_view_array(&dirac_matrices[i][0][0][0], 4, 4);
    gsl_matrix_complex_memcpy(GAMMA[i], &tmp.matrix);
  }

  // further definitions, initializations, auxiliary variables
  gsl_complex c;
  gsl_matrix_complex *aux = gsl_matrix_complex_alloc(4,4);
  gsl_matrix_complex *tm_rotation = gsl_matrix_complex_alloc(4,4);
  gsl_matrix_complex_memcpy(tm_rotation, GAMMA[5]); // define the twist rotation 1/sqrt(2) * (1 + i*g5)
  gsl_matrix_complex_scale(tm_rotation, I);         // i*g5
  gsl_matrix_complex_add(tm_rotation, GAMMA[4]);    // add identity
  GSL_SET_COMPLEX(&c, 1./sqrt(2.), 0.);
  gsl_matrix_complex_scale(tm_rotation, c);

  // initialize gsl version of projector matrices
  for (unsigned int i=0; i<no_projectors; i++)
  {
    PROJECTOR[i] = gsl_matrix_complex_alloc(4,4);
    gsl_matrix_complex_set_identity(PROJECTOR[i]); // initialize all projector to spin-1 (identity)
  }
  GSL_SET_COMPLEX(&c, .5, 0.);

  // 1: 0.25*((one+g0)*(one-i*g5*g3))^T  CHECKED
  gsl_matrix_complex_set_identity(aux);
  gsl_blas_zgemm(CblasTrans, CblasTrans, gsl_complex_negative(I), GAMMA[3], GAMMA[5], one, PROJECTOR[1]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[1]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[1]);
  gsl_matrix_complex_scale(PROJECTOR[1], c);

  // 2: 0.25*((one+g0)*(one-i*g5*(g1+g2+g3)))^T  CHECKED
  gsl_matrix_complex_memcpy(aux, GAMMA[1]);
  gsl_matrix_complex_add(aux, GAMMA[2]);
  gsl_matrix_complex_add(aux, GAMMA[3]);
  gsl_blas_zgemm(CblasTrans, CblasTrans, gsl_complex_negative(I), aux, GAMMA[5], one, PROJECTOR[2]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[2]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[2]);
  gsl_matrix_complex_scale(PROJECTOR[2], c);

  // 3: 0.25*((one+g0)*(i*g5*(g1+g2+g3)))^T  CHECKED
  gsl_matrix_complex_memcpy(aux, GAMMA[1]);
  gsl_matrix_complex_add(aux, GAMMA[2]);
  gsl_matrix_complex_add(aux, GAMMA[3]);
  gsl_blas_zgemm(CblasTrans, CblasTrans, I, aux, GAMMA[5], zero, PROJECTOR[3]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[3]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[3]);
  gsl_matrix_complex_scale(PROJECTOR[3], c);

  // 4: 0.25*((one+g0)*(i*g5*g1))^T  CHECKED
  gsl_matrix_complex_memcpy(aux, GAMMA[1]);
  gsl_blas_zgemm(CblasTrans, CblasTrans, I, aux, GAMMA[5], zero, PROJECTOR[4]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[4]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[4]);
  gsl_matrix_complex_scale(PROJECTOR[4], c);

  // 5: 0.25*((one+g0)*(i*g5*g2))^T  CHECKED
  gsl_matrix_complex_memcpy(aux, GAMMA[2]);
  gsl_blas_zgemm(CblasTrans, CblasTrans, I, aux, GAMMA[5], zero, PROJECTOR[5]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[5]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[5]);
  gsl_matrix_complex_scale(PROJECTOR[5], c);

  // 6: 0.25*((one+g0)*(i*g5*g3))^T  CHECKED
  gsl_matrix_complex_memcpy(aux, GAMMA[3]);
  gsl_blas_zgemm(CblasTrans, CblasTrans, I, aux, GAMMA[5], zero, PROJECTOR[6]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[6]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[6]);
  gsl_matrix_complex_scale(PROJECTOR[6], c);

  // 13: 0.5*(one+g0)^T  CHECKED
  gsl_matrix_complex_memcpy(aux,PROJECTOR[13]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[13]);

  // 14: 0.5*((one+g0)*g5)^T ("NON-STANDARD" !!)  CHECKED
  gsl_matrix_complex_set_identity(aux);
  gsl_matrix_complex_add(aux, GAMMA[0]);
  gsl_blas_zgemm(CblasTrans, CblasTrans, c, GAMMA[5], aux, zero, PROJECTOR[14]);

  // 15: 0.25*((one+g0)*(i*g5*g2))^T (same as 5)  CHECKED
  gsl_matrix_complex_memcpy(aux, GAMMA[2]);
  gsl_blas_zgemm(CblasTrans, CblasTrans, I, aux, GAMMA[5], zero, PROJECTOR[15]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[15]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[15]);
  gsl_matrix_complex_scale(PROJECTOR[15], c);

  // 16: 0.25*((one+g0)*(i*g5*g1))^T (same as 4)  CHECKED
  gsl_matrix_complex_memcpy(aux, GAMMA[1]);
  gsl_blas_zgemm(CblasTrans, CblasTrans, I, aux, GAMMA[5], zero, PROJECTOR[16]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[16]);
  gsl_blas_zgemm(CblasNoTrans, CblasTrans, c, aux, GAMMA[0], c, PROJECTOR[16]);
  gsl_matrix_complex_scale(PROJECTOR[16], c);

  // 20: 0.5*(one+(g1+g2+g3))^T  CHECKED
  gsl_matrix_complex_memcpy(aux, PROJECTOR[20]);
  gsl_matrix_complex_add(aux, GAMMA[1]);
  gsl_matrix_complex_add(aux, GAMMA[2]);
  gsl_matrix_complex_add(aux, GAMMA[3]);
  gsl_blas_zgemm(CblasTrans, CblasNoTrans, c, aux, GAMMA[4], zero, PROJECTOR[20]);


  // 21: 0.25*(one)^T  CHECKED
  gsl_matrix_complex_memcpy(aux, PROJECTOR[21]);
  gsl_blas_zgemm(CblasTrans, CblasNoTrans, c, aux, GAMMA[4], zero, PROJECTOR[21]);
  gsl_matrix_complex_scale(PROJECTOR[21], c);

  // 22: 0.25*((one+g0)*(i*g5*g1)) NOT in tm-basis!  CHECKED
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, I, GAMMA[5], GAMMA[1], zero, PROJECTOR[22]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[22]);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, c, GAMMA[0], aux, c, PROJECTOR[22]);
  gsl_matrix_complex_scale(PROJECTOR[22], c);

  // 23: 0.25*((one+g0)*(i*g5*g2)) NOT in tm-basis!  CHECKED
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, I, GAMMA[5], GAMMA[2], zero, PROJECTOR[23]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[23]);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, c, GAMMA[0], aux, c, PROJECTOR[23]);
  gsl_matrix_complex_scale(PROJECTOR[23], c);

  // 24: 0.25*((one+g0)*(i*g5*g3)) NOT in tm-basis!  CHECKED
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, I, GAMMA[5], GAMMA[3], zero, PROJECTOR[24]);
  gsl_matrix_complex_memcpy(aux, PROJECTOR[24]);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, c, GAMMA[0], aux, c, PROJECTOR[24]);
  gsl_matrix_complex_scale(PROJECTOR[24], c);

  // 25: 0.25*g0  CHECKED
  gsl_matrix_complex_memcpy(PROJECTOR[25], GAMMA[0]);

  // 26: g1  CHECKED
  gsl_matrix_complex_memcpy(PROJECTOR[26], GAMMA[1]);

  // 27: g2  CHECKED
  gsl_matrix_complex_memcpy(PROJECTOR[27], GAMMA[2]);

  // 28: g3  CHECKED
  gsl_matrix_complex_memcpy(PROJECTOR[28], GAMMA[3]);

  // 29: id = 0.25*g4  CHECKED
  gsl_matrix_complex_memcpy(PROJECTOR[29], GAMMA[4]);

  // 30: g5  CHECKED
  gsl_matrix_complex_memcpy(PROJECTOR[30], GAMMA[5]);

  for (unsigned int i=0; i<22; i++) // rotate to tm-basis
  {
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, tm_rotation, PROJECTOR[i], zero, aux);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, aux, tm_rotation, zero, PROJECTOR[i]);
  }
  for (unsigned int i=25; i<no_projectors; i++) // rotate to tm-basis, leave out [22..24]
  {
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, tm_rotation, PROJECTOR[i], zero, aux);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, aux, tm_rotation, zero, PROJECTOR[i]);
  }
  gsl_matrix_complex_free(aux);
  gsl_matrix_complex_free(tm_rotation);
  return;
}


// Destroys the gsl matrix projectors
// STATUS: done
void destroy_projectors()
{
  for (unsigned int i=0; i<6; i++)
  {
    gsl_matrix_complex_free(GAMMA[i]);
    GAMMA[i] = NULL;
  }
  for (unsigned int i=0; i<no_projectors; i++)
  {
    gsl_matrix_complex_free(PROJECTOR[i]);
    PROJECTOR[i] = NULL;
  }
  return;
}


// Print a projector to stdout
// STATUS: done
void print_projector(const unsigned int index)
{
  cout << endl;
  if (index<no_projectors)
  {
    if (PROJECTOR[index]!=NULL)
    {
      for (unsigned int i=0; i<4; i++)
      {
        for (unsigned int j=0; j<4; j++)
        {
          cout << showpos << setw(5) << GSL_REAL(gsl_matrix_complex_get(PROJECTOR[index], i, j)) << " " << setw(5) << GSL_IMAG(gsl_matrix_complex_get(PROJECTOR[index], i, j)) << "  ";
        }
        cout << endl;
      }
    }
    else
    {
      cout << "NULL (unititalized)" << endl;
    }
  }
  else
  {
    cout << "FAIL: Illegal projector index " << index << ". Valid indices are [0..." << no_projectors-1 << "]." << endl;
  }
  return; 
}

