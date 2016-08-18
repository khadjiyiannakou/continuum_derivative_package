// AUTHOR: Konstantin Ottnad
//
// DATE: 20140707

#include "../include/resampling.hh"

using namespace std;

// Version for real gsl vectors
// STATUS: done
int bootstrap(gsl_vector *IN, gsl_vector *&OUT, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed)
{
  if (IN!=NULL)
  {
    const unsigned int no_blocks = IN->size / blocksize;
    unsigned int sample[no_blocks];
    double BLOCKS[no_blocks];
    srand(seed);
    OUT = gsl_vector_calloc(no_samples); // initialize to zero!
    for (unsigned int block=0; block<no_blocks; block++) // perform blocking
    {
      BLOCKS[block] = 0.; // initialize to zero!
      for (unsigned int B=0; B<blocksize; B++) BLOCKS[block] += gsl_vector_get(IN, block * blocksize + B);
      BLOCKS[block]/=blocksize;
    }
    for (unsigned int boot=0; boot<no_samples; boot++) // perform bootstrapping
    {
      if (boot==0)
      {
        for (unsigned int block=0; block<no_blocks; block++) sample[block]=block; // use original data for zeroth sample
      }
      else
      {
        for (unsigned int block=0; block<no_blocks; block++) sample[block]=rand() % no_blocks; // use random sampling with replacement otherwise
      }
      for (unsigned int block=0; block<no_blocks; block++) gsl_vector_set(OUT, boot, gsl_vector_get(OUT, boot) + BLOCKS[sample[block]]);
    }
    gsl_vector_scale(OUT, 1./no_blocks); // gauge averaging -- need to divide each sample by the number of blocks
  }
  else
  {
    return -1; // no input data
  }
  return 0;
}


int bootstrap(gsl_vector *IN, gsl_vector *&OUT, double &ERR, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed)
{
  if (bootstrap(IN, OUT, no_samples, blocksize, seed)) return -1;
  ERR = gsl_stats_sd(OUT->data, OUT->stride, no_samples); // calculate the error
  return 0;
}


// Version for complex gsl vectors
// STATUS: done
int bootstrap(gsl_vector_complex *IN, gsl_vector_complex *&OUT, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed)
{
  if (IN!=NULL)
  {
    const unsigned int no_blocks = IN->size / blocksize;
    unsigned int sample[no_blocks];
    gsl_complex BLOCKS[no_blocks];
    srand(seed);
    OUT = gsl_vector_complex_calloc(no_samples); // initialize to zero!
    for (unsigned int block=0; block<no_blocks; block++) // perform blocking
    {
      BLOCKS[block] = gsl_complex_rect(0.,0.); // initialize to zero!
      for (unsigned int B=0; B<blocksize; B++) BLOCKS[block] = gsl_complex_add(BLOCKS[block], gsl_vector_complex_get(IN, block * blocksize + B));
      BLOCKS[block] = gsl_complex_div(BLOCKS[block], gsl_complex_rect(blocksize,0.));
    }
    for (unsigned int boot=0; boot<no_samples; boot++) // perform bootstrapping
    {
      if (boot==0)
      {
        for (unsigned int block=0; block<no_blocks; block++) sample[block]=block; // use original data for zeroth sample
      }
      else
      {
        for (unsigned int block=0; block<no_blocks; block++) sample[block]=rand() % no_blocks; // use random sampling with replacement otherwise
      }
      for (unsigned int block=0; block<no_blocks; block++) gsl_vector_complex_set(OUT, boot, gsl_complex_add(gsl_vector_complex_get(OUT, boot), BLOCKS[sample[block]]));
    }
    gsl_vector_complex_scale(OUT, gsl_complex_rect(1./no_blocks,0.)); // gauge averaging -- need to divide each sample by the number of blocks
  }
  else
  {
    return -1; // no input data
  }
  return 0;
}


int bootstrap(gsl_vector_complex *IN, gsl_vector_complex *&OUT, gsl_complex &ERR, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed)
{
  if (bootstrap(IN, OUT, no_samples, blocksize, seed)) return -1;
  gsl_vector_const_view Re = gsl_vector_complex_const_real(OUT);
  gsl_vector_const_view Im = gsl_vector_complex_const_imag(OUT);
  ERR = gsl_complex_rect(gsl_stats_sd(Re.vector.data, Re.vector.stride, no_samples), gsl_stats_sd(Im.vector.data, Im.vector.stride, no_samples)); // calculate the error
  return 0;
}


// Version for real gsl matrices
// STATUS: done
int bootstrap(gsl_matrix *IN, gsl_matrix *&OUT, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed)
{
  if (IN!=NULL)
  {
    const unsigned int no_blocks = IN->size1 / blocksize;
    unsigned int sample[no_blocks];
    gsl_vector *BLOCKS[no_blocks];
    srand(seed);
    OUT = gsl_matrix_calloc(no_samples, IN->size2); // initialize to zero!
    for (unsigned int block=0; block<no_blocks; block++) // perform blocking
    {
      BLOCKS[block] = gsl_vector_calloc(IN->size2); // initialize to zero!
      for (unsigned int B=0; B<blocksize; B++)
       {
         gsl_vector_const_view tmp = gsl_matrix_const_row(IN, block * blocksize + B);
         gsl_vector_add(BLOCKS[block], &tmp.vector);
       }
      gsl_vector_scale(BLOCKS[block], 1./blocksize);
    }
    for (unsigned int boot=0; boot<no_samples; boot++) // perform bootstrapping
    {
      if (boot==0)
      {
        for (unsigned int block=0; block<no_blocks; block++) sample[block]=block; // use original data for zeroth sample
      }
      else
      {
        for (unsigned int block=0; block<no_blocks; block++) sample[block]=rand() % no_blocks; // use random sampling with replacement otherwise
      }
      for (unsigned int block=0; block<no_blocks; block++)
      {
        gsl_vector_view tmp = gsl_matrix_row(OUT, boot);
        gsl_vector_add(&tmp.vector, BLOCKS[sample[block]]);
      }
    }
    gsl_matrix_scale(OUT, 1./no_blocks); // gauge averaging -- need to divide each sample by the number of blocks    
    for (unsigned int block=0; block<no_blocks; block++) gsl_vector_free(BLOCKS[block]); // clean up the BLOCKS
  }
  else
  {
    return -1; // no input data
  }
  return 0; // everything is ok.
}


int bootstrap(gsl_matrix *IN, gsl_matrix *&OUT, gsl_vector *&ERR, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed)
{
  if (bootstrap(IN, OUT, no_samples, blocksize, seed)) return -1;
  ERR = gsl_vector_calloc(IN->size2);
  for (unsigned int i=0; i<IN->size2; i++) // calculate the errors
  {
    gsl_vector_set(ERR, i, gsl_stats_sd(gsl_matrix_column(OUT, i).vector.data, gsl_matrix_column(OUT, i).vector.stride, no_samples));
  }
  return 0;
}



// Version for complex gsl matrices
// STATUS: done
int bootstrap(gsl_matrix_complex *IN, gsl_matrix_complex *&OUT, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed)
{
  if (IN!=NULL)
  {
    const unsigned int no_blocks = IN->size1 / blocksize;
    unsigned int sample[no_blocks];
    gsl_vector_complex *BLOCKS[no_blocks];
    srand(seed);
    OUT = gsl_matrix_complex_calloc(no_samples, IN->size2); // initialize to zero!
    for (unsigned int block=0; block<no_blocks; block++) // perform blocking
    {
      BLOCKS[block] = gsl_vector_complex_calloc(IN->size2); // initialize to zero!
      for (unsigned int B=0; B<blocksize; B++)
      {
        gsl_vector_complex_const_view tmp = gsl_matrix_complex_const_row(IN, block * blocksize + B);
        gsl_vector_complex_add(BLOCKS[block], &tmp.vector);
      }
      gsl_vector_complex_scale(BLOCKS[block], gsl_complex_rect(1./blocksize, 0.));
    }

    for (unsigned int boot=0; boot<no_samples; boot++) // perform bootstrapping
    {
      if (boot==0)
      {
        for (unsigned int block=0; block<no_blocks; block++) sample[block]=block; // use original data for zeroth sample
      }
      else
      {
        for (unsigned int block=0; block<no_blocks; block++) sample[block]=rand() % no_blocks; // use random sampling with replacement otherwise
      }
      for (unsigned int block=0; block<no_blocks; block++)
      {
        gsl_vector_complex_view tmp = gsl_matrix_complex_row(OUT, boot);
        gsl_vector_complex_add(&tmp.vector, BLOCKS[sample[block]]);
      }
    }
    gsl_matrix_complex_scale(OUT, gsl_complex_rect(1./no_blocks,0.)); // gauge averaging -- need to divide each sample by the number of blocks
    for (unsigned int block=0; block<no_blocks; block++) gsl_vector_complex_free(BLOCKS[block]); // clean up the BLOCKS 
  }
  else
  {
    return -1; // no input data
  }
  return 0; // everything is ok.
}


int bootstrap(gsl_matrix_complex *IN, gsl_matrix_complex *&OUT, gsl_vector_complex *&ERR, const unsigned int no_samples, const unsigned int blocksize, const unsigned int seed)
{
  if (bootstrap(IN, OUT, ERR, no_samples, blocksize, seed)) return -1;
  ERR = gsl_vector_complex_calloc(IN->size2);
  for (unsigned int i=0; i<IN->size2; i++) // calculate the errors
  {
    gsl_vector_complex_const_view tmp = gsl_matrix_complex_const_column(OUT, i);
    gsl_vector_const_view Re = gsl_vector_complex_const_real(&tmp.vector);
    gsl_vector_const_view Im = gsl_vector_complex_const_imag(&tmp.vector);
    gsl_vector_complex_set(ERR, i, gsl_complex_rect(gsl_stats_sd(Re.vector.data, Re.vector.stride, no_samples), gsl_stats_sd(Im.vector.data, Im.vector.stride, no_samples)));
  }
  return 0;
}


// Version for real gsl vectors
// STATUS: done
int jackknife(gsl_vector *IN, gsl_vector *&OUT, const unsigned int blocksize)
{
  if (IN!=NULL)
  {
    const unsigned int no_blocks = IN->size / blocksize;
    const unsigned int no_samples = no_blocks + 1; // need one more for the original data
    double BLOCKS[no_blocks];
    OUT = gsl_vector_calloc(no_samples); // initialize to zero!
    for (unsigned int block=0; block<no_blocks; block++) // perform blocking
    {
      BLOCKS[block] = 0; // initialize to zero!
      for (unsigned int B=0; B<blocksize; B++) BLOCKS[block] += gsl_vector_get(IN, block * blocksize + B);
      BLOCKS[block]/=blocksize;
    }
    for (unsigned int boot=0; boot<no_samples; boot++) // perform jackknifing
    {
      if (boot==0)
      {
        for (unsigned int block=0; block<no_blocks; block++) gsl_vector_set(OUT, boot, gsl_vector_get(OUT, boot) + BLOCKS[block]); // use original data for zeroth sample
        gsl_vector_set(OUT, boot, gsl_vector_get(OUT, boot) / no_blocks); // gauge averaging -- need to divide each sample by the number of blocks in it -- this is different for the original data!!
      }
      else
      {
        for (unsigned int block=0; block<no_blocks; block++) 
        {
          if (block!=(boot-1)) gsl_vector_set(OUT, boot, gsl_vector_get(OUT, boot) + BLOCKS[block]); // mind the shift in boot, because we use boot==0 for the oringal data. -> leave the boot-th block out
        }
        gsl_vector_set(OUT, boot, gsl_vector_get(OUT, boot) / (no_blocks-1)); // gauge averaging -- need to divide each sample by the number of blocks in it
      }
    }
  }
  else
  {
    return -1; // no input data
  }
  return 0;
}


int jackknife(gsl_vector *IN, gsl_vector *&OUT, double &ERR, const unsigned int blocksize)
{
  if (jackknife(IN, OUT, blocksize)) return -1;
  const unsigned int no_samples = IN->size / blocksize + 1; // need one more for the original data
  ERR = sqrt((no_samples-1)-1) * gsl_stats_sd(gsl_vector_subvector(OUT, 1, no_samples-1).vector.data, gsl_vector_subvector(OUT, 1, no_samples-1).vector.stride, no_samples-1); // calculate the error
  return 0;
}


// Version for complex gsl vectors
// STATUS: done
int jackknife(gsl_vector_complex *IN, gsl_vector_complex *&OUT, const unsigned int blocksize)
{
  if (IN!=NULL)
  {
    const unsigned int no_blocks = IN->size / blocksize;
    const unsigned int no_samples = no_blocks + 1; // need one more for the original data
    gsl_complex BLOCKS[no_blocks];
    OUT = gsl_vector_complex_calloc(no_samples); // initialize to zero!    
    for (unsigned int block=0; block<no_blocks; block++) // perform blocking
    {
      BLOCKS[block] = gsl_complex_rect(0.,0.); // initialize to zero!
      for (unsigned int B=0; B<blocksize; B++) BLOCKS[block] = gsl_complex_add(BLOCKS[block], gsl_vector_complex_get(IN, block * blocksize + B));
      BLOCKS[block] = gsl_complex_div(BLOCKS[block], gsl_complex_rect(blocksize,0.));
    }
    
    for (unsigned int boot=0; boot<no_samples; boot++) // perform jackknifing
    {
      if (boot==0)
      {
        for (unsigned int block=0; block<no_blocks; block++) gsl_vector_complex_set(OUT, boot, gsl_complex_add(gsl_vector_complex_get(OUT, boot), BLOCKS[block])); // use original data for zeroth sample
        gsl_vector_complex_set(OUT, boot, gsl_complex_div(gsl_vector_complex_get(OUT, boot), gsl_complex_rect(1./no_blocks, 0.))); // gauge averaging -- need to divide each sample by the number of blocks in it -- this is different for the original data!!
      }
      else
      {
        for (unsigned int block=0; block<no_blocks; block++) 
        {
          if (block!=(boot-1)) gsl_vector_complex_set(OUT, boot, gsl_complex_add(gsl_vector_complex_get(OUT, boot), BLOCKS[block]));   // mind the shift in boot, because we use boot==0 for the oringal data. -> leave the boot-th block out
        }
        gsl_vector_complex_set(OUT, boot, gsl_complex_div(gsl_vector_complex_get(OUT, boot), gsl_complex_rect(1./(no_blocks-1), 0.))); // gauge averaging -- need to divide each sample by the number of blocks in it
      }
    }
  }
  else
  {
    return -1; // no input data
  }
  return 0;
}


int jackknife(gsl_vector_complex *IN, gsl_vector_complex *&OUT, gsl_complex &ERR, const unsigned int blocksize)
{
  if (jackknife(IN, OUT, blocksize)) return -1;
  const unsigned int no_samples = IN->size / blocksize + 1; // need one more for the original data
  gsl_vector_complex_const_view tmp = gsl_vector_complex_const_subvector(OUT, 1, no_samples-1);
  gsl_vector_const_view Re = gsl_vector_complex_const_real(&tmp.vector);
  gsl_vector_const_view Im = gsl_vector_complex_const_imag(&tmp.vector);
  ERR = gsl_complex_rect(gsl_stats_sd(Re.vector.data, Re.vector.stride, no_samples-1), gsl_stats_sd(Im.vector.data, Im.vector.stride, no_samples-1)); // calculate the error
  ERR = gsl_complex_mul_real(ERR, sqrt((no_samples-1)-1));
  return 0;
}


// Version for real gsl matrices
// STATUS: done
int jackknife(gsl_matrix *IN, gsl_matrix *&OUT, const unsigned int blocksize)
{
  if (IN!=NULL)
  {
    const unsigned int no_blocks = IN->size1 / blocksize;
    const unsigned int no_samples = no_blocks + 1; // need one more for the original data
    gsl_vector *BLOCKS[no_blocks];
    OUT = gsl_matrix_calloc(no_samples, IN->size2); // initialize to zero!
    for (unsigned int block=0; block<no_blocks; block++) // perform blocking
    {
      BLOCKS[block] = gsl_vector_calloc(IN->size2); // initialize to zero!
      for (unsigned int B=0; B<blocksize; B++)
      {
        gsl_vector_const_view tmp = gsl_matrix_const_row(IN, block * blocksize + B);
        gsl_vector_add(BLOCKS[block], &tmp.vector);
      }
      gsl_vector_scale(BLOCKS[block], 1./blocksize);
    }
    
    for (unsigned int boot=0; boot<no_samples; boot++) // perform jackknifing
    {
      if (boot==0)
      {
        gsl_vector_view tmp = gsl_matrix_row(OUT, boot);
        for (unsigned int block=0; block<no_blocks; block++) gsl_vector_add(&tmp.vector, BLOCKS[block]); // use original data for zeroth sample
        gsl_vector_scale(&tmp.vector, 1./no_blocks); // gauge averaging -- need to divide each sample by the number of blocks in it -- this is different for the original data!!
      }
      else
      {
        gsl_vector_view tmp = gsl_matrix_row(OUT, boot);
        for (unsigned int block=0; block<no_blocks; block++) if (block!=(boot-1)) gsl_vector_add(&tmp.vector, BLOCKS[block]);   // mind the shift in boot, because we use boot==0 for the oringal data. -> leave the boot-th block out
        gsl_vector_scale(&tmp.vector, 1./(no_blocks-1)); // gauge averaging -- need to divide each sample by the number of blocks in it
      }
    }
    for (unsigned int block=0; block<no_blocks; block++) gsl_vector_free(BLOCKS[block]); // clean up the BLOCKS
  }
  else
  {
    return -1; // no input data
  }
  return 0;
}


int jackknife(gsl_matrix *IN, gsl_matrix *&OUT, gsl_vector *&ERR, const unsigned int blocksize)
{
  if (jackknife(IN, OUT, blocksize)) return -1;
  const unsigned int no_samples = IN->size1 / blocksize + 1; // need one more for the original data
  ERR = gsl_vector_calloc(IN->size2);
  for (unsigned int i=0; i<IN->size2; i++) // calculate the errors
  {
    gsl_vector_set(ERR, i, sqrt((no_samples-1)-1) * gsl_stats_sd(gsl_matrix_subcolumn(OUT, i, 1, no_samples-1).vector.data, gsl_matrix_subcolumn(OUT, i, 1, no_samples-1).vector.stride, no_samples-1));
  }
  return 0;
}


// Version for complex gsl matrices
// STATUS: done
int jackknife(gsl_matrix_complex *IN, gsl_matrix_complex *&OUT, const unsigned int blocksize)
{
  if (IN!=NULL)
  {
    const unsigned int no_blocks = IN->size1 / blocksize;
    const unsigned int no_samples = no_blocks + 1; // need one more for the original data
    gsl_vector_complex *BLOCKS[no_blocks];
    OUT = gsl_matrix_complex_calloc(no_samples, IN->size2); // initialize to zero!
    for (unsigned int block=0; block<no_blocks; block++) // perform blocking
    {
      BLOCKS[block] = gsl_vector_complex_calloc(IN->size2); // initialize to zero!
      for (unsigned int B=0; B<blocksize; B++)
      {
        gsl_vector_complex_const_view tmp = gsl_matrix_complex_const_row(IN, block * blocksize + B);
        gsl_vector_complex_add(BLOCKS[block], &tmp.vector);
      }
      gsl_vector_complex_scale(BLOCKS[block], gsl_complex_rect(1./blocksize, 0.));
    }
    
    for (unsigned int boot=0; boot<no_samples; boot++) // perform jackknifing
    {
      gsl_vector_complex_view tmp = gsl_matrix_complex_row(OUT, boot);
      if (boot==0)
      {
        for (unsigned int block=0; block<no_blocks; block++) gsl_vector_complex_add(&tmp.vector, BLOCKS[block]); // use original data for zeroth sample
        gsl_vector_complex_scale(&tmp.vector, gsl_complex_rect(1./(no_blocks),0.)); // gauge averaging -- need to divide each sample by the number of blocks in it -- this is different for the original data!!
      }
      else
      {
        for (unsigned int block=0; block<no_blocks; block++)
        {
          if (block!=(boot-1)) gsl_vector_complex_add(&tmp.vector, BLOCKS[block]);   // mind the shift in boot, because we use boot==0 for the oringal data. -> leave the boot-th block out
        }
        gsl_vector_complex_scale(&tmp.vector, gsl_complex_rect(1./(no_blocks-1),0.)); // gauge averaging -- need to divide each sample by the number of blocks in it
      }
    }
    for (unsigned int block=0; block<no_blocks; block++) gsl_vector_complex_free(BLOCKS[block]); // clean up the BLOCKS
  }
  else
  {
    return -1; // no input data
  }
  return 0;
}


int jackknife(gsl_matrix_complex *IN, gsl_matrix_complex *&OUT, gsl_vector_complex *&ERR, const unsigned int blocksize)
{
  if (jackknife(IN, OUT, blocksize)) return -1;
  const unsigned int no_samples = IN->size1 / blocksize + 1; // need one more for the original data
  ERR = gsl_vector_complex_calloc(IN->size2);
  for (unsigned int i=0; i<IN->size2; i++) // calculate the errors
  {
    gsl_vector_complex_const_view tmp = gsl_matrix_complex_const_subcolumn(OUT, i, 1, no_samples-1);
    gsl_vector_const_view Re = gsl_vector_complex_const_real(&tmp.vector);
    gsl_vector_const_view Im = gsl_vector_complex_const_imag(&tmp.vector);
    gsl_vector_complex_set(ERR, i, gsl_complex_rect(gsl_stats_sd(Re.vector.data, Re.vector.stride, no_samples-1), gsl_stats_sd(Im.vector.data, Im.vector.stride, no_samples-1)));
  }
  gsl_vector_complex_scale(ERR, gsl_complex_rect(sqrt((no_samples-1)-1), 0.));
  return 0;
}


// Returns a pair of independent, normal distributed random numbers
// STATUS: done
void polar_method(double *z)
{
  double u1, u2, sum;
  do
  {
    u1 = 1. - 2. * double(random()) / double(RAND_MAX);
    u2 = 1. - 2. * double(random()) / double(RAND_MAX);
    sum = u1 * u1 + u2 * u2;
  }
  while ((sum==0.0)||(sum>=1.0));
  z[0] = sqrt(-2. * log(sum) / sum);
  z[1] = u2 * z[0];
  z[0] *= u1;
  return;
}


// Performs resampling (normal distributed)
// STATUS : done
void resampling(const double mean, const double err, const unsigned int no_samples, double * results)
{
  double z[2] = {0.,0.};
  for (unsigned int b=1; b<no_samples; b++)
  {
    if (!(b%2))  // we also generated two random numbers by each call of polar_method()
    {
      z[0] = err * z[1];
    }
    else
    {
      polar_method(z);
      z[0] *= err; 
    }
    results[b] = z[0] + mean;
  }
  results[0]=mean; // zeroth sample contains original data !!!!
  return;
}
