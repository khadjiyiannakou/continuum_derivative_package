// AUTHOR: Konstantin Ottnad
//
// DATE: 20131212
//
// Implements:
// - classes to read (in position space), write (in momentum space) and manipulate correlation function
// - perform Fourier trafo for given momenta list, calculate derivative with respect to dq_i (i.e. spatial components of momentum)
//
// Open issues:
// - in Fourier transform: if not running over the entire lattice e.g. r<L: how to define the momenta, i.e. is it still 2*pi*n/L or rather 2*pi*n/r ?

#include "../include/correlator.hh"
#include "../include/io.hh"
#include "../include/projectors.hh"

using namespace std;


// translate a nucleon isospin string into a corresponding 'isospin_nucleon' enum value (-1 if invalid)
// STATUS: done
int get_isospin_nucleon(string isospin)
{
  if (isospin=="ppm")
  {
    return ppm;
  }
  else
  {
    if (isospin=="pmm")
    {
      return pmm;
    }
  }
  return -1; // return negative value for further error handling if the input string does not denote a valid nucleon isospin
}


// translate a nucleon spin string into a corresponding 'spin_nucleon' enum value (-1 if invalid)
// STATUS: done
int get_spin_nucleon(string operators)
{
  if (operators=="1_1")
  {
    return spin_1_1;
  }
  else
  {
    if (operators=="1_2")
    {
      return spin_1_2;
    }
    else
    {
      if (operators=="2_1")
      {
	return spin_2_1;
      }
      else
      {
	if (operators=="2_2")
	{
	  return spin_2_2;
	}
      }
    }
  }
  return -1; // return negative value for further error handling if the input string does not denote valid interpolating operators
}


// constructor for base class 'xcorr'
// STATUS: done
xcorr::xcorr(string xfname, string pfname, string momentalist_fname)
{
  error = 0;
  L = 0;
  T = 0;
  mode = undefined;
  no_t = 0;
  no_momenta = 0;
  header_lines = 0;
  t_values = NULL;
  for (int i=0; i<4; i++) x_source[i] = 0;
  filename = pfname;
  corr = NULL;
  x_corr = NULL;
  y_corr = NULL;
  z_corr = NULL;
  cout << endl << "Reading momenta list ... ";
  flush(cout);
  mreader = new mixed_text_reader (momentalist_fname, 1, 0, format_momenta); // read momenta list from file ...
  if (mreader->get_error_id()==0)
  {
    no_momenta = mreader->get_lines();
    cout << "done";
    flush(cout);
  }
  else
  {
    cout << mreader->get_error() << " Projecting to zero momentum only ... ";
    flush(cout);
    no_momenta = 1; // ... or project to zero momentum if there is any problem
    delete mreader->index;
    mreader->index = new int [3];
    mreader->index[0] = 0;
    mreader->index[1] = 0;
    mreader->index[2] = 0;
    error = 1;
  }
  momenta = mreader->index;  // this is just for convenience...
  reader = new binary_reader (xfname, "begin-header;", "end-header;"); // read the correlator input data in position space and header info ...
  header_lines = (*reader).read_header();          // read header -> parse it separately for each derived class !!
  int64_t buffer_size = (* reader).read_data();    // read data
  if (buffer_size==0) // check for any error
  {
    error = 2; 
    cout << (* reader).get_error();
  }
  if (get_endianess()!=bigEndian) // check for correct endianess
  {
    cout << endl << "swapping data endianess to little endian ... ";
    swap_byte_order_8((* reader).data, buffer_size/8);
    cout << "done";
  }
  return;
}


// destructor for base class 'xcorr'
// STATUS: done 
xcorr::~xcorr()
{
  delete reader;
  delete mreader;
  if (t_values!=NULL) delete [] t_values;
  if (corr!=NULL) delete [] corr;
  if (x_corr!=NULL) delete [] x_corr;
  if (y_corr!=NULL) delete [] y_corr;
  if (z_corr!=NULL) delete [] z_corr;
  cout << "done";
  return;
}


// returns the current value of L
// STATUS: done
int xcorr::get_L()
{
  return L;
}

// returns the current value of T
// STATUS: done
int xcorr::get_T()
{
  return T;
}


// Converts fftw_complex * -> complex<double>. Both arrays of size 'size' must be allocated by the caller. No error checking is performed
// STATUS: done
inline void xcorr::convert_fftw_complex(fftw_complex *in, complex<double> *out, int64_t size)
{
  for (int64_t i=0; i<size; i++) out[i] = complex<double>(in[i].re, in[i].im);
  return;
}


// Converts complex<double> * -> fftw_complex * . Both arrays of size 'size' must be allocated by the caller. No error checking is performed
// STATUS: done
inline void xcorr::convert_fftw_complex(complex<double> *in, fftw_complex *out, int64_t size)
{
  for (int64_t i=0; i<size; i++)
  {
    out[i].re = in[i].real();
    out[i].im = in[i].imag();
  }
  return;
}


// corr file header parser for base class xcorr; parses for string that are always present, regardless of the actual corr being 2pt or 3pt
// STATUS: done
int xcorr::parse_header()
{
  cout << endl << "Parsing header for file " << (*reader).get_filename() << " ... "; flush(cout);
  char dummy;
  stringstream stream;
  stream.str((*reader).parse_header("lattice dimensions (x,y,z,t) = "));
  if (stream.str()!="")
  {
    for (unsigned int i=0; i<3; i++)
    {
      stream >> L;
      stream >> dummy;
    }
    stream >> T;
  }
  else
  {
    cout << endl << "FAIL: could not determine lattice dimension from header in correlator file";
    return -1;
  }

  stream.str((*reader).parse_header("source position (x,y,z,t) = "));
  if (stream.str()!="")
  {
    for (int i=0; i<4; i++) // first three are spatial, last one is t
    {
      stream >> x_source[i];
      stream >> dummy;
    }
  }
  else
  {
    cout << endl << "FAIL: could not determine source position from header in correlator file";
    return -1;
  }
  return 0;
}


// calculates Fourier trafo x->q for given momentalist and continuum derivative with respect to q_i (i.e. spatial components, factor "x_i" in the FT)
// STATUS: done
void xcorr::FT_Dq_cont()
{
  FT_Dq_cont(L); // call the class specific FT method for the special case of radius=L
  return;
}


// calculates Fourier trafo x->q for given momentalist and lattice derivative with respect to q_i (i.e. spatial components, factor "sin(x_i)" in the FT)
// STATUS: done
void xcorr::FT_Dq_latt()
{
  FT_Dq_latt(L); // call the class specific FT method for the special case of radius=L
  return;
}


// deletes the corr, x_corr, y_corr and z_corr and resets the FT mode
// STATUS: done
void xcorr::reset()
{
  if (corr!=NULL) delete [] corr;
  if (x_corr!=NULL) delete [] x_corr;
  if (y_corr!=NULL) delete [] y_corr;
  if (z_corr!=NULL) delete [] z_corr;
  corr=NULL; // this is required to avoid double-free corruption
  x_corr=NULL;
  y_corr=NULL;
  z_corr=NULL;
  mode = undefined; // reset the FT mode -> discard any radius information
  radius = L; // just to be safe; this is the default value, that also tells every FT() method to re-allocate and start from scratch
  return;
}


// writes the correlators in momentum space to the file given by pfilename, if pfilename="" write to stdout. 'prefix' allows to give different file names
// to 2pt and 3pt functions. Calls the class specific, virtual write_file() method.
// STATUS: done
int xcorr::write_results(string prefix)
{
  cout << endl << "Writing results for r = " << radius << " ... "; flush(cout);
  stringstream stream;
  stream.str("");
  stream << prefix << ".r" << radius << "." << filename;
  if (corr!=NULL) write_file(stream.str(), corr);
  if (mode==cont) // check if FT mode is set for continuum derivatives -- not additional tag
  {
    if (x_corr!=NULL) write_file("Dq_x." + stream.str(), x_corr);
    if (y_corr!=NULL) write_file("Dq_y." + stream.str(), y_corr);
    if (z_corr!=NULL) write_file("Dq_z." + stream.str(), z_corr);
  }
  if (mode==latt) // check if FT mode is set for lattice derivatives -- append _latt to Dq_{x,y,z}
  {
    if (x_corr!=NULL) write_file("Dq_x_latt." + stream.str(), x_corr);
    if (y_corr!=NULL) write_file("Dq_y_latt." + stream.str(), y_corr);
    if (z_corr!=NULL) write_file("Dq_z_latt." + stream.str(), z_corr);
  }
  if ((mode==blockxyz)||(mode==blockxyzsingle)||(mode==FFTxyz)||(mode==FFTxyzsingle))
  {
    // need to create (fake) momentum lists (i.e. (x, 0, 0), (0, y, 0), (0, 0, z))
    int temp_no_momenta = no_momenta;
    no_momenta = get_L();        // also change the number of momenta to L
    int *temp_momenta = momenta; // save original momentum list
    int *momenta_x = new int [3*no_momenta];
    int *momenta_y = new int [3*no_momenta];
    int *momenta_z = new int [3*no_momenta];
    for (int i=0; i<3*no_momenta; i++)
    {
      momenta_x[i] = (i%3 ? 0 : i/3);
      momenta_y[i] = ((i%3==1) ? i/3 : 0);
      momenta_z[i] = ((i%3==2) ? i/3 : 0);
    }
    momenta = momenta_x;
    stringstream source_pos;
    string prefix;
    source_pos.str("");
    if ((mode==blockxyz)||(mode==blockxyzsingle)) 
    {
      source_pos << "sx" << x_source[0] << "sy" << x_source[1]<< "sz" << x_source[2]<< "st" << x_source[3];
      prefix = "block_x_";
    }
    else
    {
      prefix = "p_x";
    }
    if (x_corr!=NULL) write_file(prefix + source_pos.str() + "." + stream.str(), x_corr); // still need the source position for further use of the data
    momenta = momenta_y;
    if ((mode==blockxyz)||(mode==blockxyzsingle)) prefix = "block_y_"; else prefix = "p_y";
    if (y_corr!=NULL) write_file(prefix + source_pos.str() + "." + stream.str(), y_corr);
    momenta = momenta_z;
    if ((mode==blockxyz)||(mode==blockxyzsingle)) prefix = "block_z_"; else prefix = "p_z";
    if (z_corr!=NULL) write_file(prefix + source_pos.str() + "." + stream.str(), z_corr);
    momenta = temp_momenta; // restore original (real) momentum list for possible further applications
    no_momenta = temp_no_momenta;
    delete [] momenta_x; // clean up
    delete [] momenta_y;
    delete [] momenta_z;
  }
  cout << "done"; flush(cout);
  return 0;
}


// constructor for class xcorr_2pt
// STATUS: done
xcorr_2pt::xcorr_2pt(string xfname, string pfname, string momentalist_fname) : xcorr(xfname, pfname, momentalist_fname)
{
  if (parse_header()!=0) // parse the file header; requires virtual method (different for 2pt and 3pt functions)
  {
    cout << endl << "Error parsing header of correlator file: '" << xfname << "'.";
    error = 3;
  }
  radius = L;    // sum over entire spatial volume by default
  return;
}


// destructor for class xcorr_2pt
// STATUS: done
xcorr_2pt::~xcorr_2pt()
{
  cout << endl << "Cleaning up ... ";
  for (int i=0; i<3; i++)
  {
    if (tags[i]!=NULL) delete [] tags[i];
  }
  return;
}


// corr file header parser for class xcorr_2pt
// STATUS:
int xcorr_2pt::parse_header()
{
  if (xcorr::parse_header()!=0) return -1;
  stringstream stream;
  char dummy;

  stream.str((*reader).parse_header("number of timeslices = "));
  if (stream.str()!="")
  {
    stream >> no_t;
  }
  else
  {
    cout << endl << "FAIL: could not determine number of t-values from header in correlator file";
    return -1;
  }

  if (no_t>0)
  {
    if (t_values!=NULL) delete [] t_values;
    t_values = new int [no_t]; // allocate memory for t-value list
  }
  else
  {
    cout << endl << "FAIL: number of t-values has to be larger than zero";
    return -1;
  }

  stream.str((*reader).parse_header("timeslices = "));
  if (stream.str()!="")
  {
    for (int t=0; t<no_t; t++)
    {
      stream >> t_values[t];
      stream >> dummy;
    }
  }
  else
  {
    cout << endl << "FAIL: could not determine list of t-values from header in correlator file";
    return -1;
  }

  for (int i=0; i<3; i++) // read the information about the three open indices
  {
    stream.str("");
    stream.clear();
    stream << "index" << i << " size = ";
    stream.str((*reader).parse_header(stream.str()));
    if (stream.str()!="")
    {
      stream >> index_sizes[i];
    }
    else
    {
      cout << endl << "FAIL: could not determine size of index " << i << " from header in correlator file";
      return -1;
    }
    stream.str("");
    stream.clear();
    stream << "index" << i << " tags = ";
    stream.str((*reader).parse_header(stream.str()));
    if (index_sizes[i]>0)
    {
      tags[i] = new string [index_sizes[i]];
      for (int j=0; j<index_sizes[i]; j++)
      {
        stream >> tags[i][j];
        tags[i][j].erase(tags[i][j].length()-1, 1); // remove the separator char
      }
    }
    else
    {
      tags[i] = NULL;
    }
  }
  cout << "done"; flush(cout);
  return 0;
}


// calculates all three block corrs with unsummed x,y or z component, while summing over the two other orthogonal components, i.e. projecting them to zero. 
// Sets mode=blockxyz.
// 2pt function version
// STATUS: done
void xcorr_2pt::block_xyz()
{
  cout << endl << "Running FT to generate block corrs for x,y,z direction (orthogonal components projected to zero momentum) ... "; flush(cout);
  double *data = (double *)(* reader).data;
  // clean up in any case:
  if (corr!=NULL) delete [] corr;
  corr=NULL; // here we do not calculate the standard 2pt functions
  if (x_corr!=NULL) delete [] x_corr;
  if (y_corr!=NULL) delete [] y_corr;
  if (z_corr!=NULL) delete [] z_corr; 
  x_corr = new complex<double> [get_L() * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]]; // here we only need L instead of no_momenta
  y_corr = new complex<double> [get_L() * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
  z_corr = new complex<double> [get_L() * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
  for (int t=0; t<(get_L() * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]); t++) // init to zero
  {
    x_corr[t] = complex<double>(0.,0.);
    y_corr[t] = complex<double>(0.,0.);
    z_corr[t] = complex<double>(0.,0.);
  }
  for (int t=0; t<no_t; t++)
  {
    for (int z=0; z<L; z++)
    {
      for (int y=0; y<L; y++)
      {
        for (int x=0; x<L; x++)
        {
          for (int chan=0; chan<index_sizes[0]; chan++)
          {
            for (int isosp=0; isosp<index_sizes[1]; isosp++)
            {
              for (int mu=0; mu<index_sizes[2]; mu++)
              {
                const int64_t ix = 2*(mu + index_sizes[2]*(isosp + index_sizes[1]*(chan + index_sizes[0]*(x + L*(y + L*(z + L*t))))));
                const int ivx = mu + index_sizes[2] * (chan + index_sizes[0] * (x + L * (t + no_t *isosp)));
                const int ivy = mu + index_sizes[2] * (chan + index_sizes[0] * (y + L * (t + no_t *isosp)));
                const int ivz = mu + index_sizes[2] * (chan + index_sizes[0] * (z + L * (t + no_t *isosp)));
                x_corr[ivx] += complex<double>(data[ix], data[ix+1]);
                y_corr[ivy] += complex<double>(data[ix], data[ix+1]);
                z_corr[ivz] += complex<double>(data[ix], data[ix+1]);
              }
            }
          }
        }
      }
    }
  }
  radius = L;      // pointless but consistent
  mode = blockxyz; // set the mode -- important for possible call to FFT_xyz()
  cout << "done"; flush(cout);
  return;
}


// This version applies a given projector and returns only data in the given isospin and inerpolating operator channel.
// Sets mode==blockxyz_single
// STATUS: done
void xcorr_2pt::block_xyz(const unsigned int proj, const isospin_nucleon isospin, const spin_nucleon spin)
{
  cout << endl << "Running FT to generate block corrs for x,y,z direction (orthogonal components projected to zero momentum) ... "; flush(cout);
  double *data = (double *)(* reader).data;
  // clean up in any case:
  if (corr!=NULL) delete [] corr;
  corr=NULL; // here we do not calculate the standard 2pt functions
  if (x_corr!=NULL) delete [] x_corr;
  if (y_corr!=NULL) delete [] y_corr;
  if (z_corr!=NULL) delete [] z_corr; 
  x_corr = new complex<double> [get_L() * no_t]; // here we only need L instead of no_momenta; results are projected and for definite isospin and spin
  y_corr = new complex<double> [get_L() * no_t];
  z_corr = new complex<double> [get_L() * no_t];
  for (int t=0; t<(get_L() * no_t); t++) // init to zero
  {
    x_corr[t] = complex<double>(0.,0.);
    y_corr[t] = complex<double>(0.,0.);
    z_corr[t] = complex<double>(0.,0.);
  }
  // do the projection and reduction of data
  if (PROJECTOR[proj]==NULL) init_projectors(); // initialize projector if neccessary
  for (int t=0; t<no_t; t++)
  {
    for (int z=0; z<L; z++)
    {
      for (int y=0; y<L; y++)
      {
        for (int x=0; x<L; x++)
        {
          gsl_complex c;    // holds a single element of the spin trace
          gsl_complex temp = gsl_complex_rect(0., 0.); // need to store the current, projected correlator element
          const int64_t ix = 2*(0 + index_sizes[2]*(isospin + index_sizes[1]*(spin + index_sizes[0]*(x + L*(y + L*(z + L*t)))))); // mu==0 here
          for (unsigned int i=0; i<4; i++) // project and do the trace
          {
            gsl_matrix_complex_const_view tmp_matrix = gsl_matrix_complex_const_view_array(&(data[ix]), 4, 4);
            gsl_vector_complex_const_view tmp1 = gsl_matrix_complex_const_row(PROJECTOR[proj], i);
            gsl_vector_complex_const_view tmp2 = gsl_matrix_complex_const_column(&tmp_matrix.matrix, i);
            gsl_blas_zdotu(&tmp1.vector, &tmp2.vector, &c);
            temp = gsl_complex_add(temp, c);
	        }
          x_corr[x + L * t] += complex<double>(GSL_REAL(temp), GSL_IMAG(temp));
          y_corr[y + L * t] += complex<double>(GSL_REAL(temp), GSL_IMAG(temp));
          z_corr[z + L * t] += complex<double>(GSL_REAL(temp), GSL_IMAG(temp));
	      }
      }
    }
  } 
  radius = L;            // pointless but consistent
  mode = blockxyzsingle; // set the mode -- important for possible call to FFT_xyz()
  cout << "done"; flush(cout);
  return;
}


// Calculates the full lattice momentum dependence in x,y,z direction using FFT while projecting the other two "orthotogonal" componentes to zero 
// (e.g. for (p_x, 0, 0) we FFT in "x" and sum over y,z. If block_xyz() was called previously (i.e. of mode===blockxyz is set) it uses the available 
// results, otherwise it calls block_xyz() by itself before running the FFT. Sets mode==FFTxyz.
// 2pt function version
// STATUS: implement
void xcorr_2pt::FFT_xyz()
{
  if (mode!=blockxyz) block_xyz(); // check if the corresponding block_xyz() method has been called before
  cout << endl << "Running FFT in all directions ... "; 
  radius = L;    // pointless but consistent
  mode = FFTxyz; // set the mode -- important for possible call to FFT_xyz()
  cout << "done"; flush(cout);
  return;
}


// This version applies a given projector and returns only data in the given isospin and inerpolating operator channel.
// Sets mode==FTTxyzsingle
// STATUS: done
void xcorr_2pt::FFT_xyz(const unsigned int proj, const isospin_nucleon isospin, const spin_nucleon spin)
{
  if (mode!=blockxyzsingle) block_xyz(proj, isospin, spin); // check if the corresponding block_xyz() method has been called before
  cout << endl << "Running FFT (separately) in p_x, p_y, p_z directions ... ";
  const int64_t size = L*no_t;
  fftw_complex *temp = new fftw_complex [size];
  fftw_plan plan = fftw_create_plan_specific(L, FFTW_FORWARD, FFTW_IN_PLACE|FFTW_MEASURE|FFTW_USE_WISDOM, temp, 1, NULL, 1); // 2pt function uses FFTW_FORWARD, i.e. "-1" in exp()
  convert_fftw_complex(x_corr, temp, size); // need to move data to the temp array for the FFTW
  fftw(plan, no_t, temp, 1, L, NULL, 0, 0);
  convert_fftw_complex(temp, x_corr, size);
  convert_fftw_complex(y_corr, temp, size);
  fftw(plan, no_t, temp, 1, L, NULL, 0, 0);
  convert_fftw_complex(temp, y_corr, size);
  convert_fftw_complex(z_corr, temp, size);
  fftw(plan, no_t, temp, 1, L, NULL, 0, 0);
  convert_fftw_complex(temp, z_corr, size);
  fftw_destroy_plan(plan);
  delete [] temp;
  for (int t=0; t<no_t; t++)
  {
    for (int q=0; q<L; q++)
    {
      const double momentum = 2. * pi * (double)q / (double)L;
      x_corr[q+t*L] *= exp(complex<double>(0.,1.) * momentum * (double)x_source[0]); // finally we apply the phase factor for the source position
      y_corr[q+t*L] *= exp(complex<double>(0.,1.) * momentum * (double)x_source[1]); // mind the "+" (because of "-p*(-x_0)" for the 2pt function)
      z_corr[q+t*L] *= exp(complex<double>(0.,1.) * momentum * (double)x_source[2]);
    }
  }
  radius = L;          // pointless but consistent
  mode = FFTxyzsingle; // set the mode -- important for possible call to FFT_xyz()
  cout << "done"; flush(cout);
  return;
}


// calculates Fourier trafo x->q with spatial 'radius' for given momentalist and continuum derivative with respect to q_i 
// (i.e. spatial components, factor "x_i" in the FT). Sets mode=cont and requires mode==cont for iterative calculation.
// 2pt function version
// STATUS: done
void xcorr_2pt::FT_Dq_cont(int r)
{
  cout << endl << "Running FT with continuum derivative for r = " << r << " ... "; flush(cout);
  double *data = (double *)(* reader).data;
  if ((radius==L)||(r<radius)||(mode!=cont)) // delete if 1) no previous FT has been performed (indicated by radius==L)
  {                                          //           2) previous FT run over entire lattice already (-> re-run, indicated by radius==L),
    if (corr!=NULL) delete [] corr;          //           3) previous FT already run over larger region (indicated by r<radius)
    if (x_corr!=NULL) delete [] x_corr;      //           4) an imcompatible FT has been called before (prevent mixing of derivatives in the summation)
    if (y_corr!=NULL) delete [] y_corr;
    if (z_corr!=NULL) delete [] z_corr;
    corr = new complex<double> [no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
    x_corr = new complex<double> [no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
    y_corr = new complex<double> [no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
    z_corr = new complex<double> [no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
    for (int t=0; t<(no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]); t++) // init to zero
    {
      corr[t] = complex<double>(0.,0.);
      x_corr[t] = complex<double>(0.,0.);
      y_corr[t] = complex<double>(0.,0.);
      z_corr[t] = complex<double>(0.,0.);
    }
    radius = 0; // this resets the start values for the loops in (x,y,z)-direction !!
  } // ... in any other case memory is already allocated
  complex<double> ipx(0.,0.);
  complex<double> xcorr(0.,0.);
  for (int i=0; i<no_momenta; i++) // sweep through momenta list
  {
    for (int t=0; t<no_t; t++)
    {
      for (int z=L/2; z>-L/2; z--) // cutoff in x / x-Derivative
      {
        for (int y=L/2; y>-L/2; y--)
        {
          for (int x=r/2; x>-r/2; x--)
          {
            if (((x>0)&&(x>radius/2)) || ((x<=0)&&(x<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int chan=0; chan<index_sizes[0]; chan++)
              {
                for (int isosp=0; isosp<index_sizes[1]; isosp++)
                {
                  for (int mu=0; mu<index_sizes[2]; mu++)
                  {
                    // we assume that x,y,z are for source at zero and only shift the array indexing itself
                    // this allows to keep simple for-loops
                    const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                    const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                    const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                    const int64_t ix = 2*(mu + index_sizes[2]*(isosp + index_sizes[1]*(chan + index_sizes[0]*(X + L*(Y + L*(Z + L*t))))));
                    ipx = complex<double>(0., - momenta[3*i]*x - momenta[3*i+1]*y - momenta[3*i+2]*z); // -i*x*p => mind the "-" for the 2pt function!!
                    ipx /= (double)L;
                    ipx *= 2. * pi;
                    ipx = exp(ipx);                    
                    xcorr = complex<double>(data[ix], data[ix+1]);
                    const int iv = mu + index_sizes[2] * (chan + index_sizes[0] * (i + no_momenta * (t + no_t *isosp)));
                    corr[iv] += ipx * xcorr; // here we also sum the normal one -- note that the cutoff is not yet done in all directions for this one. This will be fixed / no longer required anyways once we write out the x,y,z single block corrs with other components projected to zero
                    x_corr[iv] += complex<double>(0., -1.) * (double)x * ipx * xcorr; // mind the minus for the 2pt function
                  }
                }
              }
            }
          }
        }
      }
      for (int z=L/2; z>-L/2; z--) // cutoff in y / y-Derivative
      {
        for (int y=r/2; y>-r/2; y--)
        {
          for (int x=L/2; x>-L/2; x--)
          {
            if (((y>0)&&(y>radius/2)) || ((y<=0)&&(y<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int chan=0; chan<index_sizes[0]; chan++)
              {
                for (int isosp=0; isosp<index_sizes[1]; isosp++)
                {
                  for (int mu=0; mu<index_sizes[2]; mu++)
                  {
                    // we assume that x,y,z are for source at zero and only shift the array indexing itself
                    // this allows to keep simple for-loops
                    const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                    const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                    const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                    const int64_t ix = 2*(mu + index_sizes[2]*(isosp + index_sizes[1]*(chan + index_sizes[0]*(X + L*(Y + L*(Z + L*t))))));
                    ipx = complex<double>(0., - momenta[3*i]*x - momenta[3*i+1]*y - momenta[3*i+2]*z); // -i*x*p => mind the "-" for the 2pt function!!
                    ipx /= (double)L;
                    ipx *= 2. * pi;
                    ipx = exp(ipx);                    
                    xcorr = complex<double>(data[ix], data[ix+1]);
                    const int iv = mu + index_sizes[2] * (chan + index_sizes[0] * (i + no_momenta * (t + no_t *isosp)));
                    y_corr[iv] += complex<double>(0., -1.) * (double)y * ipx * xcorr; // mind the minus for the 2pt function
                  }
                }
              }
            }
          }
        }
      }
      for (int z=r/2; z>-r/2; z--) // cutoff in z / z-Derivative
      {
        for (int y=L/2; y>-L/2; y--)
        {
          for (int x=L/2; x>-L/2; x--)
          {
            if (((z>0)&&(z>radius/2)) || ((z<=0)&&(z<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int chan=0; chan<index_sizes[0]; chan++)
              {
                for (int isosp=0; isosp<index_sizes[1]; isosp++)
                {
                  for (int mu=0; mu<index_sizes[2]; mu++)
                  {
                    // we assume that x,y,z are for source at zero and only shift the array indexing itself
                    // this allows to keep simple for-loops
                    const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                    const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                    const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                    const int64_t ix = 2*(mu + index_sizes[2]*(isosp + index_sizes[1]*(chan + index_sizes[0]*(X + L*(Y + L*(Z + L*t))))));
                    ipx = complex<double>(0., - momenta[3*i]*x - momenta[3*i+1]*y - momenta[3*i+2]*z); // -i*x*p => mind the "-" for the 2pt function!!
                    ipx /= (double)L;
                    ipx *= 2. * pi;
                    ipx = exp(ipx);                    
                    xcorr = complex<double>(data[ix], data[ix+1]);
                    const int iv = mu + index_sizes[2] * (chan + index_sizes[0] * (i + no_momenta * (t + no_t *isosp)));
                    z_corr[iv] += complex<double>(0., -1.) * (double)z * ipx * xcorr; // mind the minus for the 2pt function
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  radius=r;
  mode = cont; // set the current mode to continuum-like derivative
  cout << "done"; flush(cout);
  return;
}


// calculates Fourier trafo x->q with spatial 'radius' for given momentalist and lattice derivative with respect to q_i
// (i.e. spatial components, factor "sin(x_i)" in the FT). Sets mode=latt and requires mode==latt for iterative calculation.
// 2pt function version
// STATUS: done
void xcorr_2pt::FT_Dq_latt(int r)
{
  cout << endl << "Running FT with lattice derivative for r = " << r << " ... "; flush(cout);
  double *data = (double *)(* reader).data;
  if ((radius==L)||(r<radius)||(mode!=latt)) // delete if 1) no previous FT has been performed (indicated by radius==L)
  {                                          //           2) previous FT run over entire lattice already (-> re-run, indicated by radius==L),
    if (corr!=NULL) delete [] corr;          //           3) previous FT already run over larger region (indicated by r<radius)
    if (x_corr!=NULL) delete [] x_corr;      //           4) an imcompatible FT has been called before (prevent mixing of derivatives in the summation)
    if (y_corr!=NULL) delete [] y_corr;
    if (z_corr!=NULL) delete [] z_corr;
    corr = new complex<double> [no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
    x_corr = new complex<double> [no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
    y_corr = new complex<double> [no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
    z_corr = new complex<double> [no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]];
    for (int t=0; t<(no_momenta * no_t * index_sizes[0] * index_sizes[1] * index_sizes[2]); t++) // init to zero
    {
      corr[t] = complex<double>(0.,0.);
      x_corr[t] = complex<double>(0.,0.);
      y_corr[t] = complex<double>(0.,0.);
      z_corr[t] = complex<double>(0.,0.);
    }
    radius = 0; // this resets the start values for the loops in (x,y,z)-direction !!
  } // ... in any other case memory is already allocated
  complex<double> ipx(0.,0.);
  complex<double> xcorr(0.,0.);
  for (int i=0; i<no_momenta; i++) // sweep through momenta list
  {
    for (int t=0; t<no_t; t++)
    {
      for (int z=L/2; z>-L/2; z--) // cutoff in x / x-Derivative
      {
        for (int y=L/2; y>-L/2; y--)
        {
          for (int x=r/2; x>-r/2; x--)
          {
            if (((x>0)&&(x>radius/2)) || ((x<=0)&&(x<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int chan=0; chan<index_sizes[0]; chan++)
              {
                for (int isosp=0; isosp<index_sizes[1]; isosp++)
                {
                  for (int mu=0; mu<index_sizes[2]; mu++)
                  {
                    // we assume that x,y,z are for source at zero and only shift the array indexing itself
                    // this allows to keep simple for-loops
                    const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                    const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                    const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                    const int64_t ix = 2*(mu + index_sizes[2]*(isosp + index_sizes[1]*(chan + index_sizes[0]*(X + L*(Y + L*(Z + L*t))))));
                    complex<double> Dx(0.,(double)L/(2.*(double)pi)); // additional factors for the derivatives; initialize with the additional factors, from the sin() ("2i") and the lattice derivative itself, i.e. L/(2*pi) resulting in i*L/pi
                    Dx *= sin(-(double)x*2.*pi/(double)L);            // mind the minus for the 2pt function
                    ipx = complex<double>(0., - momenta[3*i]*x - momenta[3*i+1]*y - momenta[3*i+2]*z); // -i*x*p => mind the "-" for the 2pt function!!
                    ipx /= (double)L;
                    ipx *= 2. * pi;
                    ipx = exp(ipx);                    
                    xcorr = complex<double>(data[ix], data[ix+1]);
                    const int iv = mu + index_sizes[2] * (chan + index_sizes[0] * (i + no_momenta * (t + no_t *isosp)));
                    corr[iv] += ipx * xcorr; // here we also sum the normal one -- note that the cutoff is not yet done in all directions for this one. This will be fixed / no longer required anyways once we write out the x,y,z single block corrs with other components projected to zero
                    x_corr[iv] += Dx * ipx * xcorr;
                  }
                }
              }
            }
          }
        }
      }
      for (int z=L/2; z>-L/2; z--) // cutoff in y / y-Derivative
      {
        for (int y=r/2; y>-r/2; y--)
        {
          for (int x=L/2; x>-L/2; x--)
          {
            if (((y>0)&&(y>radius/2)) || ((y<=0)&&(y<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int chan=0; chan<index_sizes[0]; chan++)
              {
                for (int isosp=0; isosp<index_sizes[1]; isosp++)
                {
                  for (int mu=0; mu<index_sizes[2]; mu++)
                  {
                    // we assume that x,y,z are for source at zero and only shift the array indexing itself
                    // this allows to keep simple for-loops
                    const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                    const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                    const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                    const int64_t ix = 2*(mu + index_sizes[2]*(isosp + index_sizes[1]*(chan + index_sizes[0]*(X + L*(Y + L*(Z + L*t))))));
                    complex<double> Dy(0.,(double)L/(2.*(double)pi));
                    Dy *= sin(-(double)y*2.*pi/(double)L);
                    ipx = complex<double>(0., - momenta[3*i]*x - momenta[3*i+1]*y - momenta[3*i+2]*z); // -i*x*p => mind the "-" for the 2pt function!!
                    ipx /= (double)L;
                    ipx *= 2. * pi;
                    ipx = exp(ipx);                    
                    xcorr = complex<double>(data[ix], data[ix+1]);
                    const int iv = mu + index_sizes[2] * (chan + index_sizes[0] * (i + no_momenta * (t + no_t *isosp)));
                    y_corr[iv] += Dy * ipx * xcorr;
                  }
                }
              }
            }
          }
        }
      }
      for (int z=r/2; z>-r/2; z--) // cutoff in z / z-Derivative
      {
        for (int y=L/2; y>-L/2; y--)
        {
          for (int x=L/2; x>-L/2; x--)
          {
            if (((z>0)&&(z>radius/2)) || ((z<=0)&&(z<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int chan=0; chan<index_sizes[0]; chan++)
              {
                for (int isosp=0; isosp<index_sizes[1]; isosp++)
                {
                  for (int mu=0; mu<index_sizes[2]; mu++)
                  {
                    // we assume that x,y,z are for source at zero and only shift the array indexing itself
                    // this allows to keep simple for-loops
                    const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                    const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                    const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                    const int64_t ix = 2*(mu + index_sizes[2]*(isosp + index_sizes[1]*(chan + index_sizes[0]*(X + L*(Y + L*(Z + L*t))))));
                    complex<double> Dz(0.,(double)L/(2.*(double)pi)); 
                    Dz *= sin(-(double)z*2.*pi/(double)L);
                    ipx = complex<double>(0., - momenta[3*i]*x - momenta[3*i+1]*y - momenta[3*i+2]*z); // -i*x*p => mind the "-" for the 2pt function!!
                    ipx /= (double)L;
                    ipx *= 2. * pi;
                    ipx = exp(ipx);                    
                    xcorr = complex<double>(data[ix], data[ix+1]);
                    const int iv = mu + index_sizes[2] * (chan + index_sizes[0] * (i + no_momenta * (t + no_t *isosp)));
                    z_corr[iv] += Dz * ipx * xcorr;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  radius=r;
  mode = latt; // set the current mode to continuum-like derivative
  cout << "done"; flush(cout);
  return;
}


// writes the correlator in momentum space to the file given by pfilename, if pfilename="" write to stdout (2pt function version)
// STATUS: done
int xcorr_2pt::write_file(string fname, complex<double> *data)
{
  ostream *stream;
  ofstream ofs;
  if (fname!="")
  {
    ofs.open(fname.c_str());
    if (!(ofs.is_open()))
    {
      cout << "FAIL: Could not open file '" << fname << "' for writing data.";
      return -1;
    }
    stream = &ofs;
  }
  else
  {
    stream = &cout;
  }

  if ((mode!=blockxyzsingle)&&(mode!=FFTxyzsingle)) // either write unprojected 2pt function with all indices open ...
  {
    for (int isosp=0; isosp<index_sizes[1]; isosp++)
    {
      for (int t=0; t<no_t; t++)
      {
        for (int i=0; i<no_momenta; i++)
        {
          for (int chan=0; chan<index_sizes[0]; chan++)
          {
            for (int mu1=0; mu1<sqrt(index_sizes[2]); mu1++) // create n x n (4 x 4) matrix
            {
              (* stream).fill(' ');
              (* stream).width(4);
              (* stream) << t_values[t] << " " << showpos << momenta[3*i] << " " << momenta[3*i+1] << " " << momenta[3*i+2]; // write t and momenta
              for (int mu2=0; mu2<sqrt(index_sizes[2]); mu2++)
              {
                int iv = mu2 + sqrt(index_sizes[2]) * (mu1 + sqrt(index_sizes[2]) * (chan + index_sizes[0] * (i + no_momenta *(t + no_t *isosp))));
                (* stream) << scientific << showpos << " " << real(data[iv]) << " " << imag(data[iv]) << " " << noshowpos; // write Re and Im of corr data
              }
              (* stream).unsetf(ios::fixed | ios::scientific); // reset to non-scientific output
              (* stream) << tags[1][isosp] << " " << tags[0][chan] << endl; // add the tags for isospin and the different channels for gevp
            }
          }
        }
      }
    }
  }
  else // or write a projected 2pt function with all indices contracted
  {
    for (int t=0; t<no_t; t++) 
    {
      for (int i=0; i<no_momenta; i++)
      {
        (* stream).fill(' ');
        (* stream).width(4);
        (* stream) << t_values[t] << " " << showpos << momenta[3*i] << " " << momenta[3*i+1] << " " << momenta[3*i+2]; // write t and momenta
        (* stream) << scientific << showpos << " " << real(data[i + no_momenta * t]) << " " << imag(data[i + no_momenta * t]) << noshowpos << endl; // write Re and Im of corr data
        (* stream).unsetf(ios::fixed | ios::scientific); // reset to non-scientific output
      }
    }
  }
  if (fname!="") ofs.close();
  return 0;
}


// dumps a snapshot of the current results for corr and the derivative corr to disk
// STATUS: done
int xcorr_2pt::write_results()
{
  return xcorr::write_results("2pt");
}


// constructor for class xcorr_3pt
// STATUS: done
xcorr_3pt::xcorr_3pt(string xfname, string pfname, string momentalist_fname) : xcorr(xfname, pfname, momentalist_fname)
{
  if (parse_header()!=0) // parse the file header; requires virtual method (different for 2pt and 3pt functions)
  {
    cout << endl << "Error parsing header of correlator file: '" << xfname << "'.";
    error = 3;
  }
  radius = L; // sum over entire spatial volume by default
  return;
}


// destructor for class xcorr_3pt
// STATUS: done
xcorr_3pt::~xcorr_3pt()
{
  cout << endl << "Cleaning up ... "; flush(cout);
  return;
}


// corr file header parser for class xcorr_3pt
// STATUS: done
int xcorr_3pt::parse_header()
{
  if (xcorr::parse_header()!=0) return -1; // call inherited part of the parser first
  char dummy;
  stringstream stream;
  stream.str((*reader).parse_header("number of insertions = "));
  if (stream.str()!="")
  {
    stream >> no_t;
  }
  else
  {
    cout << endl << "FAIL: could not determine number of t-values from header in correlator file";
    return -1;
  }

  if (no_t>0)
  {
    if (t_values!=NULL) delete [] t_values;
    t_values = new int [no_t]; // allocate memory for t-value list
  }
  else
  {
    cout << endl << "FAIL: number of t-values has to be larger than zero";
    return -1;
  }

  stream.str((*reader).parse_header("t-insertions = "));
  if (stream.str()!="")
  {
    for (int t=0; t<no_t; t++)
    {
      stream >> t_values[t];
      stream >> dummy;
    }
  }
  else
  {
    cout << endl << "FAIL: could not determine list of t-values from header in correlator file";
    return -1;
  }

  stream.str((*reader).parse_header("lattice-site size = "));
  if (stream.str()!="")
  {
    stream >> site_size;
  }
  else
  {
    cout << endl << "FAIL: could not determine site size from header in correlator file";
    return -1;
  }
  cout << "done"; flush(cout);
  return 0;
}


// calculates all three block corrs with unsummed x,y or z component, while summing over the two other orthogonal components, i.e. projecting them to zero. 
// Sets mode=block_xyz.
// 3pt function version
// STATUS: done
void xcorr_3pt::block_xyz()
{
  cout << endl << "Running FT to generate block corrs for x,y,z direction (orthogonal components projected to zero momentum) ... "; flush(cout);
  double *data = (double *)(* reader).data;
  // clean up in any case:
  if (corr!=NULL) delete [] corr;
  corr=NULL; // here we do not calculate the standard 2pt functions
  if (x_corr!=NULL) delete [] x_corr;
  if (y_corr!=NULL) delete [] y_corr;
  if (z_corr!=NULL) delete [] z_corr; 
  x_corr = new complex<double> [get_L() * no_t * site_size]; // here we only need L instead of no_momenta
  y_corr = new complex<double> [get_L() * no_t * site_size];
  z_corr = new complex<double> [get_L() * no_t * site_size];
  for (int t=0; t<(get_L() * no_t * site_size); t++) // init to zero
  {
    x_corr[t] = complex<double>(0.,0.);
    y_corr[t] = complex<double>(0.,0.);
    z_corr[t] = complex<double>(0.,0.);
  }
  for (int t=0; t<no_t; t++)
  {
    for (int z=0; z<L; z++)
    {
      for (int y=0; y<L; y++)
      {
        for (int x=0; x<L; x++)
        {
          for (int s=0; s<site_size; s++)
          {
            const int64_t ix =2*(s + site_size*(x + L*(y + L*(z + L*t))));
            x_corr[(x + t * L)* site_size + s] += complex<double>(data[ix], data[ix+1]);
            y_corr[(y + t * L)* site_size + s] += complex<double>(data[ix], data[ix+1]);
            z_corr[(z + t * L)* site_size + s] += complex<double>(data[ix], data[ix+1]);
          }
        }
      }
    }
  }
  radius = L;      // pointless but consistent
  mode = blockxyz; // set the mode -- important for possible call to FFT_xyz()
  cout << "done"; flush(cout);
  return;
}


// Dummy to catch virtual function calls
// STATUS: done
void xcorr_3pt::block_xyz(const unsigned int, const isospin_nucleon, const spin_nucleon)
{
  cout << endl << "WARNING: pointless (possibly virtual?) call to xcorr_3pt::block_xyz(const unsigned int proj, const unsigned int isospin, const unsigned int spin). Calling standard xcorr_3pt::block_xyz() instead...";
  block_xyz();
  return;
}


// calculates the full lattice momentum dependence in x,y,z direction using FFT while projecting the other two "orthotogonal" componentes to zero 
// (e.g. for (p_x, 0, 0) we FFT in "x" and sum over y,z. If block_xyz() was called previously (i.e. of mode===blockxyz is set) it uses the available 
// results, otherwise it calls block_xyz() by itself before running the FFT. Sets mode==FFTxyz.
// 3pt function version
// STATUS: implement
void xcorr_3pt::FFT_xyz()
{
  if (mode!=blockxyz) block_xyz(); // check if the corresponding block_xyz() method has been called before
  cout << endl << "Running FFT (separately) in p_x, p_y, p_z directions ... ";
  const int64_t size = L*no_t*site_size;
  fftw_complex *temp = new fftw_complex [size];
  fftw_plan plan = fftw_create_plan_specific(L, FFTW_BACKWARD, FFTW_IN_PLACE|FFTW_MEASURE|FFTW_USE_WISDOM, temp, site_size, NULL, 0); // 3pt function uses FFTW_BACKWARD, i.e. "+1" in exp()
  convert_fftw_complex(x_corr, temp, size); // need to move data to the temp array for the FFTW
  for (int site=0; site<site_size; site++) fftw(plan, no_t, &(temp[site]), site_size, L*site_size, NULL, 0, 0);
  convert_fftw_complex(temp, x_corr, size);
  convert_fftw_complex(y_corr, temp, size);
  for (int site=0; site<site_size; site++) fftw(plan, no_t, &(temp[site]), site_size, L*site_size, NULL, 0, 0);
  convert_fftw_complex(temp, y_corr, size);
  convert_fftw_complex(z_corr, temp, size);
  for (int site=0; site<site_size; site++) fftw(plan, no_t, &(temp[site]), site_size, L*site_size, NULL, 0, 0);
  convert_fftw_complex(temp, z_corr, size);
  fftw_destroy_plan(plan);
  delete [] temp;
  for (int t=0; t<no_t; t++)
  {
    for (int q=0; q<L; q++)
    {
      for (int site=0; site<site_size; site++)
      {
        const int64_t iv = site + site_size * (q + t * L);
        const double momentum = 2. * pi * (double)q / (double)L;
        x_corr[iv] *= exp(complex<double>(0.,-1.) * momentum * ((double)x_source[0])); // finally we apply the phase factor for the source position
        y_corr[iv] *= exp(complex<double>(0.,-1.) * momentum * ((double)x_source[1])); // mind the "-" (because of "+p*(-x_0)" for the 2pt function)
        z_corr[iv] *= exp(complex<double>(0.,-1.) * momentum * ((double)x_source[2]));
      }
    }
  }
  radius = L;    // pointless but consistent
  mode = FFTxyz; // set the mode -- important for possible call to FFT_xyz()
  cout << "done"; flush(cout);
  return;
}


// Dummy to catch virtual function calls
// STATUS: done
void xcorr_3pt::FFT_xyz(const unsigned int, const isospin_nucleon, const spin_nucleon)
{
  cout << endl << "WARNING: pointless (possibly virtual?) call to xcorr_3pt::FFT_xyz(const unsigned int proj, const unsigned int isospin, const unsigned int spin). Calling standard xcorr_3pt::FFT_xyz() instead...";
  FFT_xyz();
  return;
}


// calculates Fourier trafo x->q with spatial 'radius' for given momentalist and continuum derivative with respect to q_i
// (i.e. spatial components, factor "x_i" in the FT). Sets mode=cont and requires mode==cont for iterative calculation.
// 3pt function version
// STATUS: done
void xcorr_3pt::FT_Dq_cont(int r)
{
  cout << endl << "Running FT for r = " << r << " ... "; flush(cout);
  double *data = (double *)(* reader).data;
  if ((radius==L)||(r<radius)||(mode!=cont)) // delete if 1) no previous FT has been performed (indicated by radius==L)
  {                                          //           2) previous FT run over entire lattice already (-> re-run, indicated by radius==L),
    if (corr!=NULL) delete [] corr;          //           3) previous FT already run over larger region (indicated by r<radius)
    if (x_corr!=NULL) delete [] x_corr;      //           4) an imcompatible FT has been called before (prevent mixing of derivatives in the summation)
    if (y_corr!=NULL) delete [] y_corr;
    if (z_corr!=NULL) delete [] z_corr;
    corr = new complex<double> [no_momenta * no_t * site_size];
    x_corr = new complex<double> [no_momenta * no_t * site_size];
    y_corr = new complex<double> [no_momenta * no_t * site_size];
    z_corr = new complex<double> [no_momenta * no_t * site_size];
    for (int t=0; t<(no_momenta * no_t * site_size); t++) // init to zero
    {
      corr[t] = complex<double>(0.,0.);
      x_corr[t] = complex<double>(0.,0.);
      y_corr[t] = complex<double>(0.,0.);
      z_corr[t] = complex<double>(0.,0.);
    }
    radius = 0; // this resets the start values for the loops in (x,y,z)-direction !!
  } // ... in any other case memory is already allocated
  complex<double> ipx(0.,0.);
  complex<double> xcorr(0.,0.); 
  for (int i=0; i<no_momenta; i++) // sweep through momenta list
  {
    for (int t=0; t<no_t; t++)
    {
      for (int z=L/2; z>-L/2; z--)  // cutoff in x / x-Derivative
      {
        for (int y=L/2; y>-L/2; y--)
        {
          for (int x=r/2; x>-r/2; x--)
          {
            if (((x>0)&&(x>radius/2)) || ((x<=0)&&(x<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int s=0; s<site_size; s++)
              {
                // we assume that x,y,z are for source at zero and only shift the array indexing itself
                // this allows to keep simple for-loops
                const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                const int64_t ix = 2*(s + site_size*(X + L*(Y + L*(Z + L*t))));
                ipx = complex<double>(0., +momenta[3*i]*x + momenta[3*i+1]*y + momenta[3*i+2]*z); // +i*x*p => mind the "+" for the 3pt function!!
                ipx /= (double)L;
                ipx *= 2. * pi;
                ipx = exp(ipx);
                xcorr = complex<double>(data[ix], data[ix+1]);
                corr[(i + t * no_momenta)* site_size + s] += ipx * xcorr;
                x_corr[(i + t * no_momenta)* site_size + s] += complex<double>(0., +1.) * (double)x * ipx * xcorr; //mind the plus for the 3pt function
              }
            }
          }
        }
      }
      for (int z=L/2; z>-L/2; z--)  // cutoff in y / y-Derivative
      {
        for (int y=r/2; y>-r/2; y--)
        {
          for (int x=L/2; x>-L/2; x--)
          {
            if (((y>0)&&(y>radius/2)) || ((y<=0)&&(y<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int s=0; s<site_size; s++)
              {
                // we assume that x,y,z are for source at zero and only shift the array indexing itself
                // this allows to keep simple for-loops
                const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                const int64_t ix = 2*(s + site_size*(X + L*(Y + L*(Z + L*t))));
                ipx = complex<double>(0., +momenta[3*i]*x + momenta[3*i+1]*y + momenta[3*i+2]*z); // +i*x*p => mind the "+" for the 3pt function!!
                ipx /= (double)L;
                ipx *= 2. * pi;
                ipx = exp(ipx);
                xcorr = complex<double>(data[ix], data[ix+1]);
                y_corr[(i + t * no_momenta)* site_size + s] += complex<double>(0., +1.) * (double)y * ipx * xcorr; // mind the plus for the 3pt function 
              }
            }
          }
        }
      }
      for (int z=r/2; z>-r/2; z--)  // cutoff in z / z-Derivative
      {
        for (int y=L/2; y>-L/2; y--)
        {
          for (int x=L/2; x>-L/2; x--)
          {
            if (((z>0)&&(z>radius/2)) || ((z<=0)&&(z<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int s=0; s<site_size; s++)
              {
                // we assume that x,y,z are for source at zero and only shift the array indexing itself
                // this allows to keep simple for-loops
                const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                const int64_t ix = 2*(s + site_size*(X + L*(Y + L*(Z + L*t))));
                ipx = complex<double>(0., +momenta[3*i]*x + momenta[3*i+1]*y + momenta[3*i+2]*z); // +i*x*p => mind the "+" for the 3pt function!!
                ipx /= (double)L;
                ipx *= 2. * pi;
                ipx = exp(ipx);
                xcorr = complex<double>(data[ix], data[ix+1]);
                z_corr[(i + t * no_momenta)* site_size + s] += complex<double>(0., +1.) * (double)z * ipx * xcorr; // mind the plus for the 3pt function 
              }
            }
          }
        }
      }
    }
  }
  radius=r;
  mode = cont;
  cout << "done"; flush(cout);
  return;
}


// calculates Fourier trafo x->q with spatial 'radius' for given momentalist and lattice derivative with respect to q_i
// (i.e. spatial components, factor "sin(x_i)" in the FT). Sets mode=latt and requires mode==latt for iterative calculation.
// 3pt function version
// STATUS: done
void xcorr_3pt::FT_Dq_latt(int r)
{
  cout << endl << "Running FT for r = " << r << " ... "; flush(cout);
  double *data = (double *)(* reader).data;
  if ((radius==L)||(r<radius)||(mode!=latt)) // delete if 1) no previous FT has been performed (indicated by radius==L)
  {                                          //           2) previous FT run over entire lattice already (-> re-run, indicated by radius==L),
    if (corr!=NULL) delete [] corr;          //           3) previous FT already run over larger region (indicated by r<radius)
    if (x_corr!=NULL) delete [] x_corr;      //           4) an imcompatible FT has been called before (prevent mixing of derivatives in the summation)
    if (y_corr!=NULL) delete [] y_corr;
    if (z_corr!=NULL) delete [] z_corr;
    corr = new complex<double> [no_momenta * no_t * site_size];
    x_corr = new complex<double> [no_momenta * no_t * site_size];
    y_corr = new complex<double> [no_momenta * no_t * site_size];
    z_corr = new complex<double> [no_momenta * no_t * site_size];
    for (int t=0; t<(no_momenta * no_t * site_size); t++) // init to zero
    {
      corr[t] = complex<double>(0.,0.);
      x_corr[t] = complex<double>(0.,0.);
      y_corr[t] = complex<double>(0.,0.);
      z_corr[t] = complex<double>(0.,0.);
    }
    radius = 0; // this resets the start values for the loops in (x,y,z)-direction !!
  } // ... in any other case memory is already allocated
  complex<double> ipx(0.,0.);
  complex<double> xcorr(0.,0.); 
  for (int i=0; i<no_momenta; i++) // sweep through momenta list
  {
    for (int t=0; t<no_t; t++)
    {
      for (int z=L/2; z>-L/2; z--)  // cutoff in x / x-Derivative
      {
        for (int y=L/2; y>-L/2; y--)
        {
          for (int x=r/2; x>-r/2; x--)
          {
            if (((x>0)&&(x>radius/2)) || ((x<=0)&&(x<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int s=0; s<site_size; s++)
              {
                // we assume that x,y,z are for source at zero and only shift the array indexing itself
                // this allows to keep simple for-loops
                const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                complex<double> Dx(0.,(double)L/(2.*pi)); // additional factors for the derivatives: initialize with the additional factors, from the sin() ("2i")  and the lattice derivative itself, i.e. L/(2*pi) resulting in i*L/pi
                Dx *= sin(+(double)x*2.*pi/(double)L); // mind the plus sign for the 3pt function
                const int64_t ix = 2*(s + site_size*(X + L*(Y + L*(Z + L*t))));
                ipx = complex<double>(0., +momenta[3*i]*x + momenta[3*i+1]*y + momenta[3*i+2]*z); // +i*x*p => mind the "+" for the 3pt function!!
                ipx /= (double)L;
                ipx *= 2. * pi;
                ipx = exp(ipx);
                xcorr = complex<double>(data[ix], data[ix+1]);
                corr[(i + t * no_momenta)* site_size + s] += ipx * xcorr;
                x_corr[(i + t * no_momenta)* site_size + s] += Dx * ipx * xcorr;
              }
            }
          }
        }
      }
      for (int z=L/2; z>-L/2; z--)  // cutoff in y / y-Derivative
      {
        for (int y=r/2; y>-r/2; y--)
        {
          for (int x=L/2; x>-L/2; x--)
          {
            if (((y>0)&&(y>radius/2)) || ((y<=0)&&(y<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int s=0; s<site_size; s++)
              {
                // we assume that x,y,z are for source at zero and only shift the array indexing itself
                // this allows to keep simple for-loops
                const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                const int64_t ix = 2*(s + site_size*(X + L*(Y + L*(Z + L*t))));
                complex<double> Dy(0.,(double)L/(2.*pi));
                Dy *= sin(+(double)y*2.*pi/(double)L);
                ipx = complex<double>(0., +momenta[3*i]*x + momenta[3*i+1]*y + momenta[3*i+2]*z); // +i*x*p => mind the "+" for the 3pt function!!
                ipx /= (double)L;
                ipx *= 2. * pi;
                ipx = exp(ipx);
                xcorr = complex<double>(data[ix], data[ix+1]);
                y_corr[(i + t * no_momenta)* site_size + s] += Dy * ipx * xcorr;
              }
            }
          }
        }
      }
      for (int z=r/2; z>-r/2; z--)  // cutoff in z / z-Derivative
      {
        for (int y=L/2; y>-L/2; y--)
        {
          for (int x=L/2; x>-L/2; x--)
          {
            if (((z>0)&&(z>radius/2)) || ((z<=0)&&(z<=-radius/2))) // only need to sum if we are outside the previously used FT-region
            {
              for (int s=0; s<site_size; s++)
              {
                // we assume that x,y,z are for source at zero and only shift the array indexing itself
                // this allows to keep simple for-loops
                const int X = ((x>=0 ? x : x+L) + x_source[0])%L; // 1. map x back to [0...L); i.e. same frame as x_source
                const int Y = ((y>=0 ? y : y+L) + x_source[1])%L; // 2. correct for actual source position
                const int Z = ((z>=0 ? z : z+L) + x_source[2])%L; // 3. finally enforce again periodicity in L
                const int64_t ix = 2*(s + site_size*(X + L*(Y + L*(Z + L*t))));
                complex<double> Dz(0.,(double)L/(2.*pi));
                Dz *= sin(+(double)z*2.*pi/(double)L);
                ipx = complex<double>(0., +momenta[3*i]*x + momenta[3*i+1]*y + momenta[3*i+2]*z); // +i*x*p => mind the "+" for the 3pt function!!
                ipx /= (double)L;
                ipx *= 2. * pi;
                ipx = exp(ipx);
                xcorr = complex<double>(data[ix], data[ix+1]);
                z_corr[(i + t * no_momenta)* site_size + s] += Dz * ipx * xcorr;
              }
            }
          }
        }
      }
    }
  }
  radius=r;
  mode = latt;
  cout << "done"; flush(cout);
  return;
}


// writes the correlator in momentum space to the file given by pfilename, if pfilename="" write to stdout (3pt function version)
// STATUS: done
int xcorr_3pt::write_file(string fname, complex<double> *data)
{
  ostream *stream;
  ofstream ofs;
  if (fname!="")
  {
    ofs.open(fname.c_str());
    if (!(ofs.is_open()))
    {
      cout << "FAIL: Could not open file '" << fname << "' for writing data.";
      return -1;
    }
    stream = &ofs;
  }
  else
  {
    stream = &cout;
  }

  for (int t=0; t<no_t; t++)
  {
    for (int s=0; s<site_size; s++)
    {
      for (int i=0; i<no_momenta; i++)
      {
        (* stream) << t_values[t] << " " << showpos << momenta[3*i] << " " << momenta[3*i+1] << " " << momenta[3*i+2]; // write t and momenta
        (* stream) << scientific << " " << real(data[(i + t * no_momenta)* site_size + s]) << " " << imag(data[(i + t * no_momenta)* site_size + s]) << noshowpos; // write Re and Im of corr data
        (* stream).unsetf(ios::fixed | ios::scientific); // reset to non-scientific output
        if (site_size==4) // add mu, nu indices if neccessary
        {
          (* stream) << " " << s; 
        }
        if (site_size==10)
        {
          const unsigned int mu = floor((sqrt(8*s+1)-1)/2);
          const unsigned int nu = s - (mu*mu + mu) / 2;
          (* stream) << " " << mu << " " << nu;
        }
        (* stream) << endl;
      }
    }
  }
  if (fname!="") ofs.close();
  return 0;
}


// dumps a snapshot of the current results for corr and the derivative corr to disk
// STATUS: done
int xcorr_3pt::write_results()
{
  return xcorr::write_results("3pt");
}


// Constructor for base class pcorr
// STATUS: done
pcorr::pcorr(const string fname, const string format)
{
  reader = new mixed_text_reader (fname, 0, 0, format); // initialize file reader, read file
  tmax = reader->get_index_max_value(0);   // determine tmax (may be different from T!). Note that 't' has always to be in the first column!!
  tmin = reader->get_index_min_value(0);
  no_momenta = reader->get_lines() / (tmax-tmin+1); // get no_momenta (still need to divide by additional, n-pt function specific factors)
  monotony = 0; // just to be safe...
  p_stride = 1; // ...to be set in constructors of derived classes
  return;
}


// Returns the error state of the protected file reader (this is the only part that can cause non-code related errors)
// STATUS: done
unsigned int pcorr::get_error_id()
{
  return reader->get_error_id();
}


// Same as above, but returns corresponding error message
// STATUS: done
string pcorr::get_error()
{
  return reader->get_error();
}


// Returns the (momentum list) index corresponding to a given momentum
// STATUS: done
int pcorr::get_momentum_index(const int momentum[])
{
  unsigned int p_index = 0;
  unsigned int found = 0;
  const unsigned int stride = reader->get_no_int(); // column-wise (local) stride
  while (p_index<reader->get_lines()) // find index corresponding to 'p'
  {
    if (reader->index[p_index*stride + 1]==momentum[0]) // momentum is expected in the second index column!
    {
      if (reader->index[p_index*stride + 2]==momentum[1])
      {
        if (reader->index[p_index*stride + 3]==momentum[2])
        {
          found = 1;
          break;
        }
      }
    }
    p_index++;
  }
  p_index /= p_stride;   // divide by the file format (line) stride, to get the momentum list index!
  if (!found) return -1; // return invalid index if the momentum has not been found in the index array at all
  return p_index;
}


// Returns p^2 for the given 'p_index'. Returns -1 if p_index is invalid
// STATUS: done
int pcorr::get_momentum_squared(unsigned int p_index)
{
  int p_squared = -1;
  const unsigned int stride = reader->get_no_int() * p_stride; // calculate the full (global) stride
  const int p[3] = {reader->index[p_index*stride + 1], reader->index[p_index*stride + 2], reader->index[p_index*stride + 3]};
  if (p_index<no_momenta)
  {
    p_squared = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
  }
  return p_squared;
}


// Returns true if the momentum list is ordered w.r.t. increasing absolute momenta
// STATUS: done
unsigned int pcorr::get_monotony()
{
  return monotony;
}


// Searches the momentum list for the first occurrence 'i' of a momentum with |p_i|=|p(p_index)|
// Returns invalid index if the original index is not in the list.
// STATUS: check
int pcorr::get_first_matching_p_index(unsigned int p_index)
{
  const int p_squared = get_momentum_squared(p_index);
  if (p_squared<0) return -1;
  p_index = 0;
  while (get_momentum_squared(p_index)<p_squared) // at least the original index fullfils the requirement, hence we do not need to checl for the end of list
  {
    p_index++;
  }
  return p_index;
}


// Search the momentum list for the next occurrence of a momentum p with |p|=|p[p_index]| starting at 'p_index'+1. 
// If 'monotony'==true we assume that the list is sorted w.r.t. increasing absolute momenta. Returns invalid index if the momentum was not found.
// STATUS: check
int pcorr::get_next_matching_p_index(unsigned int p_index)
{
  const int p_squared = get_momentum_squared(p_index);
  p_index++;
  if (monotony)
  {
    if (p_index>=no_momenta) return no_momenta; // special case for which we reached end of list
  }
  else
  {
    while (p_squared!=get_momentum_squared(p_index))
    {
      p_index++;
      if (p_index==no_momenta) break;
    }
  }
  if (p_squared!=get_momentum_squared(p_index)) return -1; // check if momenta was found
  return p_index;
}


// Returns no_momenta
// STATUS: done
unsigned int pcorr::get_no_momenta()
{
  return no_momenta;
}


// Returns t_min
// STATUS: done
unsigned int pcorr::get_tmin()
{
  return tmin;
}


// Returns t_max
// STATUS: done
unsigned int pcorr::get_tmax()
{
  return tmax;
}


// Destructor for base class pcorr
// STATUS: done
pcorr::~pcorr()
{
  delete reader;
  return;
}


// Constructor for class pcorr_2pt_nucleon
// STATUS: done
pcorr_2pt_nucleon::pcorr_2pt_nucleon(const string fname, const string format): pcorr(fname, format)
{
  p_stride = 4 * no_spin_nucleon;
  no_momenta /= no_isospin_nucleon * no_spin_nucleon * 4; // calculate number of momenta -- mind the '4' due to open spin indices (4x4 complex matrices) !!!
  monotony = 1; // set to true ...
  int p = 0; // zero is the smallest possible squared momentum...
  for (unsigned int i=0; i<no_momenta; i++)
  {
    if (get_momentum_squared(i)<p) // ... and check the momemtum list
    {
      monotony = 0;
      break;
    }
    else
    {
      p=get_momentum_squared(i);
    }
  }
  corr = new gsl_vector_complex *[no_momenta * no_isospin_nucleon * no_spin_nucleon * no_projectors]; // NO '4' here... this is where we store the result after contracting indices
  for (unsigned int i=0; i<(no_momenta * no_isospin_nucleon * no_spin_nucleon * no_projectors); i++)
  {
    corr[i] = NULL;
  }
  return;
}


// Performs the desired projection an returns reference to the result. Checks, if the projection was performed before to reduce computational cost
// STATUS: done
gsl_vector_complex * pcorr_2pt_nucleon::get_corr(const unsigned int proj, const unsigned int isospin, const unsigned int spin, const unsigned int p_index)
{
  const unsigned int corr_index =  p_index + no_momenta * (spin + no_spin_nucleon *(isospin + no_isospin_nucleon * proj));
  if (corr[corr_index]==NULL)
  {
    gsl_complex c;
    const unsigned int t_stride = p_stride * no_momenta;
    const unsigned int isospin_stride = t_stride * (tmax-tmin+1);
    if (PROJECTOR[proj]==NULL)
    {
      init_projectors(); // initialize projector if neccessary
    }
    corr[corr_index] = gsl_vector_complex_alloc(tmax-tmin+1);  // allocate the correlator
    gsl_vector_complex_set_zero(corr[corr_index]);
    for (unsigned int t=0; t<(tmax-tmin+1); t++)
    {
      for (unsigned int i=0; i<4; i++) // project and do the trace
      {
        gsl_matrix_complex_const_view tmp_matrix = gsl_matrix_complex_const_view_array(&(reader->data[8 * (t * t_stride + spin * 4 + p_index * p_stride + isospin * isospin_stride)]), 4, 4);
        gsl_vector_complex_const_view tmp1 = gsl_matrix_complex_const_row(PROJECTOR[proj], i);
        gsl_vector_complex_const_view tmp2 = gsl_matrix_complex_const_column(&tmp_matrix.matrix, i);
        gsl_blas_zdotu(&tmp1.vector, &tmp2.vector, &c);
        gsl_vector_complex_set(corr[corr_index], t, gsl_complex_add(gsl_vector_complex_get(corr[corr_index], t), c));
      }
    }
  }
  return corr[corr_index];
}


// Overload for 3pt case (dummy in pcorr_2pt_nucleon class)
// STATUS: done
gsl_vector_complex *pcorr_2pt_nucleon::get_corr(const unsigned int, const unsigned int)
{
  return NULL;
}


// Returns correlator at |p| corresponding to p_index, averaged over all p with same |p|
// '*result' must be allocated by the caller
// STATUS: done
void pcorr_2pt_nucleon::get_corr_momentum_averaged(gsl_vector_complex *result, const unsigned int proj, const unsigned int isospin, const unsigned int spin, int p_index)
{
  unsigned int count = 0;
  gsl_vector_complex_set_zero(result);
  p_index = get_first_matching_p_index(p_index); // find the first index that has the same |p| as p[p_index]
  do
  {
    gsl_vector_complex_add(result, get_corr(proj, isospin, spin, p_index));
    count++;
    p_index = get_next_matching_p_index(p_index); // find next matching index
  }
  while ((p_index<(int)no_momenta)&&(p_index>0));
  gsl_vector_complex_scale(result, gsl_complex_rect(1./count, 0.)); // average
  return;
}


// Returns correlator at |p| corresponding to p_index, averaged over all p with same |p|
// STATUS: done
void pcorr_2pt_nucleon::get_corr_momentum_averaged(gsl_vector_complex *result, const unsigned int, int)
{
  gsl_vector_complex_free(result);
  result = NULL; // force crash if this method is ever called accidentally
  return;
}


// Returns 0 (2pt function do not have insertions)
// STATUS: done
unsigned int pcorr_2pt_nucleon::get_insertion_lex_index(unsigned int insertion_index[])
{
  if (insertion_index!=NULL) // Warn if this method is ever called with non-NULL argument as this very likely indicates some fatal coding error
  {
    cout << endl << "Warning: called get_insertion_lex_index() from a pcorr_2pt_nucleon instance with non-NULL argument!";
  }
  return 0;
}

  
// Returns 0 (2pt functions are always scalar)
// STATUS: done
unsigned int pcorr_2pt_nucleon::get_rank()
{
  return 0;
}


// Destructor for class pcorr_2pt_nucleon
// STATUS: check
pcorr_2pt_nucleon::~pcorr_2pt_nucleon()
{
  for (unsigned int i=0; i<(no_momenta * no_isospin_nucleon * no_spin_nucleon * no_projectors); i++)
  {
    gsl_vector_complex_free(corr[i]);
  }
  delete [] corr;
  return;
}


// Constructor for class pcorr_npt
// STATUS: check
pcorr_npt_nucleon::pcorr_npt_nucleon(const string fname, const string format, const unsigned int insertion_rank): pcorr(fname, format)
{
  rank = insertion_rank;
  p_stride = 1; // for 3pt functions it should always be 1
  const unsigned int insertion_index_size = (unsigned int)gsl_sf_choose(3+rank, rank);
  no_momenta /= insertion_index_size; // calculate number of momenta, divide by number of independent spin tensor elements
  monotony = 1; // set to true ...
  int p = 0; 
  for (unsigned int i=0; i<no_momenta; i++)
  {
    if (get_momentum_squared(i)<p) // ... and check the momemtum list
    {
      monotony = 0;
      break;
    } 
    else
    { 
      p=get_momentum_squared(i);
    } 
  }
  corr = new gsl_vector_complex *[no_momenta * insertion_index_size];
  for (unsigned int i=0; i<(no_momenta * insertion_index_size); i++)
  {
    corr[i] = NULL;
  }
  return;
}


// Returns a reference to the desired correlation function (only generated on demand and stored for further calls in member 'corr')
// STATUS: check
gsl_vector_complex * pcorr_npt_nucleon::get_corr(const unsigned int insertion_index, const unsigned int p_index)
{
  const unsigned int corr_index = p_index + no_momenta * insertion_index;
  if (corr[corr_index]==NULL)
  {
    const unsigned int insertion_index_stride = no_momenta;
    const unsigned int t_stride = no_momenta * ((unsigned int) gsl_sf_choose(3+rank, rank));
    corr[corr_index] = gsl_vector_complex_alloc(tmax-tmin+1);  // allocate the correlator
    gsl_vector_complex_set_zero(corr[corr_index]);
    for (unsigned int t=0; t<(tmax-tmin+1); t++)
    {
      const unsigned int i=2*(p_index + insertion_index * insertion_index_stride + t * t_stride);
      gsl_vector_complex_set(corr[corr_index], t, gsl_complex_rect(reader->data[i], reader->data[i+1]));
    }
  }
  return corr[corr_index];
}


// Overload for 2pt case (dummy in pcorr_npt_nucleon class)
// STATUS: done
gsl_vector_complex *pcorr_npt_nucleon::get_corr(const unsigned int, const unsigned int, const unsigned int, const unsigned int)
{
  return NULL;
}


// Returns correlator at |p| corresponding to p_index, averaged over all p with same |p|
// '*result' must be allocated by the caller; 'p_index' has to be a valid index, no additional error checking is performed by this method
// STATUS: check
void pcorr_npt_nucleon::get_corr_momentum_averaged(gsl_vector_complex *result, const unsigned int insertion_index, int p_index)
{
  unsigned int count = 0;
  gsl_vector_complex_set_zero(result);
  p_index = get_first_matching_p_index(p_index); // find the first index that has the same |p| as p[p_index]
  do
  {
    gsl_vector_complex_add(result, get_corr(insertion_index, p_index));
    count++;
    p_index = get_next_matching_p_index(p_index); // find next matching index 
  }
  while ((p_index<(int)no_momenta)&&(p_index>0));
  gsl_vector_complex_scale(result, gsl_complex_rect(1./count, 0.)); // average
  return;
}


// Overload for 2pt case (dummy in pcorr_npt_nucleon class)
// STATUS: done
void pcorr_npt_nucleon::get_corr_momentum_averaged(gsl_vector_complex *result, const unsigned int, const unsigned int, const unsigned int, int)
{
  gsl_vector_complex_free(result);
  result = NULL; // force crash if this method is ever accidentally called
  return;
}


// Calculates the (lexicographical) index corresponding to a given n-tupel of insertions indices (assuming them to be arranged in a list without striding!).
// The size of the 'insertion[]' array most correspond to the tensor rank of the insertion.
// STATUS: done
unsigned int pcorr_npt_nucleon::get_insertion_lex_index(unsigned int insertion[])
{
  unsigned int index=0;
  for (unsigned int k=rank; k>0; k--)
  {
    index+=(unsigned int)gsl_sf_choose(insertion[rank-k] + k - 1, k);
  }
  return index;
}


// Returns the tensor rank of the insertion
// STATUS: done
unsigned int pcorr_npt_nucleon::get_rank()
{
  return rank;
}


// Destructor for class pcorr_npt
// STATUS: done
pcorr_npt_nucleon::~pcorr_npt_nucleon()
{
  const unsigned int insertion_index_size = (unsigned int)gsl_sf_choose(3+rank, rank);
  for (unsigned int i=0; i<(no_momenta * insertion_index_size); i++)
  {
    gsl_vector_complex_free(corr[i]);
  }
  delete [] corr;
  return;
}

