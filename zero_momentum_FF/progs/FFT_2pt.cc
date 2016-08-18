// AUTHOR: Konstantin Ottnad
//
// DATE: 20131212
//
// Calculates 2pt functions and Dq_(x,y,z) continuum (= "x" in FT) + lattice (i.e. "sin(x)" in FT) derivatives in momentum space using multiple radii for the FT.

#include <complex>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include "../include/correlator.hh"
#include "../include/io.hh"
#include "../include/projectors.hh"

using namespace std;

string input_file = "", output_file = "", momenta_list = ""; // file names
int isospin = ppm, operators = spin_1_1 , projector = 13;
int write_spatial_corrs = 0;

void usage()
{
  cout << endl << "This program calculates the full lattice momentum dependence in x,y,z directions for a 2pt function, while projecting the orthogonal components to zero momentum.";
  cout << endl << "Option are:";
  cout << endl << "-b (optional) enables writing of spatial x,y,z block corrs with orthogonal components projected to zero momentum (default = off)";
  cout << endl << "-f (optional) file name (postfix) for correlator (and Dq derivative of correlator); string is appended to 'Dq_<x,y,z>.2pt.r<radius>.' for the latter) (default: stdout)";
  cout << endl << "-i (required) file name for correlator data in position space";
  cout << endl << "-m (required) file name for momenta list";
  cout << endl << "-I (optional) set the isospin, i.e. ppm for proton and pmm for neutron                 (default = ppm)";
  cout << endl << "-O (optional) set the interpolating operators at source and sink for the 2pt functions (default = 1-1)";
  cout << endl << "-P (optional) set the projector; for a list see 'projectors.hh'                        (default = 13, i.e. 'Gamma_4')";
  cout << endl << "HINT: Do not forget to set the output filename postfix according to the isospin, interpolating operators and projector!" << endl;
  return;
}


int parse_cmdline_params(int argc, char **argv)
{
  int i = 0;
  while ((i = getopt(argc,argv,"h?bf:i:m:I:O:P:")) !=-1)
  {
    switch (i)
    {
      case 'b':
        write_spatial_corrs = 1;
      break;
      case 'f':
        output_file = optarg;
      break;
      case 'i':
        input_file = optarg;
      break;
      case 'm':
        momenta_list = optarg;
      break;
      case 'I':
        isospin = get_isospin_nucleon(optarg);
      break;
      case 'O':
        operators = get_spin_nucleon(optarg);
      break;
      case 'P':
        projector = atoi(optarg);
      break;
      case 'h':
      case '?':
      default:
        usage();
        return -1;
    }
  }
  if (output_file=="")
  {
    cout << endl << "WARNING: no output file given ('-f'); writing results to stdout.";
  }
  if (input_file=="")
  {
    cout << endl << "FATAL: no input file given ('-i').";
    usage();
    return -1;
  }
  if (momenta_list=="")
  {
    cout << endl << "FATAL: no momenta list given ('-m'). ";
    usage();
    return -2;
  }
  if (isospin<0)
  {
    cout << endl << "WARNING: illegal identifier passed for nucleon isospin ('-I'). Using default = ppm";
    isospin = ppm;
  }
  if (operators<0)
  {
    cout << endl << "WARNING: illegal identifier passed for nucleon interpolating operators ('-O'). Using default = 1-1.";
    operators = spin_1_1;
  }
  if ((projector<0)||(projector>=(int)no_projectors))
  {
    cout << endl << "WARNING: illegal projectors index passed; valid indices are [0..." << (no_projectors-1) << "]. Using default = 13.";
    projector = 13;
  }
  return 0;
}


int main(int argc, char *argv[])
{
  if (parse_cmdline_params(argc, argv)!=0) return -1;
  init_projectors();
  xcorr_2pt *corr = new xcorr_2pt(input_file, output_file, momenta_list);
  if (write_spatial_corrs)
  {
    corr->block_xyz(projector, (isospin_nucleon)isospin, (spin_nucleon)operators);
    corr->write_results();
  }
  corr->FFT_xyz(projector, (isospin_nucleon)isospin, (spin_nucleon)operators);
  corr->write_results();
  delete corr;
  destroy_projectors();
  cout << endl << endl;
  return 0;
}
