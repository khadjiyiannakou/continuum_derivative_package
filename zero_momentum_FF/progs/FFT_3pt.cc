// AUTHOR: Konstantin Ottnad
//
// DATE: 20131212
//
// Calculates 3pt functions and Dq_(x,y,z) continuum (= "x" in FT) + lattice (i.e. "sin(x)" in FT) derivatives in momentum space using multiple radii for the FT

#include <complex>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include "../include/correlator.hh"
#include "../include/io.hh"

using namespace std;

string input_file = "", output_file = "", momenta_list = ""; // file names
int write_spatial_corrs = 0;


void usage()
{
  cout << endl << "This program calculates the full lattice momentum dependence in x,y,z directions for a 3pt function, while projecting the orthogonal components to zero momentum.";
  cout << endl << "Option are:";
  cout << endl << "-b (optional) enables writing of spatial x,y,z block corrs with orthogonal components projected to zero momentum (default = off)";
  cout << endl << "-f (optional) file name (postfix) for correlator (and Dq derivative of correlator); string is appended to '(Dq_<x,y,z>.) 3pt.r<radius>.' for the latter) (default: stdout)";
  cout << endl << "-i (required) file name for correlator data in position space";
  cout << endl << "-m (required) file name for momenta list"; 
  return;
}

int parse_cmdline_params(int argc, char **argv)
{
  int i = 0;
  while ((i = getopt(argc,argv,"h?bf:i:m:")) !=-1)
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
      case 'h':
      case '?':
      default:
        usage();
        return -1;
    }
  }
  if (output_file=="")
  {
    cout << endl << "Warning: no output file given ('-f'); writing results to stdout.";
  }
  if (input_file=="")
  {
    cout << endl << "Fatal: no input file given ('-i').";
    usage();
    return -1;
  }
  if (momenta_list=="")
  {
    cout << endl << "Fatal: no momenta list given ('-m'). ";
    usage();
    return -2;
  } 
  return 0;
}


int main(int argc, char *argv[])
{
  if (parse_cmdline_params(argc, argv)!=0) return -1;
  xcorr_3pt *corr = new xcorr_3pt(input_file, output_file, momenta_list);
  if (write_spatial_corrs)
  {
    corr->block_xyz();
    corr->write_results();
  }
  corr->FFT_xyz();
  corr->write_results();
  delete corr;
  cout << endl << endl;
  return 0;
}
