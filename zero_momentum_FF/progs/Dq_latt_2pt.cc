// AUTHOR: Konstantin Ottnad
//
// DATE: 20131212
//
// Calculates 2pt functions and Dq_(x,y,z) lattice derivatives (i.e. "sin(x)" in FT) in momentum space using multiple radii for the FT.

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
unsigned int c_min = 0, c_max = 0, c_step = 1;

void usage()
{
  cout << endl << "This program calculates 2pt functions and Dq_(x,y,z) lattice derivatives (i.e. factor 'sin(x)' in FT)  in momentum space using multiple radii for the FT.";
  cout << endl << "This allows to check for convergence of the derivative correlator in spatial direction.";
  cout << endl << "Option are:";
  cout << endl << "-f (optional) file name (postfix) for correlator (and Dq derivative of correlator); string is appended to 'Dq_<x,y,z>.2pt.r<radius>.' for the latter) (default: stdout)";
  cout << endl << "-i (required) file name for correlator data in position space";
  cout << endl << "-m (required) file name for momenta list";
  cout << endl << "-n (optional) stepsize for cutoff used for the FT (default = 2) (*)";
  cout << endl << "-c (optional) minimal cutoff for FT (default = L) (**)";
  cout << endl << "-C (optional) maximal radius for FT (default = L) (**)" << endl;
  cout << endl << "(*)  Note that the FT is performed in the interval (c_min/2,c_max/2], s.t. n should be larger than one!";
  cout << endl << "(**) Need to set either both or none of the two options. This is because 'L' is only determined at runtime from input file" << endl;
  return;
}

int parse_cmdline_params(int argc, char **argv)
{
  int i = 0;
  while ((i = getopt(argc,argv,"h?f:i:m:n:c:C:")) !=-1)
  {
    switch (i)
    {
      case 'c':
        c_min = atoi(optarg);
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
      case 'n':
        c_step = atoi(optarg);
      break; 
      case 'C':
        c_max = atoi(optarg);
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
  if ((c_min>0)&&(c_max==0))
  {
    cout << endl << "Fatal: c_min ('-c') set but not c_max ('-C'). ";
    return -3;
  }
  if ((c_min==0)&&(c_max>0))
  {
    cout << endl << "Fatal: c_max ('-C') set but not c_min ('-c'). ";
    return -4 ;
  }
  if (c_min>c_max)
  {
    cout << endl << "Warning: c_min ('-c') larger than c_max ('-C'). Exchanging values ...";
    unsigned int tmp = c_min;
    c_min = c_max;
    c_max = tmp;
  }
  if (c_step<1)
  {
    cout << endl << "Warning: c_step ('-n') must be finite. Adjusting to default value c_step = 2 ...";
    c_step = 2;
  }
  if (c_step==1)
  {
    cout << endl << "Note: c_step==1 is pointless, correcting to c_step=2 ...";
    c_step = 2;
  }
  return 0;
}


int main(int argc, char *argv[])
{
  if (parse_cmdline_params(argc, argv)!=0) return -1;
  xcorr_2pt *corr = new xcorr_2pt(input_file, output_file, momenta_list);
  if (c_max>(unsigned)corr->get_L())
  {
    cout << endl << "Warning: c_max ('-R') larger than L. Adjusting c_max=L.";
    c_max = corr->get_L();
  }
  if ((c_min==0)&&(c_max==0))
  {
    c_min = corr->get_L();
    c_max = corr->get_L();
  }
  for (unsigned int cutoff=c_min; cutoff<=c_max; cutoff+=c_step) // loop over all radii for FT
  { 
    corr->FT_Dq_latt(cutoff);
    corr->write_results();
  }
  delete corr;
  cout << endl << endl;
  return 0;
}
