#include <iostream>                        
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <complex>
#include <stdio.h>

using namespace std;

/*
  what arguments are needed for this code to run
  1) number of timeslices 
  2) number of spatial Length
  3) source_x
  4) source_y
  5) source_z
  6) insertion file g0
  7) insertion file gx
  8) insertion file gy
  9) insertion file gz
  10) output file where we write the data

 */

const unsigned int littleEndian = 0;
const unsigned int bigEndian = 1;
typedef std::complex<double> Complex;
const unsigned int dataEndianess = 1; // this need to modified if data are in different endianess


unsigned int get_endianess()
{
  union
  {
    char c[4];
    int i;
  } word32;
  word32.i=1;
  if(word32.c[3]==1) return bigEndian;
  return littleEndian;
}

void swap_byte_order_8(char *block, int64_t block_length)
{
  register char *i,*j,*k;
  char buf;
  char *bound = block + (block_length*8);
  for(i=block; i<bound; i+=8)
    {
      j=i;
      k=j+7;
      buf = *j;
      *j  = *k;
      *k  = buf;
      j++;
      k--;
      buf = *j;
      *j  = *k;
      *k  = buf;
      j++;
      k--;
      buf = *j;
      *j  = *k;
      *k  = buf;
      j++;
      k--;
      buf = *j;
      *j  = *k;
      *k  = buf;
    }
  return;
}


int main(int argc, char* argv[]){

  if(argc != 11){
    cout << "correct order -> 1) number of timeslices, 2) number of spatial Length, 3) source_x, 4) source_y, 5) source_z,  6) insertion file g0, 7) insertion file gx, 8) insertion file gy, 9) insertion file gz, 10) output file where we write the data\n";
    wcerr << "wrong number of input files provided\n";
    exit(EXIT_FAILURE);
  }

  int ntslices = atoi(argv[1]);
  int Ls = atoi(argv[2]);
  int src_x = atoi(argv[3]);
  int src_y = atoi(argv[4]);
  int src_z = atoi(argv[5]);

  fstream f_g0, f_gx, f_gy, f_gz, f_out;
  f_g0.open(argv[6], fstream::in | fstream::binary );
  f_gx.open(argv[7], fstream::in | fstream::binary );
  f_gy.open(argv[8], fstream::in | fstream::binary );
  f_gz.open(argv[9], fstream::in | fstream::binary );
  f_out.open(argv[10], fstream::out);

  if(!f_g0 || !f_gx || !f_gy || !f_gz){
    wcerr << "Error cannot open files to read binary data\n";
    exit(EXIT_FAILURE);
  }

  if(!f_out){
    wcerr << "Error cannot open file to write output\n";
    exit(EXIT_FAILURE);
  }

  int64_t length_g0 = 0;
  int64_t length_gx = 0;
  int64_t length_gy = 0;
  int64_t length_gz = 0;

  f_g0.seekg(0,f_g0.end);
  f_gx.seekg(0,f_gx.end);
  f_gy.seekg(0,f_gy.end);
  f_gz.seekg(0,f_gz.end);
  
  // get the size of the binary file
  length_g0 = f_g0.tellg();
  length_gx = f_gx.tellg();
  length_gy = f_gy.tellg();
  length_gz = f_gz.tellg();
  
  // returns the file stream at the beginning  
  f_g0.seekg(0,f_g0.beg);
  f_gx.seekg(0,f_gx.beg);
  f_gy.seekg(0,f_gy.beg);
  f_gz.seekg(0,f_gz.beg);

  const int64_t expectedSize = ntslices * Ls * Ls * Ls * 2 * 8; // data is stored in double precision
  
  if(length_g0 != expectedSize){
    wcerr << "Size of the file g0 does not agree with expected size";
    exit(EXIT_FAILURE);
  }

  if(length_gx != expectedSize){
    wcerr << "Size of the file gx does not agree with expected size";
    exit(EXIT_FAILURE);
  }

  if(length_gy != expectedSize){
    wcerr << "Size of the file gy does not agree with expected size";
    exit(EXIT_FAILURE);
  }

  if(length_gz != expectedSize){
    wcerr << "Size of the file gz does not agree with expected size";
    exit(EXIT_FAILURE);
  }

  // allocate memory to read binary data
  char *data_g0 = NULL;
  char *data_gx = NULL;
  char *data_gy = NULL;
  char *data_gz = NULL;

  try{ 
    data_g0 = new char[expectedSize];
    data_gx = new char[expectedSize];
    data_gy = new char[expectedSize];
    data_gz = new char[expectedSize];
  } // check for bad allocation
  catch (bad_alloc& ba){
    wcerr << "Error caught:" << ba.what() << endl;
    exit(EXIT_FAILURE);
  }

  // read binary files
  f_g0.read(data_g0,expectedSize);
  f_gx.read(data_gx,expectedSize);
  f_gy.read(data_gy,expectedSize);
  f_gz.read(data_gz,expectedSize);


  if(get_endianess() != dataEndianess){
    cout << "Machine and data have different endianess\n";
    cout << "performing byte swap to data\n";
    swap_byte_order_8(data_g0, expectedSize/8);
    swap_byte_order_8(data_gx, expectedSize/8);
    swap_byte_order_8(data_gy, expectedSize/8);
    swap_byte_order_8(data_gz, expectedSize/8);
    cout << "finished byte swap\n";
  }

  Complex *xcorr_g0 = (Complex*) data_g0;
  Complex *xcorr_gx = (Complex*) data_gx;
  Complex *xcorr_gy = (Complex*) data_gy;
  Complex *xcorr_gz = (Complex*) data_gz;

  Complex Dq_corr[3*ntslices*4]; // 3: derivative directions , ntslices: number of time slices, 4: number of insertion operators

  for(int i = 0 ; i < 3*ntslices*4 ; i++){
    Dq_corr[i].real() = 0.;
    Dq_corr[i].imag() = 0.;
  }

  // summation over the volume with \vec{x}

  for(int it = 0 ; it < ntslices ; it++){
    for(int x=Ls/2; x>-Ls/2; x--)                  
      for(int y=Ls/2; y>-Ls/2; y--)
	for(int z=Ls/2; z>-Ls/2; z--)
	  {
	    // map data in the [0,L) range
	    const int Xmap = ((x>=0 ? x : x+Ls) + src_x)%Ls;
	    const int Ymap = ((y>=0 ? y : y+Ls) + src_y)%Ls;
	    const int Zmap = ((z>=0 ? z : z+Ls) + src_z)%Ls;
	    // get the index
	    const int64_t ix = it*Ls*Ls*Ls + Xmap*Ls*Ls + Ymap*Ls + Zmap;
	    // derivative on q in the x direction
	    Dq_corr[0*ntslices*4 + it*4 + 0] += Complex(0,+1) * ((double)x) * xcorr_g0[ix]; // in Konstantin code it has a plus sign
	    Dq_corr[0*ntslices*4 + it*4 + 1] += Complex(0,+1) * ((double)x) * xcorr_gx[ix];
	    Dq_corr[0*ntslices*4 + it*4 + 2] += Complex(0,+1) * ((double)x) * xcorr_gy[ix];
	    Dq_corr[0*ntslices*4 + it*4 + 3] += Complex(0,+1) * ((double)x) * xcorr_gz[ix];

	    // derivative on q in the y direction
	    Dq_corr[1*ntslices*4 + it*4 + 0] += Complex(0,+1) * ((double)y) * xcorr_g0[ix];
	    Dq_corr[1*ntslices*4 + it*4 + 1] += Complex(0,+1) * ((double)y) * xcorr_gx[ix];
	    Dq_corr[1*ntslices*4 + it*4 + 2] += Complex(0,+1) * ((double)y) * xcorr_gy[ix];
	    Dq_corr[1*ntslices*4 + it*4 + 3] += Complex(0,+1) * ((double)y) * xcorr_gz[ix];

	    // derivative on q in the z direction
	    Dq_corr[2*ntslices*4 + it*4 + 0] += Complex(0,+1) * ((double)z) * xcorr_g0[ix];
	    Dq_corr[2*ntslices*4 + it*4 + 1] += Complex(0,+1) * ((double)z) * xcorr_gx[ix];
	    Dq_corr[2*ntslices*4 + it*4 + 2] += Complex(0,+1) * ((double)z) * xcorr_gy[ix];
	    Dq_corr[2*ntslices*4 + it*4 + 3] += Complex(0,+1) * ((double)z) * xcorr_gz[ix];
	  } // for loop over spatial volume
  } // for loop over time slices

  ios::fmtflags old_settings = f_out.flags();
  for(int idq = 0 ; idq < 3 ; idq++)
    for(int it = 0 ; it < ntslices ; it++)
      for(int imu = 0 ; imu < 4 ; imu++){
	//	fprintf(f_out,"%d \t %+d %+d %+d \t %+e %+e %d\n",it,0,0,0,Dq_corr[idq*it*4 + it*4 + imu].real(), Dq_corr[idq*it*4 + it*4 + imu].imag(),imu);
	f_out << it << " +0 +0 +0 " << scientific << showpos << Dq_corr[idq*ntslices*4 + it*4 + imu].real() << " "<< Dq_corr[idq*ntslices*4 + it*4 + imu].imag();
	f_out.flags(old_settings);
	f_out << " " << imu << endl;
      }

  delete data_g0;
  delete data_gx;
  delete data_gy;
  delete data_gz;

  return 0;
}

/*
  for(int i = 0 ; i < 10 ; i++)
    cout << scientific << showpos << xcorr_gz[i].real() << " " << xcorr_gz[i].imag() << endl;
*/

    /*
  cout << "g0 file size=" << length_g0 << endl;
  cout << "gx file size=" << length_gx << endl;
  cout << "gy file size=" << length_gy << endl;
  cout << "gz file size=" << length_gz << endl;
    */
