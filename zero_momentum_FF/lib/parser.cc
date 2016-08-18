// AUTHOR: Konstantin Ottnad
//
// DATE: 20140425
//
// Implements:
// - simple parser class for pure text based input files (uses 'text_reader' as vase class) containing values of the type "token = values(s)".
//   It uses template methods to deal with various different types of value(s). Note, however, that only 'assign_value()' requires specializations in case
//   of non-standart types (i.e. types that cannot be handled directly by sstream)
//   If searching for a specific token we always return the first occurrence in the list and then delete this occurrence from the list. 
//   This allows to search for the same token multiple times, possibly returning different values. However, a particular value can only be read ONCE, as the 
//   corresponding line is deleted after reading it...
//
// NOTES: 
// - For implementation of template methods refer to './lib/parser.tcc'

#include "../include/parser.hh"

using namespace std;

// Constructor for 'text_parser'
// STATUS: done
text_parser::text_parser(string fname, unsigned int sskip, unsigned int eskip) : text_reader(fname, sskip, eskip)
{
  data = new string [lines];
  unsigned int line = 0;
  while ((!stream.eof())&&(line<lines))
  {
    getline(stream, data[line]); // load the entire file into main memory
    line++;
  }
  return;
}


// Destructor for 'text_parser'
// STATUS: done
text_parser::~text_parser()
{
  delete [] data;
  return;
}
