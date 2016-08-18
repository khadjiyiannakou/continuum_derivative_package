// AUTHOR: Konstantin Ottnad
//
// DATE: 20140425
//
// Implements:
// - simple parser class for pure text based input files (uses 'text_reader' as vase class) containing values of the type "token = values(s)".
//   It uses template methods to deal with various different types of value(s). Note, however, that only 'parse_token()' requires specializations in case
//   of non-standart types (i.e. types that cannot be handled directly by sstream)
//   If searching for a specific token we always return the first occurrence in the list and then delete this occurrence from the list. 
//   This allows to search for the same token multiple times, possibly returning different values. However, a particular value can only be read ONCE, as the 
//   corresponding line is deleted after reading it...

#ifndef parser_hh
#define parser_hh

#include <sstream>
#include "./io.hh"

class text_parser: public text_reader
{
  protected:
    string *data;
    template <typename Type> int assign_value(const string val, Type **value); // template method for converting strings to resulting values for a given token. Should return 1 if successful and 0 otherwise.

  public:
    text_parser(string fname, unsigned int sskip, unsigned int eskip);
    virtual ~text_parser();
    template <typename Type> int parse_token(const string token, Type **value); // searches for a single token, and returns 1 if successful and 0 otherwise. If the token is found, the method either overwrites a default 'value' allocated by the caller (this corresponds to a optional value) or allocates 'value' and returns the corresponding value in this variable (this corresponds to a required parameter). Otherwise it returns either 'value'=NULL or the default 'value'. Note that an uninitialized (i.e. required) value MUST be set to NULL by the caller!
    template <typename Type> int parse(const string token[], const unsigned int no_token, Type ***value); // searches for a list of tokens and returns the ID of the last token if no error occurred. If a token[i] is found, the method either overwrites a default 'value[i]' allocated by the caller (this corresponds to a optional value) or allocates 'value[i]' and returns the corresponding value in this variable (this corresponds to a required parameter). Otherwise it returns either 'value[i]'=NULL or the default 'value[i]'. If a required token is not found (or any other problem occurs) the method returns the ID of the last valid token.
};

// include implementation of template methods

#include "../lib/parser.tcc"
#endif
