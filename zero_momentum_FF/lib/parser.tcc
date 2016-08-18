// AUTHOR: Konstantin Ottnad
//
// DATE: 20140427
//
// Implements:
// - Template methods for class 'text_parser' defined in './include/parser.hh'. See also './lib/parser.cc'. 
//
// This file is included by './include/parser.hh'.


// General version for standard types, e.g. double. Returns always 1 (cannot fail, although result may be useless)
// STATUS: done
template <typename Type> int text_parser::assign_value(const string val, Type **value)
{
  stringstream s;
  s.str(val);
  if (*value==NULL) *value = new Type; // only allocate if it is not an optional value
  s >> **value;
  return 1;
}


// Searches for a single token, and returns 1 if successful and 0 otherwise. If the token is found, the method either overwrites a 
// default 'value' allocated by the caller (this corresponds to an optional value) or allocates 'value' and returns the corresponding 
// value in this variable (this corresponds to a required parameter). Otherwise it returns either 'value'=NULL or the default 'value'.
// Note that an uninitialized (i.e. required) value MUST be set to NULL by the caller!
// STATUS: done
template <typename Type> int text_parser::parse_token(const string token, Type **value)
{
  if (stream.is_open())
  {
    stringstream s;
    unsigned int fail = 0;
    for (int i=0; i<(int)lines; i++)
    {
      if (data[i].find(token)!=string::npos) // check for occurrence of 'token' in the list
      {
        s.str(data[i]);
        for (unsigned int j=0; j<s.str().find(token); j++)
        {
          if (s.str().at(j)!=' ') // check that there is nothing but (possibly) spaces in front of the token
          {
            fail = 1;
            break;
          }
        }
        if (fail) continue; // possibly found a different token or nonsense
        s.str(s.str().erase(0, s.str().find(token) + token.length()));  // if everything is ok, delete the 'token' substring
        if (s.str().find("=")==string::npos) continue; // enforce the existence of '=', otherwise discard the token and the value...
        for (unsigned int j=0; j<s.str().find("="); j++)  // check that there is nothing but spaces between 'token' and '=' ...
        {
          if (s.str().at(j)!=' ')
          {
            fail = 1;
            break;
          }
        }
        if (fail) continue;
        s.str(s.str().erase(0, s.str().find("=") + 1)); // delete the "=" and possible spaces
        if (assign_value(s.str(), value))
        {
	  data[i] = ""; // finally delete the corresponding line
          return 1;  // token found, value assigned -- return state successful 
        }
      }
    }
  }
  if (*value!=NULL) return 1; // if we have been looking for an optional parameter, it is not a problem that we did not find anything
  return 0; // if we reach this point either the correct token was not found, or the value could not be assigned or any other error occurred.
}


// Searches for a list of tokens and returns the ID of the last token if no error occurred. If a token[i] is found, the method either
// overwrites a default 'value[i]' allocated by the caller (this corresponds to a optional value) or allocates 'value[i]' and 
// returns the corresponding value in this variable (this corresponds to a required parameter). Otherwise it returns either 
// 'value[i]'=NULL or the default 'value[i]'.  
// If a required token is not found (or any other problem occurs) the method returns the ID of the last valid token.
// Note that uninitialized (i.e. required) values MUST be set to NULL by the caller!
// STATUS: done
template <typename Type> int text_parser::parse(const string token[], const unsigned int no_token, Type ***value)
{
  unsigned int count = 0;
  if (stream.is_open())
  {
    while (count<no_token)
    {
      if (parse_token(token[count], &((*value)[count])))
      {
        count++;
      }
      else
      {
        break; // abort if a required token was not found / an error occurred
      }
    }
  }
  return count-1; // return the number of the last successfully read token+values pairs
}
