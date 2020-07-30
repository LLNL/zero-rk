#include <errno.h>
#include <stdlib.h>

#include <limits.h>
#include <float.h> 

//#include <stdio.h> // printf for debugging

#include "string_utilities.h"

namespace zerork 
{
namespace utilities
{
std::string GetUpperCase(const std::string &input_string)
{
  std::string copy_string = input_string;
  const size_t const_length = copy_string.size();

  for(size_t j=0; j<const_length; ++j) {
    if('a' <= copy_string[j] && copy_string[j] <= 'z') {
      copy_string[j] += 'A' - 'a';
    }
  }
  return copy_string;
}

std::string GetLowerCase(const std::string &input_string)
{
  std::string copy_string = input_string;
  const size_t const_length = copy_string.size();

  for(size_t j=0; j<const_length; ++j) {
    if('A' <= copy_string[j] && copy_string[j] <= 'Z') {
      copy_string[j] += 'a' - 'A';
    }
  }
  return copy_string;
}

size_t GetCharCount(const std::string &input_string, const char c)
{
  size_t count = 0;
  const size_t const_length = input_string.size();

  for(size_t j=0; j<const_length; ++j) {
    if(input_string[j] == c) {
      ++count;
    }
  }
  return count;
}

bool StringIsInt(const std::string &input_string)
{
  const char valid_chars[]  = "+-0123456789";
  const char valid_digits[] = "0123456789";

  if(input_string.size() == 0) {
    return false;
  }

  // check if any character isn't valid 
  size_t bad_char_pos = input_string.find_first_not_of(valid_chars,0);
  if(bad_char_pos != std::string::npos) {
    return false;
  }

  // check if there is a prefix sign
  if(input_string[0] == '+' || input_string[0] == '-') {

    if(input_string.size() <= 1) {
      return false;
    }
    // make sure all the remaining characters are digits
    bad_char_pos = input_string.find_first_not_of(valid_digits,1);
    if(bad_char_pos != std::string::npos) {
      return false;
    }  
  }

  // check that the integer is in bounds
  errno = 0;
  double input_number = strtod(input_string.c_str(),
                               NULL); // do not return end of number pointer
  if(errno != 0) {
    return false;
  }
  
  return ((double)INT_MIN - 0.5 < input_number &&
          input_number < (double)INT_MAX + 0.5);
}


bool StringIsDouble(const std::string &input_string)
{
  const char valid_chars[]  = "+-.0123456789eE";
  const char valid_mantissa[] = ".0123456789";
  std::string mantissa_check;
  size_t find_pos=0;
  size_t start_pos=0;
  size_t e_pos=0;

  if(input_string.size() == 0) {
    return false;
  }

  // check if any character isn't valid 
  find_pos = input_string.find_first_not_of(valid_chars,0);
  if(find_pos != std::string::npos) {
    return false;
  }

  // Make sure there is more than just a prefix sign
  if(input_string[0] == '+' || input_string[0] == '-') {
    // set the start position of the mantissa to one
    start_pos = 1;
    if(input_string.size() <= 1) {
      return false;
    }
  }
  // Make sure there are are no more than one decimal point and exponent char
  if(GetCharCount(input_string,'.') >= 2) {
    return false;
  }

  // Make sure there are are no more than one decimal point and exponent char
  if(GetCharCount(input_string,'e')+GetCharCount(input_string,'E') >= 2) {
    return false;
  }

  // Make sure that the characters between the prefix symbol (if present) 
  // and the exponent char (if present) form a valid mantissa
  e_pos = input_string.find_first_of("Ee",0);

  if(e_pos == std::string::npos) {
    // No exponent symbol
    mantissa_check.assign(input_string,
			  start_pos,
                          std::string::npos);
  } else if(e_pos > start_pos+1) {
    mantissa_check.assign(input_string,
                          start_pos,
                          e_pos-start_pos); 
  } else {
    return false;
  }
  find_pos = mantissa_check.find_first_not_of(valid_mantissa, 0);
  if(find_pos != std::string::npos) {
    return false;
  }

  // Make sure that what follows the exponent char is an integer
  if(e_pos != std::string::npos) {

    if(e_pos+1 >= input_string.size()) {
      return false;
    }
    if(!StringIsInt(&input_string[e_pos+1])) {
      return false;
    }
  }

  // check that the double is in bounds
  errno = 0; // Macro that can be set internally by C-library functions like 
             // strtod
  double input_number = strtod(input_string.c_str(),
                               NULL); // do not return end of number pointer
  if(errno != 0) {
    // will be set if underflow |input_string| < DBL_MIN
    //          or if overflow  |input_string| > DBL_MAX
    //   DBL_MIN = 2.2250738585072013830902327173324e-308  (%36.31e)  
    //   DBL_MAX = 1.7976931348623157081452742373170e+308  (%36.31e)
    //printf("DBL_MIN = %36.31e\n", DBL_MIN);
    //printf("DBL_MAX = %36.31e\n", DBL_MAX);

    //printf("errno set out of range for %s\n",input_string.c_str());
    return false;
  }
  
  // TODO: Increase the rigor of the check if possible because an invalid
  //       conversion returns as zero. This is why there are the additional
  //       checks above.  The check below should always return TRUE
  return (-DBL_MAX <= input_number &&
          input_number <= DBL_MAX);
}


size_t SplitStringToVector(const std::string &input_string,
                           const std::string &delimiters,
                           std::vector<std::string> *split_string)
{
  size_t field_start = 0;
  size_t field_end   = 0;
  std::string field;
  split_string->clear();
  
  if(input_string.size() == 0) {
    return split_string->size();
  }
 
  // This approach will merge all delimiters. As a consequence, there are no
  // empty string fields stored in the split_string vector
  field_start = input_string.find_first_not_of(delimiters,0);
  if(field_start == std::string::npos) {
    return split_string->size();
  }
  field_end   = input_string.find_first_of(delimiters,field_start);
  // Note that field_end is the position of the first delimiter character
  // after the last field character.  This means that (field_end-field_start)
  // is the length of the field.
  while(field_end   != std::string::npos) {

    field.assign(input_string,field_start,field_end-field_start);
    split_string->push_back(field);

    field_start = input_string.find_first_not_of(delimiters,field_end);
    if(field_start == std::string::npos) {
      // no more fields
      return split_string->size();
    }
    field_end   = input_string.find_first_of(delimiters,field_start);
  }  
  // At this point, field_start should not be string::npos, below is just
  // another check to avoid assigning an empty string
  if(field_start == std::string::npos) {
    // no more fields
    return split_string->size();
  }
  field.assign(input_string,field_start,std::string::npos);
  split_string->push_back(field);

  return split_string->size();
}


size_t SplitStringToVector(const std::string &input_string,
                           const std::vector<size_t> &start_positions,
                           const std::vector<size_t> &stop_positions,
                           std::vector<std::string> *split_string)
{
  split_string->clear();
  std::string field;
  size_t input_length = input_string.size();

  if(input_length == 0) {
    return split_string->size();
  }

  // Take the minimum size of the two position vectors
  const size_t num_positions = 
    ((start_positions.size() < stop_positions.size()) ? 
     start_positions.size()  : stop_positions.size()); 

  for(size_t j=0; j<num_positions; ++j) {

    field.clear();

    if(start_positions[j] < stop_positions[j] &&
       start_positions[j] < input_length) {

      field.assign(input_string,
                   start_positions[j],
                   stop_positions[j]-start_positions[j]); // length
    }
    split_string->push_back(field);
  }

  return split_string->size();
}


size_t SplitBeforeDelimiter(const std::string &input_string,
                            const std::string &delimiters,
                            const bool match_exactly,
                            std::string *substring)
{
  size_t field_end = 0;

  if(&input_string == substring) {
    // if the address of the string objects are the same create a copy
    // of the input string and call again
    std::string copy_string = input_string;
    return SplitBeforeDelimiter(copy_string,
                                delimiters,
                                match_exactly,
                                substring);
  }

  substring->clear();

  if(input_string.size() == 0) {
    return substring->size();
  }
  if(match_exactly) {
    // Find returns the first character of the first match.  Note that, unlike
    // find_first_of, whenever the length of the string being searched for
    // (here delimiters) is greater than one, then the entire search string
    // must match.
    field_end = input_string.find(delimiters,0);
  } else {
    // find the first instance of any of the delimiter characters
    field_end = input_string.find_first_of(delimiters,0);
  }

  // If none of the delimiter characters are found than std::string::npos
  // is returned. If std::string::npos is the length specified to assign 
  // member function than the whole input string is copied.
  substring->assign(input_string,
                    0,            // start position
                    field_end);   // length    

  return substring->size();
}

size_t EraseAllCharacters(const std::string &input_string,
                          const std::string &erase_chars,
                          std::string *substring)
{
  size_t num_tokens;
  std::vector<std::string> keep_tokens;

  if(&input_string == substring) {
    std::string copy_string = input_string;

    return EraseAllCharacters(copy_string,
                              erase_chars,
                              substring);
  }
  substring->clear();

  // use the same logic to split the string as SplitStringToVector
  num_tokens = SplitStringToVector(input_string,
                                   erase_chars,
                                   &keep_tokens);
  for(size_t j =0; j<num_tokens; ++j) {
    substring->append(keep_tokens[j]);
  }
  return substring->size();
}


} // namespace utilities 
} // namespace advcomb
