#ifndef STRING_UTILITIES_H_
#define STRING_UTILITIES_H_

#include <string>
#include <vector>

namespace zerork
{
namespace utilities
{

const char WHITESPACE[]=" \f\n\r\t\v"; // GNU C standard whitespace
                                       // isspace(char) returns true

// Return a copy of the input string and convert all the letters [A,Z] to
// either upper or lower case.
std::string GetUpperCase(const std::string &input_string);
std::string GetLowerCase(const std::string &input_string);

// Return the total number of times that the character c appears in the
// input_string.
size_t GetCharCount(const std::string &input_string, const char c);


// bool StringIsInt(const std::string &input_string);
// To return the TRUE, the following must be TRUE:
// 
//   (1) non-zero length
//   (2) all characters valid for number representation - 
//       no whitespace allowed
//   (3) not just a prefix/sign symbol
//   (4) the string after the prefix sign/symbol (if it exists) 
//       only contains decimal digits
//   (5) The C-string conversion function strtod is called and the value of 
//       the C macro errno must be zero (from errno.h).
//   (6) The value of the double returned from strtod in (5) is within the
//       bounds (double)INT_MIN - 0.5 < input_number < (double)INT_MAX + 0.5
//
//       ERANGE is the only flag mentioned that can be set from the online
//       documentation for strtod.  This may be platform dependent, but
//       on ubuntu 14 LTS ERANGE is set for overflow and underflow.
//       Thus, the expectation is that overflow and underflow doubles return
//       false.
bool StringIsInt(const std::string &input_string);


// bool StringIsDouble(const std::string &input_string);
// To return the TRUE, the following must be TRUE:
// 
//   (1) non-zero length
//   (2) all characters valid for number representation - 
//       no whitespace allowed
//   (3) not just a prefix/sign symbol
//   (4) no more than one decimal point
//   (5) no more than one exponent character 'e' or 'E'
//   (6) the string (excluding the prefix) before the exponent character
//       (if it exists) is nonzero and only contains decimal digits or a
//       a decimal point
//   (7) the string after the exponent character is a valid integer based
//       on a call to StringIsInt in this library
//   (8) The C-string conversion function strtod is called and the value of 
//       the C macro errno must be zero (from errno.h).
//
//       ERANGE is the only flag mentioned that can be set from the online
//       documentation for strtod.  This may be platform dependent, but
//       on ubuntu 14 LTS ERANGE is set for overflow and underflow.
//       Thus, the expectation is that overflow and underflow doubles return
//       false.
bool StringIsDouble(const std::string &input_string);


// Split the input string into a vector of strings corresponding to the
// the fields separated by the delimiters or the end of the line.
//
// Return: split_string.size()
size_t SplitStringToVector(const std::string &input_string,
                           const std::string &delimiters,
                           std::vector<std::string> *split_string);

// Split the input string into a vector of strings corresponding to the
// based on vectors for the start and stop positions.
//
// A field is added as an element to the split_string vector only if
// stop_positions[j] > start_positions[j]. This means start_positions[j]
// must be less than the length of the input_string and can not be 
// std::string::npos.  If stop_positions[j] >= length of the input_string
// or set to std::string::npos, then the field records to the end of the
// input_string.
//
// 
//    
// Return: split_string.size()
size_t SplitStringToVector(const std::string &input_string,
                           const std::vector<size_t> &start_positions,
                           const std::vector<size_t> &stop_positions,
                           std::vector<std::string> *split_string);


// Split the input string at the first appearance of the delimiter and store
// all the characters before the delimiter in the substring output.  If the 
// delimiter is not found, then the substring is set to be equal to the input
// string.  This function is useful for text processing where a comment
// character or set of characters is used to indicate that the text that
// follows is to be ignored.
//
//   if(match_exactly == false), then the first instance of any of the
//                               delimiter characters is used to split the
//                               input string.
//   if(match_exactly == true),  then the first instance of the entire
//                               delimiter string is used to split the
//                               input string.
// Note if the input_string and the substring are the same object, then
// a copy of the input string is made to enable safe copying.
//
// Return: length of the substring
size_t SplitBeforeDelimiter(const std::string &input_string,
                            const std::string &delimiter,
                            const bool match_exactly,
                            std::string *substring);

// Erase any and all characters specified in the remove_chars string from the
// input string. Note if the input_string and the substring are the same
// object then a copy of the input string is made to enable safe copying.
//
// Return: length of the substring
size_t EraseAllCharacters(const std::string &input_string,
                          const std::string &erase_chars,
                          std::string *substring);



} // namespace utilities
} // namespace advcomb
#endif
