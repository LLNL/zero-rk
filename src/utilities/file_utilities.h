#ifndef FILE_UTILITIES_H_
#define FILE_UTILITIES_H_

#include <string>
#include <vector>
#include <iostream>

namespace zerork
{
namespace utilities 
{

extern const char* null_filename;

// Returns false if fopen(file_name.c_str(),"r") is NULL, and true otherwise.
bool FileIsReadable(const std::string &file_name);


// Returns false if fopen(file_name.c_str(),"w") is NULL, and true otherwise.
bool FileIsWritable(const std::string &file_name);

// std::istream& GetAnyLine(istream &input_stream, std::string *line);
//
// Replacement for C++ std::getline for strings to extract a line with any
// of the following end-of-line characters:
//
//   '\n'   : new line (unix/linux)
//   '\r'   : carriage return (old mac)
//   '\r\n' : new line and carriage return (windows/dos) 
//   EOF    : end-of-file
//
// Also handles the case when the last line has no end-of-line character.
// Note that the input argument is a pointer to a string instead of a
// reference to the string, which differs from the classic getline.
std::istream& GetAnyLine(std::istream &input_stream, std::string *line);


// Search the input_stream for the keyword.  If case_sensitive is true,
// then the keyword must match exactly at some point in the input_stream.
// If case_sensitive is false, then the keyword and input_stream line must
// match exactly after both are converted to lower case.
//
// On return:
//
//   line - set to the line in the input stream containing the first
//          appearance of the keyword. If no keyword is found, or the
//          keyword is of size zero then line is an empty string
//   line_num - (optional) set to the number of lines read while searching
//              for the keyword from the current position in the input_stream. 
//              Note that this is NOT the total lines read from the 
//              input_stream if the function is called repeatedly for the 
//              same input_stream.
//   return value - reference to the input_stream containing its current
//                  state  
std::istream& GetFirstLineWithKeyword(std::istream &input_stream,
                                      const std::string &keyword,
                                      const std::string &comment_chars,
                                      const bool case_sensitive,  
                                      std::string *line,
                                      size_t *line_num);
std::istream& GetFirstLineWithKeyword(std::istream &input_stream,
                                      const std::string &keyword,
                                      const std::string &comment_chars,
                                      const bool case_sensitive,  
                                      std::string *line);


// The Logger class allows a single printf or flush command to be sent to all
// log streams. It also allows for the suppression of all log streams by
// setting the log_filename to null_filename and use_stdout and use_stderr to
// false.  
class Logger
{
 public: 
  Logger(); // default constructor sends log streams to stdout

  Logger(const std::string &log_filename); // file name only constructor
                                           // no output to stdout or stderr
                                           // unless there is an error opening
                                           // log_filename

  Logger(const std::string &log_filename,
         const bool use_stdout,
         const bool use_stderr);
  ~Logger();
  int PrintF(const char *format, ...) const;
  int FFlush() const;  

 private:
  void InitializeStreams(const std::string &log_filename,
                         const bool use_stdout,
                         const bool use_stderr);
  bool any_log_stream_;
  bool use_stdout_;
  bool use_stderr_;
  FILE *log_file_ptr_;
};



} // namespace utilities
} // namespace advcomb

#endif

