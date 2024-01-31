#include <stdio.h>
#include <stdarg.h>
#include <time.h>


#include "string_utilities.h"
#include "file_utilities.h"

namespace zerork
{
namespace utilities 
{

#ifdef _WIN32
const char* null_filename = "nul";
#else
const char* null_filename = "/dev/null";
#endif

bool FileIsReadable(const std::string &file_name)
{
  FILE *fptr = fopen(file_name.c_str(),"r");
  if(fptr == NULL) {
    return false;
  }
  fclose(fptr);
  return true;
}

bool FileIsWritable(const std::string &file_name)
{
  FILE *fptr = fopen(file_name.c_str(),"w");
  if(fptr == NULL) {
    return false;
  }
  fclose(fptr);
  return true;
}

std::istream& GetAnyLine(std::istream &input_stream, std::string *line)
{
  line->clear();

  // The stream sentry is needed when extracting or inserting data in the
  // underlying buffer of the stream object.  Using the standard operators or
  // methods of an istream (or ostream) manages the sentry object internally.
  // The stream buffer is accessed directly here to speedup the character
  // by character reading of each line.
  std::istream::sentry stream_buffer_sentry(input_stream, 
                                            true); // T = don't skip whitespace
  std::streambuf *stream_buffer = input_stream.rdbuf();

  while(1) {
    // get the current character and advance to the next position
    int c = stream_buffer->sbumpc();
    switch (c) {

    case '\n':
      return input_stream;

    case '\r':
      // check if the next character is the newline indicating a windows/dos
      // line ending
      if(stream_buffer->sgetc() == '\n') {
        stream_buffer->sbumpc();
      }
      return input_stream;

    case EOF:
      // no end-of-line character found on last line
      if(line->size() == 0) {
        input_stream.setstate(std::ios::eofbit);
      }
      return input_stream;

    default:
      line->append(1,(char)c);
    }
  }
}


std::istream& GetFirstLineWithKeyword(std::istream &input_stream,
                                      const std::string &keyword,
                                      const std::string &comment_chars,
                                      const bool case_sensitive,  
                                      std::string *line,
                                      size_t *line_num)
{
  std::string key;
  std::string search_line;
  std::string file_line;
  size_t find_pos;
  size_t count = 0;

  line->clear();

  if(keyword.size() > 0) {

    if(case_sensitive) {
      key = keyword;
    } else {
      key = GetLowerCase(keyword);
    }

    while(GetAnyLine(input_stream, &file_line)) {

      ++count;

      SplitBeforeDelimiter(file_line,
                           comment_chars,
                           false, // match any comment_chars
                           &search_line);

      if(!case_sensitive) {
        search_line = GetLowerCase(search_line);
      }

      find_pos = search_line.find(key, 0);
      if(find_pos != std::string::npos) {
        line->assign(file_line);
        break;
      }
    }
  }

  if(line_num != NULL) {
    *line_num = count;
  }
  return input_stream;
}

std::istream& GetFirstLineWithKeyword(std::istream &input_stream,
                                      const std::string &keyword,
                                      const std::string &comment_chars,
                                      const bool case_sensitive,  
                                      std::string *line)
{
  return GetFirstLineWithKeyword(input_stream,
                                 keyword,
                                 comment_chars,
                                 case_sensitive,
                                 line,
                                 NULL);
}

// Default Logger constructor, use stdout
Logger::Logger()
{
  any_log_stream_  = true;
  use_stdout_      = true;
  use_stderr_      = false;
  log_file_ptr_    = NULL;
}

// Logger constructor for file only. Note if file does not exist, then
// it defaults to outputing the log streams to stderr.
Logger::Logger(const std::string &log_filename) 
{
  InitializeStreams(log_filename, false, false);
}

Logger::Logger(const std::string &log_filename,
               const bool use_stdout,
	       const bool use_stderr)
{
  InitializeStreams(log_filename, use_stdout, use_stderr);
}

Logger::~Logger()
{
  FFlush();
  if(log_file_ptr_ != NULL) {
    // record closing log time
    time_t raw_time;
    struct tm *time_info;

    time(&raw_time);
    time_info = localtime(&raw_time);
  
    fprintf(log_file_ptr_,
            "#------------------------------------------------------------------------------\n"
            "# Logger closed at %s"
            "#------------------------------------------------------------------------------\n",
            asctime(time_info));   // asctime returns C-string with '\n\0'
                                   // ending
    fflush(log_file_ptr_);
    fclose(log_file_ptr_);
  }
}

void Logger::InitializeStreams(const std::string &log_filename,
                               const bool use_stdout,
                               const bool use_stderr)
{
  any_log_stream_ = false;
  use_stdout_     = use_stdout;
  use_stderr_     = use_stderr;
  log_file_ptr_   = NULL;

  if(log_filename != null_filename && log_filename.size() > 0) {

    log_file_ptr_ = fopen(log_filename.c_str(),"a");

    if(log_file_ptr_ == NULL) {
      fprintf(stderr,"# ERROR: In Logger constructor,\n");
      fprintf(stderr,"#        could not open log file %s for write/append.\n",
              log_filename.c_str());
      fprintf(stderr,"#        Sending Logger file stream to stderr.\n");
      use_stderr_ = true;
    } else {
      // record starting log time
      time_t raw_time;
      struct tm *time_info;

      time(&raw_time);
      time_info = localtime(&raw_time);
  
      fprintf(log_file_ptr_,
              "#------------------------------------------------------------------------------\n"
              "# Logger opened at %s"
              "#------------------------------------------------------------------------------\n",
              asctime(time_info)); // asctime returns C-string with '\n\0'
                                   // ending
      fflush(log_file_ptr_);
    }
  } // end if(log_filename != null_filename)

  if(use_stdout_ || use_stderr_ || log_file_ptr_ != NULL) {
    any_log_stream_ = true;
  }
}

int Logger::FFlush() const
{
  if(any_log_stream_) {

    int return_code = 0;
    if(use_stdout_) {
      return_code = fflush(stdout);
    }
    if(use_stderr_) {
      return_code = fflush(stderr);
    }
    if(log_file_ptr_ != NULL) {
      return_code = fflush(log_file_ptr_);
    }
    return return_code;
  }
  return 0; // no error return code consistent with fflush
}

int Logger::PrintF(const char *format, ...) const
{
  if(any_log_stream_) {

    int return_code = 0;

    // C standard macros for handling variable argument lists
    va_list argument_list;
    va_start(argument_list, 
             format); // last named parameter in the argument list, just 
                      // before the start of the variable length argument 
                      // list 

    if(use_stdout_) {
      return_code = vfprintf(stdout, format, argument_list);
    }
    if(use_stderr_) {
      return_code = vfprintf(stderr, format, argument_list);
    }
    if(log_file_ptr_ != NULL) {
      return_code = vfprintf(log_file_ptr_, format, argument_list);
    }
    va_end(argument_list);
    return return_code;
  }
  return 0; // number of characters written excluding null byte, which is
            // the return value of printf
}


} // namespace utilities
} // namespace advcomb
