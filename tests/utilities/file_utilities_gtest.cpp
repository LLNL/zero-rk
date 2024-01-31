#ifdef _WIN32
#define _USE_MATH_DEFINES //for M_PI
#endif
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "file_utilities.h"


TEST (FileIsReadable, FileDoesNotExist) 
{
  ASSERT_FALSE(zerork::utilities::FileIsReadable("data/file_does_not_exist"));  
}
TEST (FileIsReadable, FileExists) 
{
  ASSERT_TRUE(zerork::utilities::FileIsReadable("data/file_exists"));  
}

TEST (FileIsWritable, DirDoesNotExist) 
{
  ASSERT_FALSE(zerork::utilities::FileIsWritable("dir_does_not_exist/file_does_not_exist"));  
}
TEST (FileIsWritable, FileExists) 
{
  ASSERT_TRUE(zerork::utilities::FileIsWritable("data/file_exists"));  
}

TEST (GetAnyLine, Unix_NewLine)
{
  std::ifstream input_file("data/unix_test_file.txt");
  std::string line;
  int line_num = 0;
  while(zerork::utilities::GetAnyLine(input_file,&line)) {
    if(line.size() > 0) {
      EXPECT_TRUE(line == "This line ends with a new line LF");
      ++line_num;
    }
    // will read empty line for EOF
  }
  EXPECT_EQ(line_num,3);
}

TEST (GetAnyLine, Mac_CarriageReturn)
{
  std::ifstream input_file("data/mac_test_file.txt");
  std::string line;
  int line_num = 0;
  while(zerork::utilities::GetAnyLine(input_file,&line)) {
    if(line.size() > 0) {
      EXPECT_TRUE(line == "This line ends with a carriage return CR");
      ++line_num;
    }
    // will read empty line for EOF
  }
  EXPECT_EQ(line_num,4);
}

TEST (GetAnyLine, DosWindows_Return)
{
  std::ifstream input_file("data/dos_test_file.txt");
  std::string line;
  int line_num = 0;
  while(zerork::utilities::GetAnyLine(input_file,&line)) {
    if(line.size() > 0) {
      EXPECT_TRUE(line == "This line ends with a carriage return and new line CRLF");
      ++line_num;
    }
    // will read empty line for EOF
  }
  EXPECT_EQ(line_num,5);
}

TEST (GetAnyLine, HybridEOL)
{
  std::ifstream input_file("data/hybrid_test_file.txt");
  // input file formed by cat dos.line mac.line unix.line dos.line dos.line 
  //                          mac.line mac.line unix.line unix.line

  std::string DOS_LINE =
    "This line ends with a carriage return and new line CRLF";
  std::string MAC_LINE =
    "This line ends with a carriage return CR";
  std::string UNIX_LINE =
    "This line ends with a new line LF";

  int line_num=0;
  std::string line;

  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == DOS_LINE);
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == MAC_LINE);
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == UNIX_LINE);
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == DOS_LINE);
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == DOS_LINE);
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == MAC_LINE);
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == MAC_LINE);
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == UNIX_LINE);
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == UNIX_LINE);

  while(zerork::utilities::GetAnyLine(input_file,&line)) {
    if(line.size() > 0) {
      ++line_num;
    }
  }
  EXPECT_EQ(line_num,0);
}

TEST(GetFirstLineWithKeyword, CaseSensitiveNoLineNumThermoFile)
{
  std::ifstream input_file("data/thermo_test.txt");
  std::string line;

  ASSERT_TRUE(input_file);
  // without returning the line number
  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "thermo",
                                     "",   // comment_chars
                                     true, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "! Idealized thermodynamics definition for constant specific heats") << "line = " << line << std::endl; 
  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "thermo",
                                     "",   // comment_chars
                                     true, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "thermo") << "line = " << line << std::endl; 
  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "ar",
                                     "",   // comment_chars
                                     true, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "ar                fixcp ar  1               g   300.00   5000.00  1000.00      1") << "line = " << line << std::endl; 
  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "end",
                                     "",   // comment_chars
                                     true, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "end") << "line = " << line << std::endl; 
}

TEST(GetFirstLineWithKeyword, CaseSensitiveLineNumThermoFile)
{
  std::ifstream input_file("data/thermo_test.txt");
  std::string line;
  size_t total_line_num = 0;
  size_t line_num;

  ASSERT_TRUE(input_file);
  // without returning the line number
  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "thermo",
                                     "",   // comment_chars
                                     true, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "! Idealized thermodynamics definition for constant specific heats") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,1);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "thermo",
                                     "",   // comment_chars
                                     true, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "thermo") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,16);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "ar",
                                     "",   // comment_chars
                                     true, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "ar                fixcp ar  1               g   300.00   5000.00  1000.00      1") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,22);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "end",
                                     "",   // comment_chars
                                     true, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "end") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,26);

}

TEST(GetFirstLineWithKeyword, NotCaseSensitiveNoLineNumThermoFile)
{
  std::ifstream input_file("data/thermo_test.txt");
  std::string line;

  ASSERT_TRUE(input_file);
  // without returning the line number
  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "THERMO",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "! Idealized thermodynamics definition for constant specific heats") << "line = " << line << std::endl; 

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "tHeRmO",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "thermo") << "line = " << line << std::endl; 

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "AR",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "! Argon") << "line = " << line << std::endl; 

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "AR",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "ar                fixcp ar  1               g   300.00   5000.00  1000.00      1") << "line = " << line << std::endl; 

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "enD",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line);
  EXPECT_TRUE(line == "end") << "line = " << line << std::endl; 
}

TEST(GetFirstLineWithKeyword, NotCaseSensitiveLineNumThermoFile)
{
  std::ifstream input_file("data/thermo_test.txt");
  std::string line;
  size_t total_line_num = 0;
  size_t line_num;

  ASSERT_TRUE(input_file);
  // without returning the line number
  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "THERMO",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "! Idealized thermodynamics definition for constant specific heats") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,1);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "ThErMO",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "thermo") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,16);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "AR",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "! Argon") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,18);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "AR",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "ar                fixcp ar  1               g   300.00   5000.00  1000.00      1") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,22);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "enD",
                                     "",   // comment_chars
                                     false, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "end") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,26);

}
TEST(GetFirstLineWithKeyword, IgnoreComments)
{
  std::ifstream input_file("data/thermo_test.txt");
  std::string line;
  size_t total_line_num = 0;
  size_t line_num;

  ASSERT_TRUE(input_file);
  // without returning the line number
  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "THERMO",
                                     "!",   // comment_chars
                                     false, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "thermo") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,16);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "AR",
                                     "!",   // comment_chars
                                     false, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "ar                fixcp ar  1               g   300.00   5000.00  1000.00      1") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,22);

  zerork::utilities::GetFirstLineWithKeyword(input_file,
                                     "enD",
                                     "!",   // comment_chars
                                     false, // case_sensitive
                                     &line,
                                     &line_num);
  total_line_num += line_num;
  EXPECT_TRUE(line == "end") << "line = " << line << std::endl; 
  EXPECT_EQ(total_line_num,26);

}


TEST (GetAnyLine, Unix_NewLineNoEOL)
{
  std::ifstream input_file("data/unix_no_end_of_line.txt");
  std::string line;
  int line_num = 0;

  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == "This is file uses standard new line character LF '\\n',");
  zerork::utilities::GetAnyLine(input_file,&line);
  EXPECT_TRUE(line == "but the last line does not.");

  while(zerork::utilities::GetAnyLine(input_file,&line)) {
    if(line.size() > 0) {
      ++line_num;
    }
  }
  EXPECT_EQ(line_num,0);
}

TEST (LoggerConstructor, AllocateDefault) 
{
  zerork::utilities::Logger *test_log = NULL;

  test_log = new zerork::utilities::Logger();
  
  ASSERT_TRUE(test_log != NULL);

  delete test_log;
}

TEST (LoggerConstructor, AllocateStdoutOnly) 
{
  zerork::utilities::Logger *test_log = NULL;

  test_log = new zerork::utilities::Logger("",
                                   true,
                                   false);
  
  ASSERT_TRUE(test_log != NULL);

  delete test_log;
}

TEST (LoggerConstructor, AllocateStderrOnly) 
{
  zerork::utilities::Logger *test_log = NULL;

  test_log = new zerork::utilities::Logger("",
                                   false,
                                   true);
  
  ASSERT_TRUE(test_log != NULL);

  delete test_log;
}

TEST (LoggerConstructor, AllocateFileOnly) 
{
  zerork::utilities::Logger *test_log = NULL;

  test_log = new zerork::utilities::Logger("test_output/allocate_file_only.log",
                                   false,
                                   false);
  
  ASSERT_TRUE(test_log != NULL);

  delete test_log;
}

TEST (LoggerConstructor, AllocateBadDirectory) 
{
  zerork::utilities::Logger *test_log = NULL;

  test_log = new zerork::utilities::Logger("bad_directory/allocate.log",
                                   false,
                                   false);
  
  ASSERT_TRUE(test_log != NULL);

  delete test_log;
}

TEST (LoggerConstructor, AllocateAllStreams) 
{
  zerork::utilities::Logger *test_log = NULL;

  test_log = new zerork::utilities::Logger("test_output/allocate_all_streams.log",
                                   true,
                                   true);
  
  ASSERT_TRUE(test_log != NULL);

  delete test_log;
}

TEST (LoggerMessage, SimpleText)
{
  char message_format[] = "Message 1: simple text with no arguments\n";

  zerork::utilities::Logger blackhole("");
  zerork::utilities::Logger stdout_only;
  zerork::utilities::Logger stderr_only("",false,true);
  zerork::utilities::Logger file_only("test_output/file_only_message.log");
  int message_length;

  message_length = blackhole.PrintF(message_format);
  EXPECT_EQ(message_length, 0);

  message_length = stdout_only.PrintF(message_format);
  EXPECT_EQ(message_length, 41);

  message_length = stderr_only.PrintF(message_format);
  EXPECT_EQ(message_length, 41);

  message_length = file_only.PrintF(message_format);
  EXPECT_EQ(message_length, 41);

}
TEST (LoggerMessage, SingleInteger)
{
  char message_format[] = "Message 2: 7*13 = %d single integer\n";

  zerork::utilities::Logger blackhole("");
  zerork::utilities::Logger stdout_only;
  zerork::utilities::Logger stderr_only("",false,true);
  zerork::utilities::Logger file_only("test_output/file_only_message.log");
  int message_length;

  message_length = blackhole.PrintF(message_format);
  EXPECT_EQ(message_length, 0);

  message_length = stdout_only.PrintF(message_format,7*13);
  EXPECT_EQ(message_length, 36);

  message_length = stderr_only.PrintF(message_format,7*13);
  EXPECT_EQ(message_length, 36);

  message_length = file_only.PrintF(message_format,7*13);
  EXPECT_EQ(message_length, 36);

}

TEST (LoggerMessage, MultipleValues)
{
  char message_format[] = "Message 3: old pi = %d, %5.3f from %s\n";

  zerork::utilities::Logger blackhole("");
  zerork::utilities::Logger stdout_only;
  zerork::utilities::Logger stderr_only("",false,true);
  zerork::utilities::Logger file_only("test_output/file_only_message.log");
  int message_length;

  message_length = blackhole.PrintF(message_format);
  EXPECT_EQ(message_length, 0);

  message_length = stdout_only.PrintF(message_format,
                                      (int)M_PI,
                                      M_PI,
                                      "M_PI");
  EXPECT_EQ(message_length, 39);

  message_length = stderr_only.PrintF(message_format,
                                      (int)M_PI,
                                      M_PI,
                                      "M_PI");
  EXPECT_EQ(message_length, 39);

  message_length = file_only.PrintF(message_format,
                                    (int)M_PI,
                                    M_PI,
                                    "M_PI");
  EXPECT_EQ(message_length, 39);

}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
