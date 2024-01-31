#include <math.h>

#include <iostream>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "string_utilities.h"

const char REACTION_COMMENT[] = "h2 + m = h + h + m ! Basic hydrogen reaction";
const char C_COMMENT[] = "  ++increment; // adds one to the increment\n";
const char ALL_COMMENT[] = "# this is only a comment\n";
const char REMOVE_SEPARATORS[] = "// space \t tab \r carriage_return \n new_line \f form_feed \v whatever_slash_v_is , plus commas!";
const char REACTION_TOKENS[] = "o2+h+(m)=> ho2 + (m) ! try to split the species";
const char DATA_LINE[] = " 1.23  4.56  7.89\t0.1234\n";
const char PANGRAM_TEST[] = "tHE qUiCk brOwN fox JuMps oVer the LaZy dog.";

const char THERMO1[] = "h                 120186h   1               g  0300.00   5000.00  1000.00      1";
const char THERMO2[] = " 0.02500000e+02 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00    2";
const char THERMO3[] = " 0.02547163e+06-0.04601176e+01 0.02500000e+02 0.00000000e+00 0.00000000e+00    3";
const char THERMO4[] = " 0.00000000e+00 0.00000000e+00 0.02547163e+06-0.04601176e+01                   4";

TEST(GetUpperCase, PangramTest)
{
  std::string upper_case = zerork::utilities::GetUpperCase(PANGRAM_TEST);
  std::string expected = "THE QUICK BROWN FOX JUMPS OVER THE LAZY DOG.";
  EXPECT_TRUE(expected == upper_case) << "expected = " << expected
				      << ", upper_case = " << upper_case
                                      << std::endl;
}

TEST(GetLowerCase, PangramTest)
{
  std::string lower_case = zerork::utilities::GetLowerCase(PANGRAM_TEST);
  std::string expected = "the quick brown fox jumps over the lazy dog.";
  EXPECT_TRUE(expected == lower_case) << "expected = " << expected
				      << ", lower_case = " << lower_case
                                      << std::endl;
}

TEST(GetCharCount, ReactionComment)
{
  EXPECT_EQ(zerork::utilities::GetCharCount(REACTION_COMMENT,'h'), 4)
    << "number of 'h' in '" << REACTION_COMMENT << "'" << std::endl;

}
TEST(GetCharCount, PangramTest)
{
  EXPECT_EQ(zerork::utilities::GetCharCount(PANGRAM_TEST,' '), 8) 
    << "number of spaces in '" << PANGRAM_TEST << "'" << std::endl;
}

TEST(StringIsInt, PlusPrefix)
{
  EXPECT_TRUE(zerork::utilities::StringIsInt("+12345"));
}
TEST(StringIsInt, NegativePrefix)
{
  EXPECT_TRUE(zerork::utilities::StringIsInt("-67890"));
}
TEST(StringIsInt, MaxInt)
{
  EXPECT_TRUE(zerork::utilities::StringIsInt("2147483647"));
}
TEST(StringIsInt, MinInt)
{
  EXPECT_TRUE(zerork::utilities::StringIsInt("-2147483648"));
}
TEST(StringIsInt, MaxIntPlusOne)
{
  EXPECT_FALSE(zerork::utilities::StringIsInt("2147483648"));
}
TEST(StringIsInt, MinIntMinusOne)
{
  EXPECT_FALSE(zerork::utilities::StringIsInt("-2147483649"));
}
TEST(StringIsInt, ManyZeros)
{
  EXPECT_TRUE(zerork::utilities::StringIsInt("+0000000000001234567890"));
  // verify that the above string actually converts to the expected value
  EXPECT_EQ(1234567890, atof("+0000000000001234567890"));
}

TEST(StringIsInt, JustPlus)
{
  EXPECT_FALSE(zerork::utilities::StringIsInt("+"));
}

TEST(StringIsInt, JustMinus)
{
  EXPECT_FALSE(zerork::utilities::StringIsInt("-"));
}

TEST(StringIsDouble, ValidBoltzmann)
{
  double Boltzmann = 1.38064852e-23;
  EXPECT_TRUE(zerork::utilities::StringIsDouble("1.38064852e-23"));
  EXPECT_TRUE(zerork::utilities::StringIsDouble("138064852.e-31"));
  EXPECT_TRUE(zerork::utilities::StringIsDouble("138064852e-31"));
  EXPECT_TRUE(zerork::utilities::StringIsDouble(".138064852e-22"));
  EXPECT_TRUE(zerork::utilities::StringIsDouble("+1.38064852E-23"));
  EXPECT_TRUE(zerork::utilities::StringIsDouble("+1.38064852E-023"));

  EXPECT_TRUE(zerork::utilities::StringIsDouble("+1.38064852E-00023"));
  // verify that the above string actually converts to the expected value
  EXPECT_EQ(Boltzmann, atof("+1.38064852E-00023"));

  EXPECT_TRUE(
    zerork::utilities::StringIsDouble("0.0000000000000000000000138064852e+000"));
  // verify that the above string actually converts to the expected value
  EXPECT_EQ(Boltzmann,
            atof("0.0000000000000000000000138064852e+000"));
}
TEST(StringIsDouble, InvalidBoltzmann)
{
  // the following are all expected to be invalid doubles
  EXPECT_FALSE(zerork::utilities::StringIsDouble("0x38064852"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+1.38064852D23"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+1.38064852d23"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+1.38064852 23"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+1.38064852E 23"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+1.38064852E--23"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+1.38064852.E-23"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+1.38064852eE-23"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("++1.38064852E-23"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+-1.38064852E-23"));
}

TEST(StringIsDouble, AlmostOverflow)
{
  EXPECT_TRUE(zerork::utilities::StringIsDouble("1.0e+308"))
    << "atof = " << strtod("1.0e+308",NULL) << std::endl;
  EXPECT_TRUE(zerork::utilities::StringIsDouble("-1.0e+308"))
    << "atof = " << strtod("-1.0e+308",NULL) << std::endl;
}

TEST(StringIsDouble, Overflow)
{
  EXPECT_FALSE(zerork::utilities::StringIsDouble("2.0e+308"))
    << "atof = " << strtod("2.0e+308",NULL) << std::endl;
  EXPECT_FALSE(zerork::utilities::StringIsDouble("-2.0e+308"))
    << "atof = " << strtod("-2.0e+308",NULL) << std::endl;
}

TEST(StringIsDouble, AlmostUnderflow)
{
  EXPECT_TRUE(zerork::utilities::StringIsDouble("3.0e-308"))
    << "strtod = " << strtod("3.0e-308",NULL) << std::endl;
  EXPECT_TRUE(zerork::utilities::StringIsDouble("-3.0e-308"))
    << "strtod = " << strtod("-3.0e-308",NULL) << std::endl;
}
#ifdef _WIN32
//For some unknown reason 'strod' in Visual Studio doesn't
//underflow at DBL_MIN and we have to push this further to
//get it to set errno. Leaving other systems at the "normal"
//value.
TEST(StringIsDouble, Underflow)
{
  EXPECT_FALSE(zerork::utilities::StringIsDouble("2.0e-324"))
    << "strtod = " << strtod("2.0e-324",NULL) << std::endl;
  EXPECT_FALSE(zerork::utilities::StringIsDouble("-2.0e-324"))
    << "strtod = " << strtod("-2.0e-324",NULL) << std::endl;
}
#else
TEST(StringIsDouble, Underflow)
{
  EXPECT_FALSE(zerork::utilities::StringIsDouble("2.0e-308"))
    << "strtod = " << strtod("2.0e-308",NULL) << std::endl;
  EXPECT_FALSE(zerork::utilities::StringIsDouble("-2.0e-308"))
    << "strtod = " << strtod("-2.0e-308",NULL) << std::endl;
}
#endif

TEST(StringIsDouble, JustPlus)
{
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+"));
}

TEST(StringIsDouble, JustMinus)
{
  EXPECT_FALSE(zerork::utilities::StringIsDouble("-"));
}

TEST(StringIsDouble, JustExponent)
{
  EXPECT_FALSE(zerork::utilities::StringIsDouble("e+002"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("e-004"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("+e3"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("-e-6"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("e5"));
  EXPECT_FALSE(zerork::utilities::StringIsDouble("e-10"));
}



TEST(SplitBeforeDelimiter, SingleCharacter) 
{
  std::string expected = "h2 + m = h + h + m ";
  std::string extracted;
  
  size_t extracted_len = 
    zerork::utilities::SplitBeforeDelimiter(std::string(REACTION_COMMENT),
                                    std::string("!"),
                                    false,
                                    &extracted);
  EXPECT_EQ(extracted_len,19);
  EXPECT_TRUE(expected == extracted) << "expected = " << expected
				     << ", extracted = " << extracted
                                     << std::endl;
}


TEST(SplitBeforeDelimiter, SplitToSelf) 
{
  std::string expected = "h2 + m = h + h + m ";
  std::string extracted = std::string(REACTION_COMMENT);
  
  size_t extracted_len = 
    zerork::utilities::SplitBeforeDelimiter(extracted,
                                    std::string("!"),
                                    false,
                                    &extracted);
  EXPECT_EQ(extracted_len,19);
  EXPECT_TRUE(expected == extracted) << "expected = " << expected
				     << ", extracted = " << extracted
                                     << std::endl;
}

TEST(SplitBeforeDelimiter, MultiCharacter) 
{
  std::string expected = "  ++increment; ";
  std::string extracted;
  
  size_t extracted_len = 
    zerork::utilities::SplitBeforeDelimiter(std::string(C_COMMENT),
                                    std::string("//"),
                                    true,
                                    &extracted);
  EXPECT_EQ(extracted_len,15);
  EXPECT_TRUE(expected == extracted) << "expected = " << expected
				     << ", extracted = " << extracted
                                     << std::endl;
}

TEST(SplitBeforeDelimiter, AllComment) 
{
  std::string expected = "";
  std::string extracted;
  
  size_t extracted_len = 
    zerork::utilities::SplitBeforeDelimiter(std::string(ALL_COMMENT),
                                    std::string("!#"),
                                    false,
                                    &extracted);
  EXPECT_EQ(extracted_len,0);
  EXPECT_TRUE(expected == extracted) << "expected = " << expected
				     << ", extracted = " << extracted
                                     << std::endl;
}

TEST(SplitStringToVector, DataLine)
{
  std::vector<std::string> token_vector;
  size_t num_tokens =
    zerork::utilities::SplitStringToVector(DATA_LINE,
                                   zerork::utilities::WHITESPACE,
                                   &token_vector);

  ASSERT_EQ(num_tokens, 4) << "input = " << DATA_LINE << std::endl;
  // check strings
  EXPECT_TRUE(token_vector[0] == "1.23") << "token_vector[0] = " 
                                         <<  token_vector[0] << std::endl;
  EXPECT_TRUE(token_vector[1] == "4.56") << "token_vector[1] = " 
                                         <<  token_vector[1] << std::endl;
  EXPECT_TRUE(token_vector[2] == "7.89") << "token_vector[2] = " 
                                         <<  token_vector[2] << std::endl;
  EXPECT_TRUE(token_vector[3] == "0.1234") << "token_vector[3] = " 
                                         <<  token_vector[3] << std::endl;
  // check converted strings
  EXPECT_EQ(atof(token_vector[0].c_str()), 1.23);
  EXPECT_EQ(atof(token_vector[1].c_str()), 4.56);
  EXPECT_EQ(atof(token_vector[2].c_str()), 7.89);
  EXPECT_EQ(atof(token_vector[3].c_str()), 0.1234);
}

TEST(SplitStringToVector, ReactantProductTokens)
{
  std::string reaction;
  std::vector<std::string> token_vector;
  zerork::utilities::SplitBeforeDelimiter(std::string(REACTION_TOKENS),
                                  "!",
                                  false,
                                  &reaction);
  size_t num_tokens =
    zerork::utilities::SplitStringToVector(reaction,
                                   "<=>",
                                   &token_vector);
  ASSERT_EQ(num_tokens, 2);
  EXPECT_TRUE(token_vector[0] == "o2+h+(m)") << "reactant (left side)"
                                             << std::endl;
  EXPECT_TRUE(token_vector[1] == " ho2 + (m) ") << "reactant (right side)"
                                             << std::endl;
}

TEST(SplitStringToVector, SpeciesTokens)
{
  std::string reaction;
  std::vector<std::string> token_vector;
  zerork::utilities::SplitBeforeDelimiter(std::string(REACTION_TOKENS),
                                  "!",
                                  false,
                                  &reaction);
  size_t num_tokens =
    zerork::utilities::SplitStringToVector(reaction,
                                   "+<=>",
                                   &token_vector);
  for(size_t j=0; j<num_tokens; ++j) {
    // A copy is made when the input and output strings are the same object
    // address
    zerork::utilities::EraseAllCharacters(token_vector[j],
                                  zerork::utilities::WHITESPACE,
                                  &token_vector[j]);
  }

  ASSERT_EQ(num_tokens, 5);
  EXPECT_TRUE(token_vector[0] == "o2")  << "1st species"
                                        << std::endl;
  EXPECT_TRUE(token_vector[1] == "h")   << "2nd species"
                                        << std::endl;
  EXPECT_TRUE(token_vector[2] == "(m)") << "3rd species"
                                        << std::endl;
  EXPECT_TRUE(token_vector[3] == "ho2") << "4th species"
                                        << std::endl;
  EXPECT_TRUE(token_vector[4] == "(m)") << "5th species"
                                        << std::endl;
}


TEST(SplitStringToVector, ValidThermo)
{
  std::vector<std::string> split_string;
  
  size_t start_position[] = { 0, 15, 30, 45, 60, 79};
  size_t stop_position[]  = {15, 30, 45, 60, 75, 80};
  std::vector<size_t> start_vector(start_position,
                                   start_position+
                                   sizeof(start_position)/sizeof(size_t));
  std::vector<size_t> stop_vector(stop_position,
                                  stop_position+
                                  sizeof(stop_position)/sizeof(size_t));

  size_t num_fields = 
    zerork::utilities::SplitStringToVector(THERMO2,
                                   start_vector,
                                   stop_vector,
                                   &split_string);
  ASSERT_EQ(num_fields, 6);
  EXPECT_EQ(atof(split_string[0].c_str()), 2.5) << "split_string[0] = "
					<< split_string[0] << std::endl;
  EXPECT_EQ(atof(split_string[1].c_str()), 0.0) << "split_string[1] = "
					<< split_string[1] << std::endl;
  EXPECT_EQ(atof(split_string[2].c_str()), 0.0) << "split_string[2] = "
					<< split_string[2] << std::endl;
  EXPECT_EQ(atof(split_string[3].c_str()), 0.0) << "split_string[3] = "
					<< split_string[3] << std::endl;
  EXPECT_EQ(atof(split_string[4].c_str()), 0.0) << "split_string[4] = "
					<< split_string[4] << std::endl;
  EXPECT_EQ(atoi(split_string[5].c_str()), 2) << "split_string[5] = "
					<< split_string[5] << std::endl;
}


TEST(EraseAllCharacters, RemoveAllSeparators) 
{
  std::string expected = "//spacetabcarriage_returnnew_lineform_feedwhatever_slash_v_ispluscommas!";
  std::string extracted;

  zerork::utilities::EraseAllCharacters(REMOVE_SEPARATORS,
                                std::string(zerork::utilities::WHITESPACE)+",",
                                &extracted);
  EXPECT_TRUE(expected == extracted) << "expected = " << expected
				     << ", extracted = " << extracted
                                     << std::endl;
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
