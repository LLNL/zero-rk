#include "stdlib.h" // exit()
#include "stdio.h"

#include <map>
#include <fstream>
#include <iostream>

#include "thermo_parser.h"

#include "utilities/file_utilities.h" //zerork::utiltiies::GetAnyLine

int ParseThermoFile(const char *file_name,
                    std::vector<std::string> &species_name,
                    std::vector<std::string> &species_definition,
                    std::vector<zerork::ElementalComposition> &species_composition,
                    std::vector<JanafThermoData> &janaf_data,
                    std::vector<double> &global_range)
{
  using zerork::utilities::GetAnyLine;
  bool found_end_key=false;
  ThermoLineType line_type;
  int num_species_records = 0;
  int num_lines_read   = 0;
  int janaf_line_count = 3;   // 3 janaf lines are required to be read before
                              // loading the next species
  int thermo_start_line = -1; // line number where THERMO keyword is found

  JanafThermoData current_janaf_data;
  std::map<std::string, int> species_position;
  std::string current_species_name;
  std::string current_species_definition;
  zerork::ElementalComposition current_species_composition;
  

  // TODO: make file check with C++ fstream  
  FILE *fptr = fopen(file_name,"r");
  if(fptr == NULL) {
    printf("ERROR: In ParseThermoFile(...),\n");
    printf("       could not open thermo file %s for read\n",file_name);
    fflush(stdout);
    return -1;
  }
  fclose(fptr);  
  // --------------------------------------
  std::ifstream thermo_file(file_name, std::ios_base::in);
  std::string line;

  species_name.clear();
  species_position.clear();
  species_definition.clear();
  species_composition.clear();
  janaf_data.clear();

  while(GetAnyLine(thermo_file,&line) && found_end_key == false) {
    ++num_lines_read;
    //printf("num_lines_read = %d\n",num_lines_read); fflush(stdout);
    line_type = GetThermoLineType(num_lines_read,
                                  thermo_start_line,
                                  janaf_line_count,
                                  line);
    switch(line_type) {
    case THERMO_START:
      thermo_start_line=num_lines_read;
      break;
    case THERMO_END:
      found_end_key = true;
      break;
    case THERMO_RANGE:
      // extract global range
      if(!SetGlobalRange(line,global_range)) {
        printf("Failed to global T range at line %d\n",num_lines_read);
        exit(-1);
      }
      break;
    case SPECIES_DATA:
      // extract species name, species definition, Tlow, Thigh, Tmatch
      if(!SetSpeciesLineData(num_lines_read,
                             line,
                             global_range,
                             current_species_name,
                             current_species_definition,
                             current_species_composition,
                             current_janaf_data)) {
        printf("Failed to set new species data at line %d\n",num_lines_read);
        exit(-1);
      }                      
      janaf_line_count = 0; // reset janaf line count
      break;
    case JANAF_DATA:
      ++janaf_line_count;
      // extract janaf line
      if(!SetJanafLineData(line,janaf_line_count,current_janaf_data)) {
        printf("Failed to set JANAF data at line %d\n",num_lines_read);
        exit(-1);
      }
      break;
    case INVALID:
      printf("Stopping at INVALID line %d\n",num_lines_read);
      exit(-1);
    default:
      break;
    }
    if(line_type == JANAF_DATA && janaf_line_count == 3) {
      // TODO: add logic to keep the first or overwrite the repeated species
      // record species record if new
      if(species_position.find(current_species_name) == 
         species_position.end()) {
        // store the species name, definition and JANAF data to the vectors
        species_name.push_back(current_species_name);
        species_definition.push_back(current_species_definition);
        species_composition.push_back(current_species_composition);
        janaf_data.push_back(current_janaf_data);
        // record the species position into the map
        species_position[current_species_name] = num_species_records;
      } else {
        printf("INFO: Ignoring duplicate thermo for species %s\n",
               current_species_name.c_str());
        printf("      at lines [%d - %d]\n",
               num_lines_read-3,num_lines_read);
      } 
    }
    //printf("read %d lines\n",num_lines_read); fflush(stdout);
  }
  return species_position.size();
}

int ParseMechFileForSpecies(const char *file_name,
                            std::vector<std::string> &species_name)
{
  using zerork::utilities::GetAnyLine;
  size_t start_position;
  int num_lines_read = 0;
  bool found_species_key = false;
  bool found_end_key = false;
  std::map<std::string, int> species_position; // map used to locate
                                               // duplicate entries

  // TODO: make file check with C++ fstream  
  FILE *fptr = fopen(file_name,"r");
  if(fptr == NULL) {
    printf("ERROR: In ParseMechFileForSpecies(...),\n");
    printf("       could not open mech file %s for read\n",file_name);
    fflush(stdout);
    return -1;
  }
  fclose(fptr);  
  // --------------------------------------
  std::ifstream thermo_file(file_name, std::ios_base::in);
  std::string line;
  std::vector<std::string> line_words;

  species_name.clear();
  species_position.clear();

  while(GetAnyLine(thermo_file,&line) && found_end_key == false) {
    ++num_lines_read;
 
    // find the first character that isn't whitespace or a return
    start_position = line.find_first_not_of(" \t\n\r");
    // returns std::string::npos if nothing is found

    if(start_position != std::string::npos) {
      // a character that isn't whitespace or a return is found

      if(line[start_position] != '!' && line[start_position] != '#') {
        // not a comment line
        if(StartsWithKey("SPECIES",start_position,line,false)) {
          found_species_key = true;
        } else if((StartsWithKey("END",start_position,line,false) ||
                  StartsWithKey("REAC",start_position,line,false) ||
                  StartsWithKey("ELEM",start_position,line,false) || 
                  StartsWithKey("THERM",start_position,line,false)) &&
                  found_species_key) {
          // determine the end of the species section by the "END" keyword
          // or one of the other section keywords 
          found_end_key = true;
        } else if(found_species_key && !found_end_key) {
          // line containing species declarations
          SplitWords(line,
                     " \t",   // characters marking the splits between words
                     "\n\r!", // character marking the end of the line to split
                     line_words);

          //printf("line number %d: %lu species\n",
          //       num_lines_read,line_words.size());
          for(size_t j=0; j<line_words.size(); ++j) {
            //printf("  %s\n",line_words[j].c_str());
            if(species_position.find(line_words[j]) == 
               species_position.end()) {
              // if species name is not found in position map, it is a new
              // species to record
              species_position[line_words[j]] = species_name.size();
              species_name.push_back(line_words[j]);

            } else {
              // found duplicate species
              printf("INFO: Ignoring duplicate thermo for species %s\n",
                     line_words[j].c_str());
              printf("      at line %d.\n", num_lines_read);
            }     
          }        
        }
      }
    }
  }
  return species_position.size();
}

ThermoLineType GetThermoLineType(const int current_line,
                                 const int thermo_start_line,
                                 const int janaf_line_count,
                                 std::string &line)
{
  bool is_column80_one=false;
  size_t found_position;

  // check if column 80 is one, will be used to distinguish a number starting
  // line as the start of a new species instead of a 5th thermo line
  if(line.size() >= 80) {
    if(line[79] == '1') {
      is_column80_one = true;
    }
  }

  // find the first character that isn't whitespace or a return
  found_position = line.find_first_not_of(" \t\n\r");
 
  if(found_position == std::string::npos) {
    return BLANK; // no letter, number or non-whitespace or return character
  }
  if(line[found_position] == '!' || line[found_position] == '#') {
    return COMMENT;
  }
  if((line[found_position] >= '0' && line[found_position] <= '9') ||
     line[found_position] == '-' || line[found_position] == '.') {

    // found a number first, line is either the global temperature range
    // or a JANAF data line, or a species that starts with a number
    if(current_line == thermo_start_line+1 && janaf_line_count == 3) {
      return THERMO_RANGE;
    } else if(janaf_line_count < 3) {
      return JANAF_DATA;
    } else if(janaf_line_count == 3 && is_column80_one) {
      return SPECIES_DATA; // note that species 
    }
    printf("ERROR: can not identify line %d.\n",
           current_line);
    printf("       It is does not seem to be a thermo range, species, \n");
    printf("       or a valid JANAF data line\n");
    return INVALID;
    
  } else if((line[found_position] >= 'a' && line[found_position] <= 'z') ||
            (line[found_position] >= 'A' && line[found_position] <= 'Z') ||
            line[found_position] == '(' || line[found_position] == '[' ||
            line[found_position] == '{' ) {

    // check for specific keys - 'false' checks for a case insensitive match
    if(StartsWithKey("THERMO",found_position,line,false)) { 
      return THERMO_START;
    }
    if(StartsWithKey("END",found_position,line,false)) {
      return THERMO_END;
    }
    // make sure it is a new species that isn't declared in the middle
    // of the janaf
    if(janaf_line_count == 3) {
      return SPECIES_DATA;
    }
    printf("ERROR: can not specify a new species in the middle of the\n");
    printf("       JANAF data at line %d\n",current_line);
    return INVALID;
  }
  printf("ERROR: can not identify line %d with starting character %c\n",
         current_line,line[found_position]);

  return INVALID;
}

bool StartsWithKey(const std::string &key, 
                   const size_t start_position, 
                   const std::string &str,
                   const bool is_case_sensitive)
{
  size_t key_length = key.size();
  std::string key_copy = key;
  //std::string substring_copy = str.substr(start_position,key_length);
  std::string substring_copy;
  SafeSubstring(str,
                start_position,
                key_length,
                substring_copy);
  
  // DEBUG
  //printf("key = %s\n",key.c_str());
  //printf("str = %s\n",str.c_str());
  //fflush(stdout);

  if(substring_copy.size() < key_copy.size()) {
    return false;
  }

  if(!is_case_sensitive) {
    for(size_t j=0; j<key_length; ++j) {
      if('A' <= key_copy[j] && key_copy[j] <= 'Z') {
        key_copy[j] += 'a' - 'A';
      }
      if('A' <= substring_copy[j] && substring_copy[j] <= 'Z') {
        substring_copy[j] += 'a' - 'A';
      }
    }
  }

  if(key_copy.compare(substring_copy) == 0) {
    return true;
  }
  return false;
}
// return true if successful and false if a zero-length field is found
bool SetJanafLineData(const std::string &str,
                      const int janaf_line,
                      JanafThermoData &data)
{
  const int num_width = 15;
  int num_fields = 5;
  double read_coef[5];
  std::string substring;

  if(janaf_line == 3) { // last line only needs 4 fields
    num_fields = 4;
  } 

  for(int j=0; j<num_fields; ++j){
    //substring = str.substr(num_width*j,num_width);
    SafeSubstring(str,
                  num_width*j,
                  num_width,
                  substring);
    if(substring.size() < 1) {
      printf("ERROR: In SetJanafLineData(),\n");
      printf("       could not find JANAF coefficient field %d in set %d\n",
             j,janaf_line);
      return false;
    }
    read_coef[j] = atof(substring.c_str());
  }
  // record data to the coefficients  
  if(janaf_line == 1) {
    // set the first five high coefficients
    data.high_coef[0] = read_coef[0];
    data.high_coef[1] = read_coef[1];
    data.high_coef[2] = read_coef[2];
    data.high_coef[3] = read_coef[3];
    data.high_coef[4] = read_coef[4];
  } else if(janaf_line == 2) {
    // set the last two high coefficients and the first three low coefficients
    data.high_coef[5] = read_coef[0];
    data.high_coef[6] = read_coef[1];
    data.low_coef[0]  = read_coef[2];
    data.low_coef[1]  = read_coef[3];
    data.low_coef[2]  = read_coef[4];
  } else if(janaf_line == 3) {
    // set the last 4 low coefficients
    data.low_coef[3] = read_coef[0];
    data.low_coef[4] = read_coef[1];
    data.low_coef[5] = read_coef[2];
    data.low_coef[6] = read_coef[3];
  } else {
    printf("ERROR: In SetJanafLineData(),\n");
    printf("       could not process JANAF coefficient set %d\n",
           janaf_line);
     return false;
  }

  return true;
}
// return true if successful and false if a zero-length field is found
bool SetGlobalRange(const std::string &str,
                    std::vector<double> &global_range)
{
  const double initial_value=-1.0e300;
  int num_filled;
  global_range.clear();
  global_range.assign(3,initial_value);
  // TODO: chemkin 2 standard is 3 fixed 10 character fields, decide if we
  //       need this check
  num_filled = sscanf(str.c_str(),
                      "%lf%lf%lf",
                      &global_range[0],
                      &global_range[1],
                      &global_range[2]);
  if(num_filled == 3) {
    return true;
  }
  return false;
}

bool SetSpeciesLineData(const int current_line,
                        const std::string &str,
                        const std::vector<double> &global_range,
                        std::string &species_name,
                        std::string &species_definition,
                        zerork::ElementalComposition &species_composition,
                        JanafThermoData &data)
{
  std::string T_min_string,T_match_string,T_max_string, comp_str;
  //size_t start_pos, end_pos,length;
  // species name in columns 1:18 this seems to be the original standard
  // but no longer strictly used as comments are allowed to fill this space
  // so it appears as the species name is the first string token in the first
  // 18 characters
  //species_name = str.substr(0,18);
  SafeSubstring(str,
                0,
                18,
                species_name);
  // remove any whitespace padding
  ExtractFirstNonWhitespace(species_name);
  if(species_name.size() < 1) {
    printf("ERROR: could not extract species name on line %d\n",
           current_line);
    return false;
  }

  // store the full species definition for later use columns 1:45
  //species_definition = str.substr(0,45);
  SafeSubstring(str,
                0,
                45,
                species_definition);
  if(species_definition.size() != 45) {
    printf("WARNING: line %d - species definition string has unexpected\n",
           current_line);
    printf("         length %lu != d\n",species_definition.size());
  }
  SafeSubstring(str,
                45,
                10,
                T_min_string);
  ExtractFirstNonWhitespace(T_min_string);
  SafeSubstring(str,
                65,
                8,
                T_match_string);
  ExtractFirstNonWhitespace(T_match_string);
  SafeSubstring(str,
                55,
                10,
                T_max_string);
  ExtractFirstNonWhitespace(T_max_string);

  if(T_min_string.size() < 1) {
    printf("WARNING: Using global T_min = %10.5f for species %s on line %d\n",
           global_range[0],
           species_name.c_str(),
           current_line);
    data.T_min = global_range[0];
  } else {
    data.T_min = atof(T_min_string.c_str());
  }

  if(T_match_string.size() < 1) {
    printf("WARNING: Using global T_match = %10.5f for species %s on line %d\n",
           global_range[1],
           species_name.c_str(),
           current_line);
    data.T_match = global_range[1];
  } else {
    data.T_match = atof(T_match_string.c_str());
  }

  if(T_max_string.size() < 1) {
    printf("WARNING: Using global T_max = %10.5f for species %s on line %d\n",
           global_range[1],
           species_name.c_str(),
           current_line);
    data.T_max = global_range[2];
  } else {
    data.T_max = atof(T_max_string.c_str());
  }

  //printf("Name = %s (string length = %lu)\n",
  //       species_name.c_str(),species_name.size());

  species_composition.clear();
  species_composition.name() = species_name;
  SafeSubstring(str, 24, 5, comp_str);
  ProcessCompString(comp_str,species_composition);

  SafeSubstring(str, 29, 5, comp_str);
  ProcessCompString(comp_str,species_composition);

  SafeSubstring(str, 34, 5, comp_str);
  ProcessCompString(comp_str,species_composition);

  SafeSubstring(str, 39, 5, comp_str);
  ProcessCompString(comp_str,species_composition);

  return true;
}

void ProcessCompString(const std::string comp_str,
                       zerork::ElementalComposition &species_composition)
{
  std::vector<std::string> comp_vec;
  SplitWords(comp_str, " ", "#!", comp_vec);
  if(comp_vec.size() == 2) {
    std::string elem = comp_vec[0];
    int n_atoms = atoi(comp_vec[1].c_str());
    species_composition[elem] += n_atoms;
  }
}

size_t SafeSubstring(const std::string &source_str,
                     const size_t start_pos,
                     const size_t length,
                     std::string &extract_str)
{
  extract_str.clear();
  if(start_pos >= source_str.size()) {
    // will throw an exemption if substr is called
    return 0; 
  }
  extract_str =source_str.substr(start_pos,length);  

  return extract_str.size();
}

size_t ExtractFirstNonWhitespace(const std::string &source_str,
                                 std::string &extract_str)
{
  size_t start_pos, end_pos;
  extract_str.clear();
  if(source_str.size() < 1) {
    return 0;
  }
  start_pos = source_str.find_first_not_of(WHITESPACE);
  if(start_pos == std::string::npos) {
    // did not find any non-whitespace
    return 0;
  }
  end_pos   = source_str.find_first_of(WHITESPACE,
                                       start_pos); // string position to 
                                                   // start the next search
  if(end_pos == std::string::npos) {
    // using length equal to string::npos extracts to the end of the string
    return SafeSubstring(source_str,
                         start_pos,
                         end_pos,
                         extract_str);
  } else {
    return SafeSubstring(source_str,
                         start_pos,
                         end_pos-start_pos,
                         extract_str);
  }
}
size_t ExtractFirstNonWhitespace(std::string &str)
{
  std::string str_copy = str;
  return ExtractFirstNonWhitespace(str_copy,
                                   str);
}

// SplitWords
//
// Inputs: source_str (const) source string to split into words
//         split_chars (const) string containing all the characters defining
//           a split between words (e.g. " \t," would treat spaces, tabs and
//           commas as delimiters)
//         stop_chars (const) string containing all the characters that
//           determine when the string should stop being parsed
//           (e.g. "!#\n\r" would stop parsing for comment characters '!' 
//           and '#', and the new line/carriage return characters).
//
// Outputs: words is a string vector that contains the "words" found in the
//          source string.  A "word" is composed of any character not found
//          in split_chars and stop_chars. Note that the words vector is
//          cleared first so the words are not appended to an existing vector.
//
// Notes on using the string class find_first* functions:
//
// If the find_first* commands fail to find a character in the string that
// matches the conditions, they return std::string::npos.
//
// npos is a static member constant value with the greatest possible value
// for an element of type size_t.  This value, when used as the value for a 
// len (or sublen) parameter in string's member functions, means "until the 
// end of the string".
// 
// As a return value, it is usually used to indicate no matches.
// This constant is defined with a value of -1, which because size_t is an
// unsigned integral type, it is the largest possible representable value
// for this type.
size_t SplitWords(const std::string &source_str,
                  const std::string &split_chars,
                  const std::string &stop_chars,
                  std::vector<std::string> &words)
{
  words.clear();
  if(source_str.size() == 0) {
    return 0;
  }

  std::string next_word;

  size_t start_pos = source_str.find_first_not_of(split_chars); 
  size_t next_pos =  source_str.find_first_of(split_chars,start_pos);
  size_t end_pos   = source_str.find_first_of(stop_chars);

  // source_str[start_pos] to source_str[next_pos-1] represents the portion of
  // the source string to split into the next word
  while(start_pos < next_pos && start_pos < end_pos) {
    SafeSubstring(source_str,
                  start_pos,
                  next_pos-start_pos,
                  next_word);
    words.push_back(next_word);
    start_pos = source_str.find_first_not_of(split_chars,next_pos);
    next_pos  = source_str.find_first_of(split_chars,start_pos);
  }

  return words.size();
}
