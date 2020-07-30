#ifndef THERMO_PARSER_H_
#define THERMO_PARSER_H_

#include <string>
#include <vector>

#include "zerork/elemental_composition.h"

#include "janaf_thermo.h"


const char WHITESPACE[] = " \t";

enum ThermoLineType {COMMENT, BLANK, THERMO_START, THERMO_RANGE, THERMO_END,
                     JANAF_DATA, SPECIES_DATA, INVALID};

int ParseThermoFile(const char *file_name,
                    std::vector<std::string> &species_name,
                    std::vector<std::string> &species_definition,
                    std::vector<zerork::ElementalComposition> &species_composition,
                    std::vector<JanafThermoData> &janaf_data,
                    std::vector<double> &global_range);
int ParseMechFileForSpecies(const char *file_name,
                            std::vector<std::string> &species_name);

void ProcessCompString(const std::string comp_str,
                       zerork::ElementalComposition &species_composition);

bool StartsWithKey(const std::string &key, 
                   const size_t start_position, 
                   const std::string &str,
                   const bool is_case_sensitive);

ThermoLineType GetThermoLineType(const int current_line,
                                 const int thermo_start_line,
                                 const int janaf_line_count,
                                 std::string &line);
bool SetJanafLineData(const std::string &str,
                      const int janaf_line,
                      JanafThermoData &data);
bool SetGlobalRange(const std::string &str,
                    std::vector<double> &global_range);
bool SetSpeciesLineData(const int current_line,
                        const std::string &str,
                        const std::vector<double> &global_range,
                        std::string &species_name,
                        std::string &species_definition,
                        zerork::ElementalComposition &species_composition,
                        JanafThermoData &data);
size_t SafeSubstring(const std::string &source_str,
                     const size_t start_pos,
                     const size_t length,
                     std::string &extract_str);

size_t ExtractFirstNonWhitespace(const std::string &source_str,
                                 std::string &extract_str);
size_t ExtractFirstNonWhitespace(std::string &str);

size_t SplitWords(const std::string &source_str,
                  const std::string &split_chars,
                  const std::string &stop_chars,
                  std::vector<std::string> &words);

#endif
