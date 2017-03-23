#ifndef STRING_OPS_H_
#define STRING_OPS_H_

#include <string>
#include <vector>

void split_by_delim(const std::string &s, char delim, 
		    std::vector<std::string>& substrings);

std::string uppercase(std::string str);

bool string_starts_with(std::string& s, std::string prefix);

bool string_ends_with(std::string& s, std::string suffix); 

bool orderByLengthAndSequence(const std::string& s1, const std::string s2);

int length_suffix_match(std::string& s1, std::string& s2);

#endif
