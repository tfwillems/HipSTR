#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "stringops.h"

void split_by_delim(const std::string &s, char delim, 
		    std::vector<std::string>& substrings){
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
    substrings.push_back(item);
}

std::string uppercase(std::string str){
  std::stringstream res;
  for (size_t i = 0; i < str.size(); i++)
    res << static_cast<char>(toupper(str[i]));
  return res.str();
}

bool string_starts_with(const std::string&s, std::string prefix){
  if (s.size() < prefix.size())
    return false;
  return s.substr(0, prefix.size()).compare(prefix) == 0;
}

bool string_ends_with(const std::string& s, std::string suffix){
  if (s.size() < suffix.size())
    return false;
  return s.substr(s.size()-suffix.size(), suffix.size()).compare(suffix) == 0;
}

bool orderByLengthAndSequence(const std::string& s1, const std::string& s2){
  if (s1.size() != s2.size())
    return s1.size() < s2.size();
  return s1.compare(s2) < 0;
}

int length_suffix_match(const std::string& s1, const std::string& s2){
  auto iter_1 = s1.rbegin();
  auto iter_2 = s2.rbegin();
  int num_matches = 0;
  for (; iter_1 != s1.rend() && iter_2 != s2.rend(); ++iter_1, ++iter_2){
    if (*iter_1 != *iter_2)
      break;
    num_matches++;
  }
  return num_matches;
}
