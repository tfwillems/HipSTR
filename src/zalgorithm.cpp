#include <iostream>
#include <sstream>

#include <stdlib.h>

#include "error.h"
#include "zalgorithm.h"

namespace ZAlgorithm{
  static void suffix_helper(const std::string& s1, const std::string& s2, int s2_left, int s2_right,
			    std::vector<int>& s1_matches, std::vector<int>& num_matches){
    num_matches  = std::vector<int>(s2_right - s2_left + 1, -1);
    int leftmost = s2_right+1, right_index = s2_right+1;
    for (int i = s2_right; i >= s2_left; i--){
      if (i <= leftmost){
	int index_a = s1.size()-1, index_b = i;
	while (index_a >= 0 && index_b >= 0 && (char)tolower(s1[index_a]) == (char)tolower(s2[index_b])){
	  index_a--;
	  index_b--;
	}
	num_matches[i-s2_left] = i - index_b;
	if (index_b < i){
	  right_index = i;
	  leftmost    = index_b + 1;
	}
      }
      else {
	int twin     = i - right_index + s1.size()-1;
	int new_left = i - s1_matches[twin] + 1;
	if (new_left > leftmost)
	  num_matches[i-s2_left] = s1_matches[twin];
	else if (new_left < leftmost)
	  num_matches[i-s2_left] = i-leftmost+1;
	else {
	  int index_a = s1.size()-2-i+leftmost, index_b = leftmost-1;
	  while (index_a >= 0 && index_b >= 0 && (char)tolower(s1[index_a]) == (char)tolower(s2[index_b])){
	    index_a--;
	    index_b--;
	  }
	  num_matches[i-s2_left] = i-index_b;
	  right_index            = i;
	  leftmost               = index_b + 1;
	}
      }
    }
  }

  static void prefix_helper(const std::string& s1, const std::string& s2, int s2_left, int s2_right,
			    std::vector<int>& s1_matches, std::vector<int>& num_matches, int offset){
    num_matches = std::vector<int>(s2_right-s2_left+1+offset, -1);
    int rightmost = -1, left_index = -1;
    const int s1_size = static_cast<int>(s1.size());
    const int s2_size = static_cast<int>(s2.size());
    for (int i = s2_left; i <= s2_right; i++){
      if (i >= rightmost){
	int index_a = 0, index_b = i;
	while (index_a < s1_size && index_b < s2_size && (char)tolower(s1[index_a]) == (char)tolower(s2[index_b])){
	  index_a++;
	  index_b++;
	}
	num_matches[i-s2_left+offset] = index_b - i;
	if (index_b > i){
	  left_index = i;
	  rightmost  = index_b - 1;
	}
      }
      else {
	int twin      = i - left_index;
	int new_right = i + s1_matches[twin] - 1;
	if (new_right < rightmost)
	  num_matches[i-s2_left+offset] = s1_matches[twin];
	else if (new_right > rightmost)
	  num_matches[i-s2_left+offset] = rightmost-i+1;
	else {
	  int index_a = rightmost+1-i, index_b = rightmost+1;
	  while (index_a < s1_size && index_b < s2_size && (char)tolower(s1[index_a]) == (char)tolower(s2[index_b])){
	    index_a++;
	    index_b++;
	  }
	  num_matches[i-s2_left+offset] = index_b - i;
	  left_index     = i;
	  rightmost      = index_b - 1;
	}
      }
    }
  }

  void GetPrefixMatchCounts(const std::string& s1, const std::string& s2, std::vector<int>& num_matches) {
    std::vector<int> s1_matches;
    prefix_helper(s1, s1, 1, s1.size()-1, s1_matches, s1_matches,  1);
    prefix_helper(s1, s2, 0, s2.size()-1, s1_matches, num_matches, 0);
  }  

  void GetSuffixMatchCounts(const std::string& s1, const std::string& s2, std::vector<int>& num_matches) {
    std::vector<int> s1_matches;
    suffix_helper(s1, s1, 0, s1.size()-2, s1_matches, s1_matches);
    suffix_helper(s1, s2, 0, s2.size()-1, s1_matches, num_matches);
  }  

  void GetPrefixMatchCounts(const std::string& s1, const std::string& s2, int s2_start, int s2_stop, std::vector<int>& num_matches) {
    if (s2_start < 0 or s2_stop >= static_cast<int>(s2.size()))
	printErrorAndDie("Invalid string indices provided to GetPrefixMatchCounts");
    std::vector<int> s1_matches;
    prefix_helper(s1, s1, 1, s1.size()-1, s1_matches, s1_matches,  1);
    prefix_helper(s1, s2, s2_start, s2_stop, s1_matches, num_matches, 0);
  }  

  void GetSuffixMatchCounts(const std::string& s1, const std::string& s2, int s2_start, int s2_stop, std::vector<int>& num_matches) {
    if (s2_start < 0 or s2_stop >= static_cast<int>(s2.size()))
	printErrorAndDie("Invalid string indices provided to GetSuffixMatchCounts");
    std::vector<int> s1_matches;
    suffix_helper(s1, s1, 0, s1.size()-2, s1_matches, s1_matches);
    suffix_helper(s1, s2, s2_start, s2_stop, s1_matches, num_matches);
  }
}

