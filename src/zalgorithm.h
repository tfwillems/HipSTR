#include <string>
#include <vector>

namespace ZAlgorithm{
  /* 
   * For each position in s2, calculates the length of the matching prefix of s1 and s2[i...]
   * and stores it in num_matches[i]. The provided vector is cleared and sized appropriately.
   * Runs in O(length_of_s1 + length_of_s2)
   */
  void GetPrefixMatchCounts(const std::string& s1, const std::string& s2, std::vector<int>& num_matches);
  
  /* 
   * For each position in s2, calculates the length of the matching suffix of s1 and s2[0...i]
   * and stores it in num_matches[i]. The provided vector is cleared and sized appropriately.
   * Runs in O(length_of_s1 + length_of_s2)
   */
  void GetSuffixMatchCounts(const std::string& s1, const std::string& s2, std::vector<int>& num_matches);

  /* 
   * For each position i in s2 in the range [s2_start, s2_stop], calculates the length of 
   * the matching prefix of s1 and s2[i...] and stores it in num_matches[i-s2_start]. 
   * The provided vector is cleared and sized appropriately.
   * Runs in O(length_of_s1 + size_of_s2_range)
   */
  void GetPrefixMatchCounts(const std::string& s1, const std::string& s2, int s2_start, int s2_stop, std::vector<int>& num_matches);
  
  /* 
   * For each position i in s2 in the range [s2_start, s2_stop], calculates the length of 
   * the matching suffix of s1 and s2[0...i] and stores it in num_matches[i-s2_start]. 
   * The provided vector is cleared and sized appropriately.
   * Runs in O(length_of_s1 + size_of_s2_range)
   */
  void GetSuffixMatchCounts(const std::string& s1, const std::string& s2, int s2_start, int s2_stop, std::vector<int>& num_matches);
}
