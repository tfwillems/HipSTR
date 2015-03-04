#ifndef STRING_OPS_H_
#define STRING_OPS_H_

#include <string>
#include <vector>

void split_by_delim(const std::string &s, char delim, 
		    std::vector<std::string>& substrings);


#endif
