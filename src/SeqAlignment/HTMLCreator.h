#ifndef HTML_CREATOR_H_
#define HTML_CREATOR_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

const char DELETION_CHAR = '-';
const char NOT_APP_CHAR  = '*';
const char SPACE_CHAR    = ' ';

void writeReferenceString(const std::string& reference_string,
			  std::ostream& output,
			  const std::string& locus_id,
			  const std::vector<bool>& within_locus,
			  bool draw_locus_id);

void writeAlignmentStrings(const std::string& reference_string,
			   std::ostream& output,
			   const std::string& locus_id,
			   const std::vector<std::string>& alignment_strings,
			   const std::vector<std::string>& alignment_samples,
			   const std::map<std::string, std::string>& sample_info,
			   bool highlight);
#endif
