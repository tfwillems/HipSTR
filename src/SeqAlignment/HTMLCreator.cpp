#include <assert.h>
#include <iostream>
#include <sstream>

#include "../error.h"
#include "HTMLCreator.h"
#include "../stringops.h"

void writeReferenceString(const std::string& reference_string,
			  std::ostream& output,
			  const std::string& locus_id,
			  const std::vector<bool>& within_locus,
			  bool draw_locus_id){
  output << locus_id << "\t" << "ALL" << "\t" << "<div>" << "\t<table class=\"reftable\">";
  if (draw_locus_id)
    output << " <caption>" << locus_id << "</caption> ";
  output << "\n";

  output << locus_id << "\t" << "ALL" << "\t" << "<tr style='font-weight: bold' class=\"reference\">" << "0 ";
  for (int i = 0; i < reference_string.size(); i++)
    output << reference_string[i];
  output << "</tr>" << std::endl;
}

void writeAlignmentStrings(const std::string& reference_string,
			   std::ostream& output,
			   const std::string& locus_id,
			   const std::vector<std::string>& alignment_strings,
			   const std::vector<std::string>& alignment_samples,
			   const std::map<std::string, std::string>& sample_info,
			   bool highlight){
  assert(alignment_strings.size() == alignment_samples.size());
  for (int i = 0; i < alignment_strings.size(); i++){
    if (i == 0 or alignment_samples[i-1].compare(alignment_samples[i]) != 0){
      std::stringstream vcf_ss;
      auto info_iter = sample_info.find(alignment_samples[i]);
      if (info_iter != sample_info.end())
	vcf_ss << info_iter->second;

      std::string label = alignment_samples[i] + ": " + vcf_ss.str();
      output << locus_id << "\t" << alignment_samples[i] << "\t"
	     << "<tr> <td class=\"samplename\" style=\"text-align:left;\" colspan=\"" << label.size() 
	     << "\"> <font color=\"red\">" << label <<  "</font> </td> </tr>\n";
    }

    output << locus_id << "\t" << alignment_samples[i] << "\t" << "<tr>";
    int j;
    for (j = 0; j < alignment_strings[i].size(); j++)
      if (alignment_strings[i][j] != SPACE_CHAR)
	break;
    output << j << " ";

    for (; j < alignment_strings[i].size(); j++){
      char c = alignment_strings[i][j]; 		
      bool snp    = (tolower(c) != tolower(reference_string[j]) && reference_string[j] != NOT_APP_CHAR && c != NOT_APP_CHAR && c != SPACE_CHAR && c != DELETION_CHAR);
      bool insert = (c != NOT_APP_CHAR && c != SPACE_CHAR && reference_string[j] == NOT_APP_CHAR);

      if (highlight && snp){
	switch(c){
	case 'A':
	  output << "H";
	  break;
	case 'C':
	  output << "I";
	  break;
	case 'G':
	  output << "J";
	  break;
	case 'T':
	  output << "K";
	  break;
	case 'N':
	  output << "L";
	  break;
	default:
	  printErrorAndDie("Invalid base for HTML creation: " + std::string(1, c));
	  break;
	}
      }
      else if (highlight && insert){
	switch(c){
	case 'A':
	  output << "a";
	  break;
	case 'C':
	  output << "c";
	  break;
	case 'G':
	  output << "g";
	  break;
	case 'T':
	  output << "t";
	  break;
	case 'N':
	  output << "n";
	  break;
	default:
	  printErrorAndDie("Invalid base for HTML creation: " + std::string(1, c));
	  break;
	}
      }
      else if (c == SPACE_CHAR)
	output << "x";
      else {
	if (c == '*')
	  output << "*";
	else
	  output << c;
      }
    }
    output << "</tr>\n";
  }
}
