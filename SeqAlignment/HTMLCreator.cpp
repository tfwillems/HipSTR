#include <assert.h>
#include <iostream>
#include <sstream>

#include "../error.h"
#include "HTMLCreator.h"
#include "../stringops.h"

std::map<char,std::string> initialize_colors(){
  std::map<char,std::string> m;
  m['A'] = "purple";
  m['a'] = "purple";
  m['C'] = "blue";
  m['c'] = "blue";
  m['G'] = "green";
  m['g'] = "green";
  m['T'] = "orange";
  m['t'] = "orange";
  m['N'] = "purple";
  m['n'] = "purple";
  m['-'] = "red";
  m['*'] = "gray";
  m[' '] = "white";
  m['x'] = "white";
  return m;
}

std::map<char, std::string> base_colors = initialize_colors();

void writeHeader(std::ostream& output){
  output << "#\t#\t#\t"
	 << "<style type='text/css'> "

	 <<  ".ref {"
	 << " color: white;"
	 << " font-family: Courier;"
	 << "} "

	 << "td {"
	 << " text-align:center;"
	 << " vertical-align:middle;"
	 << "} "

	 << ".locustd {"
	 << " font-style: italic;"
	 << " color: black;"
	 << "} "

	 << ".snptd {"
	 << "  background-color: gold;"
	 << "} "

	 << ".inserttd {"
	 << " background-color: red;"
	 << "} "

	 << ".spacertd {"
	 << " color: white;"
	 << "} "

	 << ".reftable {"
	 << " color: white;"
         << " font-family: Courier;"
	 << " font-weight: bold;"
	 << " font-size: 13px;"
	 << "} "

	 << ".readtable {"
	 << " font-family: Courier;"
	 << " font-weight: normal;"
	 << " font-size: 13px;"
	 << "} "

	 << "caption {"
	 << " background: #dbb768;"
	 << " color:black;"
	 << " font-weight: bold;"
	 << " font-size: 1.1em;"
	 << " text-align: left;"
	 << "} "

	 << ".hover {"
	 << " background-color: pink;"
	 << "} "

	 << ".Atd {"
	 << " color: purple;"
	 << "} "

	 << ".Ctd {"
	 << " color: blue;"
	 << "} "

	 << ".Gtd {"
	 << " color: green;"
	 << "} "

	 << ".Ttd {"
	 << " color: orange;"
	 << "} "

	 << ".vtd {"
	 << " color: gray;"
	 << "} "

	 << ".-td {"
	 << " color: red;"
	 << "} "

	 << ".Areftd {"
	 << " background-color: purple;"
	 << "} "

	 << ".Creftd {"
	 << " background-color: blue;"
	 << "} "

	 << ".Greftd {"
	 << " background-color: green;"
	 << "} "

	 << ".Treftd {"
	 << " background-color: orange;"
	 << "} "

	 << ".vreftd {"
	 << " background-color: gray;"
	 << "} "

	 << "</style>" << "\n";
}


void writeReferenceString(std::string& reference_string, 
			  std::ostream& output, 
			  std::string locus_id, 
			  std::vector<bool>& within_locus, 
			  bool draw_locus_id){
  output << locus_id << "\t"
	 << "<div class='reference'>"
	 << "\t<table class=\"reftable\">";

  if (draw_locus_id)
    output << " <caption>" << locus_id << "</caption> ";
  output << "<tr>";
  for (int i = 0; i < reference_string.size(); i++){
    std::string color = base_colors[reference_string[i]];
    if (within_locus[i]){
      output << "<td class=\"locustd\" style='background-color:" << color << "'>" << reference_string[i] << "</td>";
    }
    else {
      if (reference_string[i] == '*')
	output << "<td class=\"" << 'v' << "reftd\">" << reference_string[i] << "</td>";
      else
	output << "<td class=\"" << reference_string[i] << "reftd\">" << reference_string[i] << "</td>";
    }
  }
  output
    << "</tr>"
    << "</table>"
    << "</div>" << std::endl; 
}

void writeAlignmentStrings(std::string& reference_string, 
			   std::ostream& output, 
			   std::string locus_id,
			   std::vector<std::string>& alignment_strings, 
			   std::vector<std::string>& alignment_samples, 
			   std::map<std::string, std::string>& sample_info,
			   bool highlight){
  assert(alignment_strings.size() == alignment_samples.size());
  for (int i = 0; i < alignment_strings.size(); i++){
    if (i == 0 or alignment_samples[i-1].compare(alignment_samples[i]) != 0){
      std::stringstream vcf_ss;
      std::map<std::string, std::string>::iterator info_iter = sample_info.find(alignment_samples[i]);
      if (info_iter != sample_info.end())
	vcf_ss << info_iter->second;

      std::string label = alignment_samples[i] + ": " + vcf_ss.str();
      output << locus_id << "\t" << "<tr> <td style=\"text-align:left;\" colspan=\"" << label.size() << "\"> <font color=\"red\">" << label <<  "</font> </td> </tr>\n";
    }

    output << locus_id << "\t" << "<tr>";
    int j;
    for (j = 0; j < alignment_strings[i].size(); j++){
      if (alignment_strings[i][j] != SPACE_CHAR)
	break;
    }
    if (j > 0)
      output << "<td colspan=\"" << j << "\"> </td>";

    for (; j < alignment_strings[i].size(); j++){
      char c = alignment_strings[i][j]; 		
      std::string color = base_colors[c];
      bool snp    = (tolower(c) != tolower(reference_string[j]) && reference_string[j] != NOT_APP_CHAR && c != NOT_APP_CHAR && c != SPACE_CHAR && c != DELETION_CHAR);
      bool insert = (c != NOT_APP_CHAR && c != SPACE_CHAR && reference_string[j] == NOT_APP_CHAR);

      if (highlight && snp)
	output << "<td class=\"snptd\"> <font color=\"" << color << "\">" << c << "</font></td>";
      else if (highlight && insert)
	output << "<td class=\"inserttd\"> <font color=\"" << color << "\">" << c << "</font></td>";
      else if (c == SPACE_CHAR)
	output << "<td class=\"spacertd\">x</td>";
      else{
	if (c == '*')
	  output << "<td class=\"" << 'v' << "td\">" << '*' << "</td>";
	else
	  output << "<td class=\"" << c << "td\">" << c << "</td>";	  
      }
    }
    output << "</tr>\n";
  }
}
