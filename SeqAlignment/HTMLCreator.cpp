#include <assert.h>
#include <iostream>
#include <sstream>

#include "../error.h"
#include "HTMLCreator.h"
#include "../stringops.h"

// TO DO: Specify a relative file path
std::string css_file = "/Users/tfwillems/Desktop/Coding/HipSTR/css/display_style.css";

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
  output << "<style type='text/css'>"    << "\n"

	 <<  ".ref {"                    << "\n"
	 << " color: white;"             << "\n"
	 << " font-family: Courier;"     << "\n"
	 << "}"                          << "\n"

	 << "td {"                       << "\n"
	 << " text-align:center;"        << "\n"
	 << " vertical-align:middle;"    << "\n"
	 << "}"                          << "\n"

	 << ".locustd {"                 << "\n"
	 << " font-style: italic;"       << "\n"
	 << " color: black;"             << "\n"
	 << "}"                          << "\n"

	 << ".snptd {"                   << "\n"
	 << "  background-color: gold;"  << "\n"
	 << "}"                          << "\n"
    
	 << ".inserttd {"                << "\n"
	 << " background-color: red;"    << "\n"
	 << "}"                          << "\n"
  
	 << ".spacertd {"                << "\n"
	 << " color: white;"             << "\n"
	 << "}"                          << "\n"
    
	 << ".reftable {"                << "\n"
	 << " color: white;"             << "\n"
         << " font-family: Courier;"     << "\n"
	 << " font-weight: bold;"        << "\n"
	 << " font-size: 13px;"          << "\n" 
	 << "}"                          << "\n"
    
	 << ".readtable {"               << "\n"
	 << " font-family: Courier;"     << "\n"
	 << " font-weight: normal;"      << "\n"
	 << " font-size: 13px;"          << "\n"
	 << "}"                          << "\n"
    
	 << "caption {"                  << "\n"
	 << " background: #dbb768;"      << "\n"
	 << " color:black;"              << "\n"
	 << " font-weight: bold;"        << "\n"
	 << " font-size: 1.1em;"         << "\n"
	 << " text-align: left;"         << "\n"
	 << "}"                          << "\n"
    
	 << ".hover {"                   << "\n"
	 << " background-color: pink;"   << "\n"
	 << "}"                          << "\n"
    
	 << ".Atd {"                     << "\n"
	 << " color: purple;"            << "\n"
	 << "}"                          << "\n"
    
	 << ".Ctd {"                     << "\n"
	 << " color: blue;"              << "\n"
	 << "}"                          << "\n"
    
	 << ".Gtd {"                     << "\n"
	 << " color: green;"             << "\n"
	 << "}"                          << "\n"
    
	 << ".Ttd {"                     << "\n"
	 << " color: orange;"            << "\n"
	 << "}"                          << "\n"
    
	 << ".vtd {"                     << "\n"
	 << " color: gray;"              << "\n"
	 << "}"                          << "\n"
    
	 << ".-td {"                     << "\n"
	 << " color: red;"               << "\n"
	 << "}"                          << "\n"
    
	 << ".Areftd {"                  << "\n"
	 << " background-color: purple;" << "\n"
	 << "}"                          << "\n"
    
	 << ".Creftd {"                  << "\n"
	 << " background-color: blue;"   << "\n"
	 << "}"                          << "\n"
    
	 << ".Greftd {"                  << "\n"
	 << " background-color: green;"  << "\n"
	 << "}"                          << "\n"
    
	 << ".Treftd {"                  << "\n"
	 << " background-color: orange;" << "\n"
	 << "}"                          << "\n"
    
	 << ".vreftd {"                  << "\n"
	 << " background-color: gray;"   << "\n"
	 << "}"                          << "\n"

	 << "</style>"                   << "\n";
}


void writeReferenceString(std::string& reference_string, 
			  std::ostream& output, 
			  std::string locus_id, 
			  std::vector<bool>& within_locus, 
			  bool draw_locus_id){
  output 
    << "<div class='reference'>\n"
    << "\t<table class=\"reftable\">" << "\n";

  if (draw_locus_id)
    output << "\t\t<caption>" << locus_id << "</caption>" << "\n";
  output << "\t\t<tr>" << "\n";
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
    << "\t\t</tr>" << "\n"
    << "\t</table>" << "\n"
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
      output << "\t\t<tr> <td style=\"text-align:left;\" colspan=\"" << label.size() << "\"> <font color=\"red\">" << label <<  "</font> </td> </tr>\n";
    }

    output << "<tr>";
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
    output << "\t\t</tr>" << "\n";
  }
}
