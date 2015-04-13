#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "error.h"
#include "pstream.h"
#include "RepeatRegion.h"

using namespace redi;

const int trf_score_thresh []     = {24, 22, 28, 28, 32, 34};
const std::string trf_cmd         = "trf - 2 7 7 80 10 20 6 -d -ngs -h";
const pstreams::pmode all3streams = pstreams::pstdin|pstreams::pstdout|pstreams::pstderr;

void get_repeat_regions(std::string& seq, int32_t seq_start, std::vector<RepeatRegion>& regions){
  regions.clear();
  
  //Create TRF command
  pstream ps(trf_cmd, all3streams);

  // Pipe input to TRF 
  ps << ">0" << std::endl
     << seq  << std::endl;
  ps.rdbuf()->peof(); // Signals end of stdin

  // Process results
  std::string res;
  int res_index = -1;
  while (std::getline(ps.out(), res)) {
    if (res.size() == 0)
      break;
    else if (res[0] == '@')
      res_index = atoi(res.substr(1).c_str());
    else {
      std::stringstream data; data << res;
      int rep_start, rep_end, score; std::string motif, sequence;
      int t_int; std::string t_str; float t_float;
      if (! (data 
	     >> rep_start >> rep_end 
	     >> t_int   >> t_float >> t_int >> t_int 
	     >> t_int   >> score   >> t_int >> t_int >> t_int >> t_int
	     >> t_float >> motif >> sequence))
	printErrorAndDie("Improperly formatted TRF resutls file");
      
      if (motif.size() > 1 &&  score >= trf_score_thresh[motif.size()-1])
	regions.push_back(RepeatRegion(seq_start+rep_start-1, seq_start+rep_end, motif, sequence));
    }
  }
  ps.close();
  sortRepeatRegions(regions);
}


bool compareRepeatRegions(const RepeatRegion& r1, const RepeatRegion& r2){
  if (r1.get_start() != r2.get_start())
    return r1.get_start() < r2.get_start();
  else
    return r1.get_stop()  > r2.get_stop();
}

void sortRepeatRegions(std::vector<RepeatRegion>& repeat_regions){
  std::sort(repeat_regions.begin(), repeat_regions.end(), compareRepeatRegions);
}
