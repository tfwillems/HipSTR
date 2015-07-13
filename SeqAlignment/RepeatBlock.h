#ifndef REPEAT_BLOCK_H_
#define REPEAT_BLOCK_H_

#include <algorithm>

#include "../error.h"
#include "../stutter_model.h"
#include "HapBlock.h"
#include "RepeatStutterInfo.h"

class RepeatBlock : public HapBlock {
 private:
    RepeatStutterInfo* repeat_info_;

 public:
 RepeatBlock(int32_t start, int32_t end, std::string ref_seq, int period, StutterModel* stutter_model): HapBlock(start, end, ref_seq){
      repeat_info_ = new RepeatStutterInfo(period, ref_seq, stutter_model); 
    }
    
    ~RepeatBlock(){ delete repeat_info_; }

    void add_alternate(std::string& alt){
      HapBlock::add_alternate(alt);
      repeat_info_->add_alternate_allele(alt);
    }

    RepeatStutterInfo* get_repeat_info(){ return repeat_info_; }


    HapBlock* reverse(){
      std::string rev_ref_seq = ref_seq_;
      std::reverse(rev_ref_seq.begin(), rev_ref_seq.end());
      RepeatBlock* rev_block = new RepeatBlock(start_, end_, rev_ref_seq, repeat_info_->get_period(), repeat_info_->get_stutter_model());
      for (unsigned int i = 0; i < alt_seqs_.size(); i++) {
	std::string alt = alt_seqs_[i];
	std::reverse(alt.begin(), alt.end());
	rev_block->add_alternate(alt);
      }
      rev_block->initialize();
      return rev_block;
    }
};


#endif
