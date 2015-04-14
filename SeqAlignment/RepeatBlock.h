#ifndef REPEAT_BLOCK_H_
#define REPEAT_BLOCK_H_

#include <algorithm>

#include "error.h"
#include "HapBlock.h"
#include "RepeatInfo.h"

class RepeatBlock : HapBlock {
 private:
    RepeatInfo* repeat_info_;

 public:
   RepeatBlock(int32_t start, int32_t end, std::string ref_seq, int period): HapBlock(start, end, ref_seq){
      repeat_info_ = new RepeatInfo(period, ref_seq); 
    }
    
    ~RepeatBlock(){
      delete repeat_info_;
    }

    void add_alternate(std::string alt){
      HapBlock::add_alternate(alt);
      repeat_info_->add_alternate_allele(alt);
    }

    RepeatInfo* get_repeat_info(){ return repeat_info_; }
};


#endif
