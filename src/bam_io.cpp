#include <sstream>

#include "bam_io.h"
#include "error.h"
#include "stringops.h"

void BamAlignment::ExtractSequenceFields(){
  int32_t length = b_->core.l_qseq;
  bases_         = std::string(length, ' ');
  qualities_     = std::string(length, ' '); 
  if (length == 0)
    return;
  
  // Rebuild the quality string
  // 33 is the reference point for the quality encoding
  uint8_t* quals = bam_get_qual(b_);
  for (int32_t i = 0; i < length; ++i)
    qualities_[i] = (char)(quals[i] + 33);
  
  // Rebuild the sequenced bases string
  uint8_t *bases = bam_get_seq(b_);
  for (int32_t i = 0; i < length; ++i)
    bases_[i] = HTSLIB_INT_TO_BASE[bam_seqi(bases, i)];
  
  // Rebuild the CIGAR operations
  int32_t num_cigar_ops = b_->core.n_cigar;
  uint32_t* cigars      = bam_get_cigar(b_);
  cigar_ops_.clear();
  for (int32_t i = 0; i < num_cigar_ops; ++i)
    cigar_ops_.push_back(CigarOp(bam_cigar_opchr(cigars[i]), bam_cigar_oplen(cigars[i])));
  
  built_ = true;
}


void BamHeader::parse_read_groups(){
  assert(read_groups_.empty());
  std::stringstream ss; ss << header_->text;
  std::string line;
  std::vector<std::string> tokens;
  while (std::getline(ss, line)){
    if (string_starts_with(line, "@RG")){
      split_by_delim(line, '\t', tokens);
      ReadGroup rg;
      for (int i = 1; i < tokens.size(); i++){
	if (string_starts_with(tokens[i], "ID:"))
	  rg.SetID(tokens[i].substr(3));
	else if (string_starts_with(tokens[i], "SM:"))
	  rg.SetSample(tokens[i].substr(3));
	else if (string_starts_with(tokens[i], "LB:"))
	  rg.SetLibrary(tokens[i].substr(3));
      }
      tokens.clear();
      read_groups_.push_back(rg);
    }
  }
}

bool BamCramReader::SetRegion(const std::string& chrom, int32_t start, int32_t end){
  bool reuse_offset = (min_offset_ != 0 && chrom.compare(chrom_) == 0 && start >= start_);
  if (reuse_offset && first_aln_.GetEndPosition() > start && first_aln_.Position() < end)
    reuse_offset = false;

  std::stringstream region;
  region << chrom << ":" << start+1 << "-" << end;
  std::string region_str = region.str();
  iter_ = sam_itr_querys(idx_, hdr_, region_str.c_str());

  if (iter_ == NULL && region_str.size() > 3 && region_str.substr(0, 3).compare("chr") == 0)
    iter_ = sam_itr_querys(idx_, hdr_, region_str.substr(3).c_str());

  if (iter_ != NULL){
    chrom_ = chrom;
    start_ = start;

    if (reuse_offset)
      if (iter_->n_off == 1 && min_offset_ >= iter_->off[0].u && min_offset_ <= iter_->off[0].v)
	iter_->off[0].u = min_offset_;

    min_offset_ = 0;
    return true;
  }
  else {
    chrom_      = "";
    start_      = -1;
    min_offset_ = 0;
    return false;
  }
}

bool BamCramReader::GetNextAlignment(BamAlignment& aln){
  if (iter_ == NULL) return false;

  if (sam_itr_next(in_, iter_, aln.b_) < 0){
    hts_itr_destroy(iter_);
    iter_ = NULL;
    return false;
  }

  if (min_offset_ == 0){
    first_aln_  = aln;
    min_offset_ = iter_->curr_off;
  }

  // Set up alignment instance variables
  aln.built_   = false;
  aln.file_    = path_;
  aln.length_  = aln.b_->core.l_qseq;
  aln.pos_     = aln.b_->core.pos;
  aln.end_pos_ = bam_endpos(aln.b_);
  return true;
}



bool BamCramMultiReader::SetRegion(const std::string& chrom, int32_t start, int32_t end){
  aln_heap_.clear();
  for (int32_t reader_index = 0; reader_index < bam_readers_.size(); reader_index++){
    if (!bam_readers_[reader_index]->SetRegion(chrom, start, end))
      return false;
    if (bam_readers_[reader_index]->GetNextAlignment(cached_alns_[reader_index])){
      if (merge_type_ == ORDER_ALNS_BY_POSITION)
	aln_heap_.push_back(std::pair<int32_t, int32_t>(-cached_alns_[reader_index].Position(), reader_index));
      else if (merge_type_ == ORDER_ALNS_BY_FILE)
	aln_heap_.push_back(std::pair<int32_t, int32_t>(-reader_index, reader_index));
      else
	printErrorAndDie("Invalid merge order in SetRegion()");
    }
  }
  std::make_heap(aln_heap_.begin(), aln_heap_.end());
  return true;
}

bool BamCramMultiReader::GetNextAlignment(BamAlignment& aln){
  if (aln_heap_.empty())
    return false;
  std::pop_heap(aln_heap_.begin(), aln_heap_.end());
  int32_t reader_index = aln_heap_.back().second;
  aln_heap_.pop_back();

  // Assign optimal alignment to provided reference
  aln = cached_alns_[reader_index];

  // Add reader's next alignment to the cache
  if (bam_readers_[reader_index]->GetNextAlignment(cached_alns_[reader_index])){
    if (merge_type_ == ORDER_ALNS_BY_POSITION)
      aln_heap_.push_back(std::pair<int32_t, int32_t>(-cached_alns_[reader_index].Position(), reader_index));
    else if (merge_type_ == ORDER_ALNS_BY_FILE)
      aln_heap_.push_back(std::pair<int32_t, int32_t>(-reader_index, reader_index));
    else
      printErrorAndDie("Invalid merge order in GetNextAlignment()");

    std::push_heap(aln_heap_.begin(), aln_heap_.end());
  }
  return true;
}




void compare_bam_headers(const BamHeader* hdr_a, const BamHeader* hdr_b, const std::string& file_a, const std::string& file_b){
  std::stringstream error_msg;
  if (hdr_a->num_seqs() != hdr_b->num_seqs()){
    error_msg << "BAM header mismatch issue. BAM headers for files " << file_a << " and " << file_b << " must have the same number of reference sequences";
    printErrorAndDie(error_msg.str());
  }

  for (int32_t i = 0; i < hdr_a->num_seqs(); ++i){
    if (hdr_a->ref_name(i).compare(hdr_b->ref_name(i)) != 0){
      error_msg << "BAM header mismatch issue. Order of reference sequences in BAM headers for files " << file_a << " and " << file_b << " must match";
      printErrorAndDie(error_msg.str());
    }
    if (hdr_a->ref_length(i) != hdr_b->ref_length(i)){
      error_msg << "BAM header mismatch issue. Length of reference sequences in BAM headers for files " << file_a << " and " << file_b << " must match";
      printErrorAndDie(error_msg.str());
    }
  }
}


void BamAlignment::TrimAlignment(int32_t min_read_start, int32_t max_read_stop, char min_base_qual){
  if (!built_)
    ExtractSequenceFields();
  assert(bases_.size() == qualities_.size());

  int ltrim = 0;
  int32_t start_pos = pos_;
  while ((start_pos < min_read_start) && cigar_ops_.size() > 0){
    // Check if we should stop trimming b/c the quality score is above the threshold
    bool qual_above_thresh = false;
    switch (cigar_ops_.front().Type){
    case 'M': case '=': case 'X': case 'I': case 'S':
      qual_above_thresh = (qualities_[ltrim] > min_base_qual);
      break;
    case 'D': case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in trimAlignment");
      break;
    }
    if (qual_above_thresh)
      break;

    switch(cigar_ops_.front().Type){
    case 'M': case '=': case 'X':
      ltrim++;
      start_pos++;
      break;
    case 'D':
      start_pos++;
      break;
    case 'I': case 'S':
      ltrim++;
      break;
    case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in TrimAlignment");
      break;
    }
    if (cigar_ops_.front().Length == 1)
      cigar_ops_.erase(cigar_ops_.begin(), cigar_ops_.begin()+1);
    else
      cigar_ops_.front().Length--;
  }

  int rtrim = 0, qual_string_len = qualities_.size()-1;
  int32_t end_pos = end_pos_;
  while ((end_pos > max_read_stop) && cigar_ops_.size() > 0){
    // Check if we should stop trimming b/c the quality score is above the threshold
    bool qual_above_thresh = false;
    switch(cigar_ops_.back().Type){
    case 'M': case '=': case 'X': case 'I': case 'S':
      qual_above_thresh = (qualities_[qual_string_len-rtrim] > min_base_qual);
      break;
    case 'D': case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in TrimAlignment");
      break;
    }
    if (qual_above_thresh)
      break;

    switch(cigar_ops_.back().Type){
    case 'M': case '=': case 'X':
      rtrim++;
      end_pos--;
      break;
    case 'D':
      end_pos--;
      break;
    case 'I': case 'S':
      rtrim++;
      break;
    case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in trimAlignment");
      break;
    }
    if (cigar_ops_.back().Length == 1)
      cigar_ops_.pop_back();
    else
      cigar_ops_.back().Length--;
  }

  assert(ltrim+rtrim <= bases_.size());
  bases_     = bases_.substr(ltrim, bases_.size()-ltrim-rtrim);
  qualities_ = qualities_.substr(ltrim, qualities_.size()-ltrim-rtrim);
  length_   -= (ltrim + rtrim);
  pos_       = start_pos;
  end_pos_   = end_pos;
}

void BamAlignment::TrimLowQualityEnds(char min_base_qual){
  return TrimAlignment(end_pos_+1, pos_-1, min_base_qual);
}
