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
  for (size_t i = 0; i < bam_readers_.size(); i++){
    if (!bam_readers_[i]->SetRegion(chrom, start, end))
      return false;
    if (bam_readers_[i]->GetNextAlignment(cached_alns_[i]))
      aln_heap_.push_back(std::pair<int32_t,size_t>(-cached_alns_[i].Position(), i));
  }
  std::make_heap(aln_heap_.begin(), aln_heap_.end());
  return true;
}

bool BamCramMultiReader::GetNextAlignment(BamAlignment& aln){
  if (aln_heap_.empty())
    return false;
  std::pop_heap(aln_heap_.begin(), aln_heap_.end());
  size_t reader_index = aln_heap_.back().second;
  aln_heap_.pop_back();

  //Assign optimal alignment to provided reference
  aln = cached_alns_[reader_index];

  // Add reader's next alignment to the cache
  if (bam_readers_[reader_index]->GetNextAlignment(cached_alns_[reader_index])){
    aln_heap_.push_back(std::pair<int32_t, size_t>(-cached_alns_[reader_index].Position(), reader_index));
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
