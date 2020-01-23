#ifndef BAM_IO_H_
#define BAM_IO_H_

#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/cram/cram.h"
#include "htslib/sam.h"

#include "error.h"

// htslib encodes each base using a 4 bit integer
// This array converts each integer to its corresponding base
static const char HTSLIB_INT_TO_BASE[16] = {' ', 'A', 'C', ' ',
					    'G', ' ', ' ', ' ', 
					    'T', ' ', ' ', ' ', 
					    ' ', ' ', ' ', 'N'};

class CigarOp {
public:
  char Type;
  int32_t Length;
  
  CigarOp(char type, int32_t length){
    Type   = type;
    Length = length;
  }
};




class BamAlignment {
private:
  std::string bases_;
  std::string qualities_;
  std::vector<CigarOp> cigar_ops_;

  void ExtractSequenceFields();

  template<typename Iterator>
    int TrimHelper(Iterator begin, Iterator end, char base_qual_cutoff);

public:
  bam1_t *b_;
  std::string file_;
  std::string ref_, mate_ref_;
  bool built_;
  int32_t length_;
  int32_t pos_, end_pos_;

  BamAlignment(){
    b_       = bam_init1();
    built_   = false;
    length_  = -1;
    pos_     = 0;
    end_pos_ = -1;
  }

  BamAlignment(const BamAlignment &aln)
    : bases_(aln.bases_), qualities_(aln.qualities_), cigar_ops_(aln.cigar_ops_), file_(aln.file_), ref_(aln.ref_), mate_ref_(aln.mate_ref_){
    b_ = bam_init1();
    assert(bam_copy1(b_, aln.b_) != NULL);
    built_     = aln.built_;
    length_    = aln.length_;
    pos_       = aln.pos_;
    end_pos_   = aln.end_pos_;
  }

  BamAlignment& operator=(const BamAlignment& aln){
    assert(bam_copy1(b_, aln.b_) != NULL);
    file_      = aln.file_;
    ref_       = aln.ref_;
    mate_ref_  = aln.mate_ref_;
    built_     = aln.built_;
    length_    = aln.length_;
    pos_       = aln.pos_;
    end_pos_   = aln.end_pos_;
    bases_     = aln.bases_;
    qualities_ = aln.qualities_;
    cigar_ops_ = aln.cigar_ops_;
    return *this;
  }

  ~BamAlignment(){
    bam_destroy1(b_);
  }

  /* Number of bases */
  int32_t Length()              const { return length_;  }

  /* 0-based position where alignment starts*/
  int32_t Position()            const { return pos_;     }

  /* Non-inclusive end position of alignment */
  int32_t GetEndPosition()      const { return end_pos_; }

  /* Name of the read */
  std::string Name()            const { return std::string(bam_get_qname(b_)); }

  /* Name of the read's reference sequence */
  const std::string& Ref()      const { return ref_;              }

  /* Mapping quality score*/
  uint16_t MapQuality()         const { return b_->core.qual;     }

  const std::string& MateRef()  const { return mate_ref_;         }

  /* 0-based position where mate's alignment starts */
  int32_t MatePosition()        const { return b_->core.mpos;     }

  /* Name of file from which the alignment was read */
  const std::string& Filename() const { return file_;             }
  
  /* Sequenced bases */
  const std::string& QueryBases(){
    if (!built_) ExtractSequenceFields();
    return bases_;
  }

  /* Quality score for each base */
  const std::string& Qualities(){
    if (!built_) ExtractSequenceFields();
    return qualities_;
  }

  const std::vector<CigarOp>& CigarData(){
    if (!built_) ExtractSequenceFields();
    return cigar_ops_;
  }

  bool RemoveTag(const char tag[2]) const {
    uint8_t* tag_data = bam_aux_get(b_, tag);
    if (tag_data == NULL)
      return false;
    return (bam_aux_del(b_, tag_data) == 0);
  }

  bool HasTag(const char tag[2]) const { return bam_aux_get(b_, tag) != NULL; }

  /*
  bool AddCharTag(const char tag[2], char& value){
    if (HasTag(tag))
      return false;
    return (bam_aux_append(b_, tag, 'A', 1, (uint8_t*)&value) == 0);
  }

  bool AddIntTag(const char tag[2], int64_t& value){
    if (HasTag(tag))
      return false;
    return (bam_aux_append(b_, tag, 'i', ___, (uint8_t*)&value) == 0);
  }

  bool AddFloatTag(const char tag[2], double& value){
    if (HasTag(tag))
      return false;
    return (bam_aux_append(b_, tag, 'f', ___, (uint8_t*)&value) == 0);
  }
  */

  bool AddStringTag(const char tag[2], const std::string& value){
    if (HasTag(tag))
      return false;
    return (bam_aux_append(b_, tag, 'Z', value.size()+1, (uint8_t*)value.c_str()) == 0);
  }

  bool GetCharTag(const char tag[2], char& value) const {
    uint8_t* tag_data = bam_aux_get(b_, tag);
    if (tag_data == NULL)
      return false;
    value = bam_aux2A(tag_data);
    return true; // TO DO: Check errno
  }

  bool GetIntTag(const char tag[2], int64_t& value) const {
    uint8_t* tag_data = bam_aux_get(b_, tag);
    if (tag_data == NULL)
      return false;
    value = bam_aux2i(tag_data);
    return true; // TO DO: Check errno
  }

  bool GetFloatTag(const char tag[2], double& value) const {
    uint8_t* tag_data = bam_aux_get(b_, tag);
    if (tag_data == NULL)
      return false;
    value = bam_aux2f(tag_data);
    return true; // TO DO: Check errno
  }
  
  bool GetStringTag(const char tag[2], std::string& value) const{
    uint8_t* tag_data = bam_aux_get(b_, tag);
    if (tag_data == NULL)
      return false;
    char* ptr = bam_aux2Z(tag_data);
    value = std::string(ptr);
    return true; // TO DO: Check errno
  }  

  bool IsDuplicate()         const { return (b_->core.flag & BAM_FDUP)         != 0;}
  bool IsFailedQC()          const { return (b_->core.flag & BAM_FQCFAIL)      != 0;}
  bool IsMapped()            const { return (b_->core.flag & BAM_FUNMAP)       == 0;}
  bool IsMateMapped()        const { return (b_->core.flag & BAM_FMUNMAP)      == 0;}
  bool IsReverseStrand()     const { return (b_->core.flag & BAM_FREVERSE)     != 0;}
  bool IsMateReverseStrand() const { return (b_->core.flag & BAM_FMREVERSE)    != 0;}
  bool IsPaired()            const { return (b_->core.flag & BAM_FPAIRED)      != 0;}
  bool IsProperPair()        const { return (b_->core.flag & BAM_FPROPER_PAIR) != 0;}
  bool IsFirstMate()         const { return (b_->core.flag & BAM_FREAD1)       != 0;}
  bool IsSecondMate()        const { return (b_->core.flag & BAM_FREAD2)       != 0;}

  bool StartsWithSoftClip(){
    if (!built_) ExtractSequenceFields();
    if (cigar_ops_.empty())
      return false;
    return cigar_ops_.front().Type == 'S';
  }

  bool EndsWithSoftClip(){
    if (!built_) ExtractSequenceFields();
    if (cigar_ops_.empty())
      return false;
    return cigar_ops_.back().Type == 'S';
  }

  bool StartsWithHardClip(){
    if (!built_) ExtractSequenceFields();
    if (cigar_ops_.empty())
      return false;
    return cigar_ops_.front().Type == 'H';
  }

  bool EndsWithHardClip(){
    if (!built_) ExtractSequenceFields();
    if (cigar_ops_.empty())
      return false;
    return cigar_ops_.back().Type == 'H';
  }

  bool MatchesReference(){
    if (!built_) ExtractSequenceFields();
    for (auto cigar_iter = cigar_ops_.begin(); cigar_iter != cigar_ops_.end(); cigar_iter++)
      if (cigar_iter->Type != 'M' && cigar_iter->Type != '=')
	return false;
    return true;
  }

  void SetIsDuplicate(bool ok){
    if (ok) b_->core.flag |= BAM_FDUP;
    else    b_->core.flag &= (~BAM_FDUP);
  }

  void SetIsFailedQC(bool ok){
    if (ok) b_->core.flag |= BAM_FQCFAIL;
    else    b_->core.flag &= (~BAM_FQCFAIL);
  }

  void SetIsMapped(bool ok){
    if (ok) b_->core.flag &= (~BAM_FUNMAP);
    else    b_->core.flag |= BAM_FUNMAP;
  }

  void SetIsMateMapped(bool ok){
    if (ok) b_->core.flag &= (~BAM_FMUNMAP);
    else    b_->core.flag |= BAM_FMUNMAP;
  }

  void SetIsReverseStrand(bool ok){
    if (ok) b_->core.flag |= BAM_FREVERSE;
    else    b_->core.flag &= (~BAM_FREVERSE);
  }

  void SetIsMateReverseStrand(bool ok){
    if (ok) b_->core.flag |= BAM_FMREVERSE;
    else    b_->core.flag &= (~BAM_FMREVERSE);
  }

  void SetIsPaired(bool ok){
    if (ok) b_->core.flag |= BAM_FPAIRED;
    else    b_->core.flag &= (~BAM_FPAIRED);
  }

  void SetIsProperPair(bool ok){
    if (ok) b_->core.flag |= BAM_FPROPER_PAIR;
    else    b_->core.flag &= (~BAM_FPROPER_PAIR);
  }

  void SetIsFirstMate(bool ok){
    if (ok) b_->core.flag |= BAM_FREAD1;
    else    b_->core.flag &= (~BAM_FREAD1);
  }

  void SetIsSecondMate(bool ok){
    if (ok) b_->core.flag |= BAM_FREAD2;
    else    b_->core.flag &= (~BAM_FREAD2);
  }

  /*
   *  Trim an alignment that extends too far upstream or downstream of the provided region or has low base qualities on the ends
   *  Trims until either i) the base quality exceeds the provided threshold or ii) the alignment is fully within the provided region bounds
   *  Modifies the alignment such that subsequent calls to each function reflect the trimmming
   *  However, if the aligment is written to a new BAM file, the original alignment will be output
   */
  void TrimAlignment(int32_t min_read_start, int32_t max_read_stop);

  void TrimNumBases(int left_trim, int right_trim);

  void TrimLowQualityEnds(char min_base_qual);
};


std::string BuildCigarString(const std::vector<CigarOp>& cigar_data);


class ReadGroup {
 private:
  std::map<std::string, std::string> tag_dict_;

 public:
  ReadGroup(){}

  ReadGroup(const std::string& id, const std::string& sample, const std::string& library){
    SetID(id);
    SetSample(sample);
    SetLibrary(library);
  }

  bool HasTag(const std::string& tag) const {
    return tag_dict_.find(tag) != tag_dict_.end();
  }
  bool HasID()      const { return HasTag("ID"); }
  bool HasSample()  const { return HasTag("SM"); }
  bool HasLibrary() const { return HasTag("LB"); }

  const std::string& GetTag(const std::string& tag) const {
    auto iter = tag_dict_.find(tag);
    if (iter == tag_dict_.end())
      printErrorAndDie("Read group does not contain a " + tag + " tag");
    return iter->second;
  }
  const std::string& GetID()      const { return GetTag("ID"); }
  const std::string& GetSample()  const { return GetTag("SM"); }
  const std::string& GetLibrary() const { return GetTag("LB"); }

  void SetTag(const std::string& tag, const std::string& value){
    tag_dict_[tag] = value;
  }
  void SetID(const std::string& id)          { SetTag("ID", id);      }
  void SetSample(const std::string& sample)  { SetTag("SM", sample);  }
  void SetLibrary(const std::string& library){ SetTag("LB", library); }
};


class BamHeader {
 protected:
  std::map<std::string, int32_t> seq_indices_;
  std::vector<std::string> seq_names_;
  std::vector<uint32_t> seq_lengths_;
  std::vector<ReadGroup> read_groups_;

  void parse_read_groups(const char *text);

 public:
  bam_hdr_t *header_;

  explicit BamHeader(bam_hdr_t *header){
    header_ = header;
    for (int32_t i = 0; i < header_->n_targets; i++){
      seq_names_.push_back(std::string(header_->target_name[i]));
      seq_lengths_.push_back(header_->target_len[i]);
      seq_indices_.insert(std::pair<std::string, int32_t>(seq_names_.back(), i));
    }
    parse_read_groups(header_->text);
  }

  const std::vector<uint32_t>& seq_lengths()  const { return seq_lengths_; }
  const std::vector<std::string>& seq_names() const { return seq_names_;   }
  virtual const std::vector<ReadGroup>& read_groups(int file_index) const {
    assert(file_index == 0);
    return read_groups_;
  }

  int32_t num_seqs() const { return header_->n_targets; }
  int32_t ref_id(const std::string& ref) const {
    auto iter = seq_indices_.find(ref);
    if (iter == seq_indices_.end())
      return -1;
    return iter->second;
  }

  bool has_sequence(const std::string& chrom) const {
    return ref_id(chrom) != -1;
  }
  
  std::string ref_name(int32_t ref_id) const {
    if (ref_id == -1)
      return "*";
    if (ref_id >= 0 && ref_id < seq_names_.size())
      return seq_names_[ref_id];
    printErrorAndDie("Invalid reference ID provided to ref_name() function");
  }

  uint32_t ref_length(int32_t ref_id) const {
    if (ref_id >= 0 && ref_id < seq_lengths_.size())
      return seq_lengths_[ref_id];
    printErrorAndDie("Invalid reference ID provided to ref_length() function");
  }
};

void compare_bam_headers(const BamHeader* hdr_a, const BamHeader* hdr_b, const std::string& file_a, const std::string& file_b);


class BamMultiHeader : public BamHeader {
 private:
  std::string base_file_name_;
  std::vector< std::vector<ReadGroup> > read_groups_by_file_;

 public:
 BamMultiHeader(const BamHeader* header, const std::string& filename) : BamHeader(header->header_){
    base_file_name_ = filename;
    read_groups_by_file_.push_back(read_groups_);
    read_groups_.clear();
  }

  void add_header(const BamHeader* header, const std::string& file_name){
    // Ensure that the header sequences are consistent
    compare_bam_headers(this, header, base_file_name_, file_name);

    // Add the read groups
    parse_read_groups(header->header_->text);
    read_groups_by_file_.push_back(read_groups_);
    read_groups_.clear();
  }

  const std::vector<ReadGroup>& read_groups(int file_index) const { return read_groups_by_file_[file_index]; }
};




class BamCramReader {
private:
  samFile   *in_;
  bam_hdr_t *hdr_;
  hts_idx_t *idx_;
  std::string path_;
  BamHeader*  header_;
  bool shared_header_;
  bool cram_done_;

  // Instance variables for the most recently set region
  hts_itr_t *iter_;        // Iterator
  std::string chrom_;      // Chromosome
  int32_t     start_;      // Start position
  int32_t     end_;        // End position
  uint64_t    min_offset_; // Offset after first alignment. For BAMs, this is a memory offset, while for CRAMs its
                           // the index of the first alignment in the CRAM slice
  BamAlignment first_aln_; // First alignment
  bool reuse_first_aln_;

  // Private unimplemented copy constructor and assignment operator to prevent operations
  BamCramReader(const BamCramReader& other);
  BamCramReader& operator=(const BamCramReader& other);

  bool file_exists(const std::string& path){
    return (access(path.c_str(), F_OK) != -1);
  }

  void clear_cram_data_structures();

public:
  BamCramReader(const std::string& path, std::string fasta_path = "");

  const BamHeader* bam_header() const { return header_; }
  const std::string& path()     const { return path_;   }
  
  ~BamCramReader();

  bool GetNextAlignment(BamAlignment& aln);

  // Prepare the BAM/CRAM for reading the entire chromosome
  bool SetChromosome(const std::string& chrom);
  
  // Prepare the BAM/CRAM for reading all alignments overlapping the provided region
  bool SetRegion(const std::string& chrom, int32_t start, int32_t end);

  void use_shared_header(BamHeader* header){
    if (!shared_header_){
      bam_hdr_destroy(hdr_);
      delete header_;
    }

    shared_header_ = true;
    header_        = header;
    hdr_           = header_->header_;
  }
};






class BamCramMultiReader {
 private:
  std::vector<BamCramReader*> bam_readers_;
  std::vector<bool> reader_unset_;
  std::vector<BamAlignment> cached_alns_;
  std::vector<std::pair<int32_t, int32_t> > aln_heap_;
  int merge_type_;
  BamMultiHeader* multi_header_;

  // Instance variables for the most recently set region
  std::string chrom_;      // Chromosome
  int32_t     start_;      // Start position
  int32_t     end_;        // End position

  // Private unimplemented copy constructor and assignment operator to prevent operations
  BamCramMultiReader(const BamCramMultiReader& other);
  BamCramMultiReader& operator=(const BamCramMultiReader& other);

 public:
  const static int ORDER_ALNS_BY_POSITION = 0;
  const static int ORDER_ALNS_BY_FILE     = 1;

  BamCramMultiReader(const std::vector<std::string>& paths, std::string fasta_path = "", int merge_type = ORDER_ALNS_BY_POSITION, bool share_headers = true){
    if (paths.empty())
      printErrorAndDie("Must provide at least one file to BamCramMultiReader constructor");
    if (merge_type != ORDER_ALNS_BY_POSITION && merge_type != ORDER_ALNS_BY_FILE)
      printErrorAndDie("Invalid merge type provided to BamCramMultiReader constructor");
    for (size_t i = 0; i < paths.size(); i++){
      cached_alns_.push_back(BamAlignment());
      bam_readers_.push_back(new BamCramReader(paths[i], fasta_path));
      if (i == 0)
	multi_header_ = new BamMultiHeader(bam_readers_[i]->bam_header(), paths[i]);
      else {
	multi_header_->add_header(bam_readers_[i]->bam_header(), paths[i]);
	if (share_headers)
	  bam_readers_[i]->use_shared_header(multi_header_);
      }
    }
    merge_type_   = merge_type;
    reader_unset_ = std::vector<bool>(bam_readers_.size(), false);
    chrom_        = "";
    start_        = -1;
    end_          = -1;
  }

  ~BamCramMultiReader(){
    delete multi_header_;
    for (size_t i = 0; i < bam_readers_.size(); i++)
      delete bam_readers_[i];
  }

  int get_merge_type() const { return merge_type_; }
  const BamHeader* bam_header() const { return multi_header_; }

  bool SetRegion(const std::string& chrom, int32_t start, int32_t end);

  bool GetNextAlignment(BamAlignment& aln);
};








class BamWriter {
 private:
  BGZF* output_;

  // Private unimplemented copy constructor and assignment operator to prevent operations
  BamWriter(const BamWriter& other);
  BamWriter& operator=(const BamWriter& other);

 public:
  BamWriter(const std::string& path, const BamHeader* bam_header){
    std::string mode = "w";
    output_ = bgzf_open(path.c_str(), mode.c_str());
    if (output_ == NULL)
      printErrorAndDie("Failed to open BAM output file");
    if (bam_hdr_write(output_, bam_header->header_) == -1)
      printErrorAndDie("Failed to write the BAM header to the output file");
  }

  void Close(){
    if (bgzf_close(output_) != 0)
      printErrorAndDie("Failed to close BAM output file");
    output_ = NULL;
  }

  bool SaveAlignment(BamAlignment& aln){
    if (output_ == NULL)
      return false;
    return (bam_write1(output_, aln.b_) != -1);
  }

  ~BamWriter(){
    if (output_ != NULL)
      bgzf_close(output_);
  }
};

#endif
