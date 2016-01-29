#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

#include "../error.h"		    
#include "NWNoRefEndPenalty.h"

class IndelTracker {
private:
  const static unsigned int MAX_BITS = 64;
  const static unsigned int BITS_PER_INDEL = 9;
  uint64_t val_;
  unsigned int al_length_;
  unsigned int num_indels_;

public:
  IndelTracker(uint64_t val, unsigned int num_indels, unsigned int al_length){
    val_        = val;
    num_indels_ = num_indels;
    al_length_  = al_length;
  }

  IndelTracker(){
    val_        = 0;
    num_indels_ = 0;
    al_length_  = 0;
  }

  IndelTracker copy(){
    return IndelTracker(val_, num_indels_, al_length_+1);
  }

  IndelTracker add(bool die_on_overflow){
    unsigned int loc = al_length_+1;
    if ((num_indels_+1)*BITS_PER_INDEL > MAX_BITS){
      if (die_on_overflow){
	uint64_t tmp  = val_;
	uint64_t mask = (1 << BITS_PER_INDEL)-1;
	for (unsigned int i = 0; i < num_indels_; i++){
	  std::cerr << (tmp & mask) << " ";
	  tmp >>= BITS_PER_INDEL; 
	}
	
	std::cerr << num_indels_+1 << std::endl;
	printErrorAndDie("Maximum number of indels exceeded in IndelTracker");
      }
      else {
	set_max_val();
	val_ -= 1;
      }
    }
    if (static_cast<unsigned int>((1 << BITS_PER_INDEL)-1) < loc)
      printErrorAndDie("Location is too large for IndelTracker: " + std::to_string(loc));
    return IndelTracker((val_ << BITS_PER_INDEL) | loc, num_indels_+1, al_length_+1);
  }

  unsigned int num_indels(){
    return num_indels_;
  }

  void set_max_val(){
    val_        = -1;
    num_indels_ = MAX_BITS/BITS_PER_INDEL;
    al_length_  = 0;
  }

  bool less_than(IndelTracker& other){
    return val_ < other.val_;
  }

  bool greater_than(IndelTracker& other){
    return val_ > other.val_;
  }
  
  static unsigned int max_indels(){
    return MAX_BITS/BITS_PER_INDEL;
  }
};



namespace NWNoRefEndPenalty {
  // Create alignment scoring matrix
  const float a           =  2.0; // Match
  const float b           = -2.0; // Mismatch
  const float GAPOPEN     =  5.0;
  const float GAPEXTEND   =  0.125;
  
  // Match and mismatch as expected for A, C, G and T
  // N character matches all characters
  const float s[5][5]     = {{ a, b, b, b, a },
			     { b, a, b, b, a },
			     { b, b, a, b, a },
			     { b, b, b, a, a },
			     { a, a, a, a, a }};

  // Large value used to penalize impossible configurations
  const float LARGE = 1000000;

  // Convert character to integer representing the base's
  // index in the scoring matrix
  int base_to_int(char c){
    c = toupper(c);
    switch(c){
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    case 'N':
      return 4;
    default:
      printErrorAndDie("Invalid character '" + std::string(1, c) + "' in read");
    }
    return -1;
  }

  float bestIndex(float s1, float s2, float s3, int* ptr){
    if (s2 > s1){
      if (s2 > s3){
	*ptr = 1;
	return s2;
      }
      else {
	*ptr = 2;
	return s3;
      }
    }
    else {
      if (s3 > s1){
	*ptr = 2;
	return s3;
      }
      else {
	*ptr = 0;
	return s1;
      }
    }
  }

  void findOptimalStop(int& L1, int& L2, 
		       std::vector<float>& M, 
		       std::vector<float>& Iref, 
		       std::vector<float>& Iread, 
		       float& best_val, int& best_col, int& best_type){
    best_val   = -LARGE;
    best_col   = -1;
    best_type  = -1;
    int column = 0;
    for (int i = L2*(L1+1); i < (L2+1)*(L1+1); i++){
      if (M[i] >= best_val){
	best_val  = M[i];
	best_col  = column;
	best_type = 0;
      }
    
      if (Iref[i] > best_val){
	best_val  = Iref[i];
	best_col  = column;
	best_type = 1;
      }
    
      if (Iread[i] > best_val){
	best_val  = Iread[i];
	best_col  = column;
	best_type = 2;
      }
      column++;
    }
  }

  void nw_helper(std::vector<float>& M,    std::vector<float>& Iref,    std::vector<float>& Iread, 
		 std::vector<int>& traceM, std::vector<int>& traceIref, std::vector<int>& traceIread,
		 const std::string& refseq, const std::string& readseq){
    int L1 = refseq.length();
    int L2 = readseq.length();

    // Various variables used in the matrix calculations
    int ref_base, read_base, oindex, nindex;
    float s1, s2, s3;
    int c;

    std::vector<int> ref_base_ints, read_base_ints;
    ref_base_ints.reserve(L1); read_base_ints.reserve(L2);
    for (unsigned int i = 0; i < refseq.size(); i++)
      ref_base_ints[i] = base_to_int(refseq[i]);
    for (unsigned int i = 0; i < readseq.size(); i++)
      read_base_ints[i] = base_to_int(readseq[i]);

    // Fill in the 3 matrices using dynamic programming
    for (int i = 1; i <= L2; i++){
      read_base = read_base_ints[i-1];

      for (int j = 1; j <= L1; j++){
	nindex    = i*(L1+1)+j;
	ref_base  = ref_base_ints[j-1];

	// Update M matrix (examine (i-1, j-1))
	oindex          = (i-1)*(L1+1)+(j-1);
	s1              = M[oindex];
	s2              = Iref[oindex];
	s3              = Iread[oindex];
	M[nindex]       = bestIndex(s1, s2, s3, &c) + s[ref_base][read_base];
	traceM[nindex]  = c;

	// Update Iref matrix (examine (i,j-1))
	oindex             = i*(L1+1) + (j-1);
	s1                 = M[oindex]     - GAPOPEN;
	s2                 = Iref[oindex]  - GAPEXTEND;
	s3                 = Iread[oindex] - GAPOPEN;
	Iref[nindex]       = bestIndex(s1, s2, s3, &c);
	traceIref[nindex]  = c;

	// Update Iread matrix (examine (i-1,j))
	oindex              = (i-1)*(L1+1) + j;
	s1                  = M[oindex]     - GAPOPEN;
	s2                  = Iref[oindex]  - GAPOPEN;
	s3                  = Iread[oindex] - GAPEXTEND;
	Iread[nindex]       = bestIndex(s1, s2, s3, &c);
	traceIread[nindex]  = c;
      }
    }
  }


  void traceAlignment(int best_col, 
		      int best_type, 
		      int L1, int L2,
		      std::vector<int>& traceM, 
		      std::vector<int>& traceIref, 
		      std::vector<int>& traceIread,
		      const std::string& refseq, 
		      const std::string& readseq,
		      std::string& ref_seq_al, 
		      std::string& read_seq_al,
		      std::vector<BamTools::CigarOp>& cigar_list){
    cigar_list.clear();
    std::stringstream refseq_ss, readseq_ss, cigar_ss;
  
    // Handle trailing gaps
    for(int i = L1; i > best_col; i--){
      refseq_ss  << refseq.at(i-1);
      readseq_ss << "-";
    }

    // Traceback the optimal alignment
    int best_row = L2;
    std::string raw_cigar;
    int index;
    while (best_row > 0){
      index = best_row*(L1+1) + best_col;
      if (best_type == 0){
	// M
	refseq_ss  << refseq.at(best_col-1);
	readseq_ss << readseq.at(best_row-1);

	if (base_to_int(refseq.at(best_col-1)) == base_to_int(readseq.at(best_row-1)))
	  cigar_ss << "=";
	else
	  cigar_ss << "X";

	best_type   = traceM[index];
	best_row--;
	best_col--;
      } 
      else if (best_type == 1){
	//Iref
	refseq_ss  << refseq.at(best_col-1);
	readseq_ss << "-";
	cigar_ss   << "D";
	best_type   = traceIref[index];
	best_col--;
      } 
      else if (best_type == 2){
	// Iread
	refseq_ss  << "-";
	readseq_ss << readseq.at(best_row-1);
	cigar_ss   << "I";
	best_type   = traceIread[index];
	best_row--;
      } 
      else
	printErrorAndDie("Invalid matrix type in Needleman-Wunsch alignment");
    }

    // Handle leading gaps
    for (int i = best_col; i > 0; i--){
      refseq_ss  << refseq.at(i-1);
      readseq_ss << "-";
    }
  
    // Order alignment front to back
    ref_seq_al  = refseq_ss.str();
    read_seq_al = readseq_ss.str();
    raw_cigar   = cigar_ss.str();
    reverse(ref_seq_al.begin(),  ref_seq_al.end());
    reverse(read_seq_al.begin(), read_seq_al.end());
    reverse(raw_cigar.begin(),   raw_cigar.end());

    // Simplify cigar string
    char cigar_char = raw_cigar[0];
    int  num        = 1;
    char new_cigar_char;
    for(unsigned int i = 1; i < raw_cigar.length(); i++){
      new_cigar_char = raw_cigar[i];
      if (new_cigar_char != cigar_char){
	cigar_list.push_back(BamTools::CigarOp(cigar_char, num));
	num = 1;
	cigar_char = new_cigar_char;
      }
      else
	num += 1;
    }
    cigar_list.push_back(BamTools::CigarOp(cigar_char, num));
    /*
    if (cigar_list.front().Type == 'I')
      cigar_list.front().Type = 'S';
    if (cigar_list.back().Type == 'I')
      cigar_list.back().Type = 'S';
    */
  }

 

  void initMatrices(std::vector<float>& M,    std::vector<float>& Iref,    std::vector<float>& Iread,
		    std::vector<int>& traceM, std::vector<int>& traceIref, std::vector<int>& traceIread,
		    int L1, int L2){
    M[0]     = 0.0;
    Iref[0]  = -LARGE; // Impossible
    Iread[0] = -LARGE; // Impossible

    // Fill in row (0,n)
    for(int i = 1; i < L1+1; i++){
      // No penalty for leading affine gap in reference seq 
      Iref[i]       = 0.0;  
      traceIref[i]  = 1;

      // Impossible
      Iread[i]      = -LARGE;
      traceIread[i] = -1;

      // Impossible
      M[i]      = -LARGE;
      traceM[i] = -1;
    }

    // Fill in column (n, 0)
    for(int i = 1; i < L2+1; i++){
       int index = i*(L1+1);

      // Penalty for leading affine gap in read sequence
      Iread[index]      = -GAPOPEN-(i-1)*GAPEXTEND;
      traceIread[index] = 2;

      // Impossible
      Iref[index]       = -LARGE; 
      traceIref[index]  = -1;

      // Impossible
      M[index]      = -LARGE;      
      traceM[index] = -1;
    }
  }


  void Align(const std::string& ref_seq, const std::string& read_seq,
	     std::string& ref_seq_al, std::string& read_seq_al,
	     float* score, std::vector<BamTools::CigarOp>& cigar_list) {
    int L1       = ref_seq.length();
    int L2       = read_seq.length();
    int mat_size = (L1+1)*(L2+1);

    // Scoring matrices:
    //  M:     Ref and read bases aligned
    //  Iref:  Ref base aligned with gap
    //  Iread: Read base aligned with gap
    std::vector<float> M(mat_size), Iref(mat_size), Iread(mat_size);

    // Traceback matrices
    std::vector<int> traceM(mat_size), traceIref(mat_size), traceIread(mat_size);

    // Initialize matrices
    initMatrices(M, Iref, Iread, traceM, traceIref, traceIread, L1, L2);

    // Fill out scoring and traceback matrices using variant of NW algorithm
    nw_helper(M, Iref, Iread, traceM, traceIref, traceIread, ref_seq, read_seq);

    // Find the best ending point for the alignment
    float best_val;
    int best_col, best_type;
    findOptimalStop(L1, L2, M, Iref, Iread, best_val, best_col, best_type);
    *score = best_val;

    // Construct the alignment strings and CIGAR string using the traceback 
    // matrices and the optimal end position
    traceAlignment(best_col, best_type, L1, L2, traceM, traceIref, traceIread,
		   ref_seq, read_seq, ref_seq_al, read_seq_al, cigar_list);
  }

  
  float bestIndex(float s1, float s2, float s3, IndelTracker* t1, IndelTracker* t2, IndelTracker* t3, int* best_type){
    IndelTracker max_val; 
    max_val.set_max_val();
    float best_val           = std::max(s1, std::max(s2, s3));
    IndelTracker* best_track = &max_val;

    if (s1 == best_val){
      if (t1->less_than(*best_track)){
	*best_type = 0;
	best_track = t1;
      }
    }

    if (s2 == best_val){
      if (t2->less_than(*best_track)){
	*best_type = 1;
	best_track = t2;
      }
    }

    if (s3 == best_val){
      if (t3->less_than(*best_track)){
	*best_type = 2;
	best_track = t3;
      }
    }

    if (!best_track->less_than(max_val))
      printErrorAndDie("Logical error in bestIndex()");
    return best_val;
  }



   void left_align_helper(std::vector<float>& M,     std::vector<float>& Iref,    std::vector<float>& Iread, 
			  std::vector<int>& traceM,  std::vector<int>& traceIref, std::vector<int>& traceIread,
			  const std::string& refseq, const std::string& readseq,
			  int start_col, int end_col, int max_indels){
    int L1 = refseq.length();
    int L2 = readseq.length();

    // Indel records
    int ntracks = end_col - start_col + 2;
    std::vector<IndelTracker> prev_M_tracks(ntracks), prev_Iref_tracks(ntracks), prev_Iread_tracks(ntracks);
    std::vector<IndelTracker> M_tracks(ntracks), Iref_tracks(ntracks), Iread_tracks(ntracks);

    // Various variables used in the matrix calculations
    int ref_base, read_base, oindex, nindex;
    int otindex, ntindex;
    float s1, s2, s3;
    int c;
    IndelTracker t1, t2, t3;

    std::vector<int> ref_base_ints, read_base_ints;
    ref_base_ints.reserve(L1); read_base_ints.reserve(L2);
    for (unsigned int i = 0; i < refseq.size(); i++)
      ref_base_ints[i] = base_to_int(refseq[i]);
    for (unsigned int i = 0; i < readseq.size(); i++)
      read_base_ints[i] = base_to_int(readseq[i]);

    // Fill in requested subsection of the 3 matrices using dynamic programming
    for (int i = 1; i <= L2; i++){
      read_base = read_base_ints[i-1];

      for (int j = start_col; j <= end_col; j++){
	nindex    = i*(L1+1)+j;
	ntindex   = j-start_col+1;
	ref_base  = ref_base_ints[j-1];

	// Update M matrix (examine (i-1, j-1))
	oindex          = (i-1)*(L1+1)+(j-1);
	otindex         = j-start_col;
	s1              = M[oindex];
	s2              = Iref[oindex];
	s3              = Iread[oindex];
	M[nindex]       = bestIndex(s1, s2, s3, 
				    &(prev_M_tracks[otindex]), 
				    &(prev_Iref_tracks[otindex]), 
				    &(prev_Iread_tracks[otindex]), &c) + s[ref_base][read_base];
	traceM[nindex]  = c;
	if (c == 0)
	  M_tracks[ntindex] = prev_M_tracks[otindex].copy();
	else if (c == 1)
	  M_tracks[ntindex] = prev_Iref_tracks[otindex].copy();
	else if (c == 2)
	  M_tracks[ntindex] = prev_Iread_tracks[otindex].copy();
	else
	  printErrorAndDie("Invalid matrix type");
		
       	// Update Iref matrix (examine (i,j-1))
	oindex             = i*(L1+1) + (j-1);
	otindex            = j-start_col;
	s1                 = M[oindex]     - GAPOPEN;
	s2                 = Iref[oindex]  - GAPEXTEND;
	s3                 = Iread[oindex] - GAPOPEN;

	t1 = M_tracks[otindex].add(false);
	t2 = Iref_tracks[otindex].copy();
	t3 = Iread_tracks[otindex].add(false);
	Iref[nindex] = bestIndex(s1, s2, s3, &t1, &t2, &t3, &c);
	traceIref[nindex] = c;
	if (c == 0)
	  Iref_tracks[ntindex] = t1;
	else if (c == 1)
	  Iref_tracks[ntindex] = t2;
	else if (c == 2)
	  Iref_tracks[ntindex] = t3;
	else
	  printErrorAndDie("Invalid matrix type");

	// Update Iread matrix (examine (i-1,j))
	oindex              = (i-1)*(L1+1) + j;
	otindex             = j-start_col+1;
	s1                  = M[oindex]     - GAPOPEN;
	s2                  = Iref[oindex]  - GAPOPEN;
	s3                  = Iread[oindex] - GAPEXTEND;
	
	t1 = prev_M_tracks[otindex].add(false);
	t2 = prev_Iref_tracks[otindex].add(false);
	t3 = prev_Iread_tracks[otindex].copy();
	Iread[nindex] = bestIndex(s1, s2, s3, &t1, &t2, &t3, &c);
	traceIread[nindex] = c;
	if (c == 0)
	  Iread_tracks[ntindex] = t1;
	else if (c == 1)
	  Iread_tracks[ntindex] = t2;
	else if (c == 2)
	  Iread_tracks[ntindex] = t3;
	else
	  printErrorAndDie("Invalid matrix type");
      }

      // Move current tracks to previous tracks
      M_tracks.swap(prev_M_tracks);
      Iref_tracks.swap(prev_Iref_tracks);
      Iread_tracks.swap(prev_Iread_tracks);
    }
  }


  float ScoreAlignment(std::string& ref_seq_al, std::string& read_seq_al){
    if (ref_seq_al.size() != read_seq_al.size())
      printErrorAndDie("Alignment strings must have the same length");

    unsigned int start = 0; 
    while(start < read_seq_al.size() && read_seq_al[start] == '-')
      start++;
    int stop = read_seq_al.size()-1;
    while (stop >= 0 && read_seq_al[stop] == '-')
      stop--;

    float score   = 0;
    bool ref_gap  = false, read_gap = false;
    for(int i = start; i <= stop; i++){
      if (ref_seq_al[i] == '-'){
	if (ref_gap)
	  score -= GAPEXTEND;
	else {
	  score -= GAPOPEN;
	  ref_gap = true;
	}
      }
      else if (read_seq_al[i] == '-'){
	if(read_gap)
	  score -= GAPEXTEND;
	else {
	  score -= GAPOPEN;
	  read_gap = true;
	}
      }
      else { 
	ref_gap = read_gap = false;
	score  += s[base_to_int(ref_seq_al[i])][base_to_int(read_seq_al[i])];
      }
    }
    return score;
  }

  
  bool LeftAlign(const std::string& ref_seq, const std::string& read_seq,
		 std::string& ref_seq_al, std::string& read_seq_al,
		 float* score, std::vector<BamTools::CigarOp>& cigar_list) {
    int L1       = ref_seq.length();
    int L2       = read_seq.length();
    int mat_size = (L1+1)*(L2+1);

    // Scoring matrices:
    // M:     Ref and read bases aligned
    // Iref:  Ref base aligned with gap
    // Iread: Read base aligned with gap
    std::vector<float> M(mat_size), Iref(mat_size), Iread(mat_size);

    // Traceback matrices
    std::vector<int> traceM(mat_size), traceIref(mat_size), traceIread(mat_size);

    // Initialize matrices
    initMatrices(M, Iref, Iread, traceM, traceIref, traceIread, L1, L2);

    // Fill out scoring and traceback matrices using variant of NW algorithm
    nw_helper(M, Iref, Iread, traceM, traceIref, traceIread, ref_seq, read_seq);

    // Find the best ending point for the alignment
    float best_val;
    int best_col, best_type;
    findOptimalStop(L1, L2, M, Iref, Iread, best_val, best_col, best_type);
    *score = best_val;

    // Construct the alignment strings and CIGAR string using the traceback 
    // matrices and the optimal end position
    traceAlignment(best_col, best_type, L1, L2, traceM, traceIref, traceIread,
		   ref_seq, read_seq, ref_seq_al, read_seq_al, cigar_list);

    // Don't proceed if the read sequence extends past the reference boundaries
    if (cigar_list.front().Type == 'S' || cigar_list.back().Type == 'S')
      return false;
    
    // Determine start column index in matrix for optimal alignment
    unsigned int start_col = 0;
    while (start_col < read_seq_al.size() && read_seq_al[start_col] == '-')
      start_col++;
    start_col++;
    
    // Determine maximum number of indels
    int num_indels = 0;
    for (std::vector<BamTools::CigarOp>::iterator cigar_iter = cigar_list.begin(); cigar_iter != cigar_list.end(); cigar_iter++){
      if (cigar_iter->Type == 'I' || cigar_iter->Type == 'D')
	num_indels++;
    }

    // Check that the indel tracker can support this number of indels
    if (num_indels > IndelTracker::max_indels())
      return false;

    if (num_indels > 0){
      // Recalculate portion of matrices corresponding to optimal alignment
      // using dynamic programming that tracks indel locations
      left_align_helper(M, Iref, Iread, traceM, traceIref, traceIread,
			ref_seq, read_seq, start_col, best_col, num_indels);

      // Construct the alignment strings and CIGAR string using the fixed matrices
      traceAlignment(best_col, best_type, L1, L2, traceM, traceIref, traceIread,
		     ref_seq, read_seq, ref_seq_al, read_seq_al, cigar_list);
    }
    return true;
  }
}
