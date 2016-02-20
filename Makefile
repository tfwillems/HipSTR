##
## Makefile for all executables
##

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXXFLAGS= -O3 -g -D_FILE_OFFSET_BITS=64 -std=c++0x -DMACOSX -pthread -pedantic #-Wunreachable-code -Weverything

## Source code files, add new files to this list
SRC_COMMON  = base_quality.cpp error.cpp region.cpp stringops.cpp seqio.cpp zalgorithm.cpp alignment_filters.cpp extract_indels.cpp mathops.cpp pcr_duplicates.cpp fastahack/Fasta.cpp fastahack/split.cpp
SRC_SIEVE   = filter_main.cpp filter_bams.cpp insert_size.cpp
SRC_HIPSTR  = hipstr_main.cpp bam_processor.cpp stutter_model.cpp snp_phasing_quality.cpp snp_tree.cpp em_stutter_genotyper.cpp seq_stutter_genotyper.cpp snp_bam_processor.cpp genotyper_bam_processor.cpp vcf_input.cpp read_pooler.cpp version.cpp process_timer.cpp
SRC_SEQALN  = SeqAlignment/AlignmentData.cpp SeqAlignment/HapAligner.cpp SeqAlignment/RepeatStutterInfo.cpp SeqAlignment/AlignmentModel.cpp SeqAlignment/AlignmentOps.cpp SeqAlignment/HapBlock.cpp SeqAlignment/NWNoRefEndPenalty.cpp SeqAlignment/Haplotype.cpp SeqAlignment/RepeatBlock.cpp SeqAlignment/HaplotypeGenerator.cpp SeqAlignment/HTMLCreator.cpp SeqAlignment/AlignmentViz.cpp SeqAlignment/AlignmentTraceback.cpp SeqAlignment/StutterAlignerClass.cpp
SRC_RNASEQ  = exploratory/filter_rnaseq.cpp exploratory/exon_info.cpp

# For each CPP file, generate an object file
OBJ_COMMON  := $(SRC_COMMON:.cpp=.o)
OBJ_SIEVE   := $(SRC_SIEVE:.cpp=.o)
OBJ_HIPSTR  := $(SRC_HIPSTR:.cpp=.o)
OBJ_SEQALN  := $(SRC_SEQALN:.cpp=.o)
OBJ_RNASEQ  := $(SRC_RNASEQ:.cpp=.o)

BAMTOOLS_ROOT=bamtools
VCFLIB_ROOT=vcflib

LIBS              = -L./ -lm -lhts -L$(BAMTOOLS_ROOT)/lib -L$(VCFLIB_ROOT)/tabixpp/ -Lvcflib/tabixpp/htslib/ -lz
INCLUDE           = -I$(BAMTOOLS_ROOT)/src -I$(VCFLIB_ROOT)/ -I$(VCFLIB_ROOT)/tabixpp/htslib/
BAMTOOLS_LIB      = $(BAMTOOLS_ROOT)/lib/libbamtools.a
VCFLIB_LIB        = vcflib/libvcflib.a
PHASED_BEAGLE_JAR = PhasedBEAGLE/PhasedBEAGLE.jar

.PHONY: all
all: version BamSieve HipSTR test/fast_ops_test test/haplotype_test test/read_vcf_alleles_test test/read_vcf_priors_test test/snp_tree_test test/vcf_snp_tree_test $(PHASED_BEAGLE_JAR) exploratory/RNASeq exploratory/Clipper exploratory/10X exploratory/Mapper
	rm version.cpp
	touch version.cpp

version:
	git describe --abbrev=7 --dirty --always --tags | awk '{print "#include \"version.h\""; print "const std::string VERSION = \""$$0"\";"}' > version.cpp

# Clean the generated files of the main project only (leave Bamtools/vcflib alone)
.PHONY: clean
clean:
	rm -f *.o *.d BamSieve HipSTR test/allele_expansion_test test/fast_ops_test test/haplotype_test test/read_vcf_alleles_test test/read_vcf_priors_test test/snp_tree_test test/vcf_snp_tree_test SeqAlignment/*.o ${PHASED_BEAGLE_JAR} exploratory/RNASeq exploratory/Clipper exploratory/Mapper exploratory/10X

# Clean all compiled files, including bamtools/vcflib
.PHONY: clean-all
clean-all: clean
	if test -d bamtools/build ; then \
		$(MAKE) -C bamtools/build clean ; \
		rm -rf bamtools/build ; \
	fi

# The GNU Make trick to include the ".d" (dependencies) files.
# If the files don't exist, they will be re-generated, then included.
# If this causes problems with non-gnu make (e.g. on MacOS/FreeBSD), remove it.
include $(subst .cpp,.d,$(SRC))

# The resulting binary executable

BamSieve: $(OBJ_COMMON) $(OBJ_SIEVE) $(BAMTOOLS_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

HipSTR: $(OBJ_COMMON) $(OBJ_HIPSTR) $(BAMTOOLS_LIB) $(VCFLIB_LIB) $(OBJ_SEQALN)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

exploratory/RNASeq: $(OBJ_COMMON) $(OBJ_RNASEQ) $(BAMTOOLS_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/haplotype_test: test/haplotype_test.cpp SeqAlignment/Haplotype.cpp SeqAlignment/HapBlock.cpp SeqAlignment/NWNoRefEndPenalty.cpp SeqAlignment/RepeatBlock.cpp error.cpp stringops.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/em_stutter_test: test/em_stutter_test.cpp em_stutter_genotyper.cpp genotyper_bam_processor.cpp error.cpp mathops.cpp stringops.cpp stutter_model.cpp $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/fast_ops_test: test/fast_ops_test.cpp mathops.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^

test/read_vcf_alleles_test: test/read_vcf_alleles_test.cpp error.cpp region.cpp vcf_input.cpp $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/read_vcf_priors_test: test/read_vcf_priors_test.cpp error.cpp region.cpp vcf_input.cpp $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/snp_tree_test: snp_tree.cpp error.cpp test/snp_tree_test.cpp $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

test/vcf_snp_tree_test: test/vcf_snp_tree_test.cpp error.cpp snp_tree.cpp $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)


exploratory/Clipper: exploratory/count_trimmed_bases.cpp error.cpp zalgorithm.cpp $(BAMTOOLS_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

exploratory/10X: base_quality.cpp exploratory/calc_10x_barcode_phasings.cpp error.cpp mathops.cpp snp_phasing_quality.cpp snp_tree.cpp $(BAMTOOLS_LIB) $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

exploratory/Mapper: error.cpp seqio.cpp stringops.cpp exploratory/mapping_efficiency.cpp $(BAMTOOLS_LIB) $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Build each object file independently
%.o: %.cpp $(BAMTOOLS_LIB)  $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

# Auto-Generate header dependencies for each CPP file.
%.d: %.cpp $(BAMTOOLS_LIB)  $(VCFLIB_LIB)
	$(CXX) -c -MP -MD $(CXXFLAGS) $(INCLUDE) $< > $@

# Rebuild BAMTools if needed
$(BAMTOOLS_LIB):
	git submodule update --init --recursive bamtools
	git submodule update --recursive bamtools
	( cd bamtools && mkdir build && cd build && cmake .. && $(MAKE) )

# Rebuild VCFLIB if needed               
$(VCFLIB_LIB):
	git submodule update --init --recursive vcflib
	git submodule update --recursive vcflib
	cd vcflib && $(MAKE)

# Rebuild PhasedBEAGLE if needed
$(PHASED_BEAGLE_JAR):
	git submodule update --init --recursive PhasedBEAGLE
	git submodule update --recursive PhasedBEAGLE
	cd PhasedBEAGLE && $(MAKE)
