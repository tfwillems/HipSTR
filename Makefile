##
## Makefile for STRImputer
##

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXXFLAGS= -O3 -g -D_FILE_OFFSET_BITS=64 -std=c++0x -DMACOSX
#CXXFLAGS= -O0 -g -D_FILE_OFFSET_BITS=64 -std=c++0x

## Source code files, add new files to this list
SRC = main.cpp error.cpp factor_builder.cpp stutter_model.cpp snp_phasing_quality.cpp snp_tree.cpp

# For each CPP file, generate an object file
OBJ := $(SRC:.cpp=.o)

BAMTOOLS_ROOT=bamtools
LIBDAI_ROOT=/Users/tfwillems/Downloads/libDAI-0.3.1
VCFLIB_ROOT=vcflib

# -lgmp -lgmpxx needed for libDAI linking 
LIBS = -L./ -lz -lm -lgmp -lgmpxx -L$(BAMTOOLS_ROOT)/lib -L$(VCFLIB_ROOT)/tabixpp/ 
INCLUDE = -I$(LIBDAI_ROOT)/include/ -I/usr/local/opt/boost149/include -I$(VCFLIB_ROOT)/ -I$(BAMTOOLS_ROOT)/src
LIBDAI_LIB = $(LIBDAI_ROOT)/lib/libdai.a
ARGWEAVER_LIB = argweaver/lib/libargweaver.a
BAMTOOLS_LIB = $(BAMTOOLS_ROOT)/lib/libbamtools.a
VCFLIB_LIB = vcflib/libvcflib.a

.PHONY: all
all: str-imputer snp_tree_test vcf_test

# Clean the generated files
.PHONY: clean
clean:
	rm -f *.o *.d str-imputer snp_tree_test vcf_test

# Clean all compiled files, including bamtools/vcflib
.PHONY: clean-all
clean-all: clean

# The GNU Make trick to include the ".d" (dependencies) files.
# If the files don't exist, they will be re-generated, then included.
# If this causes problems with non-gnu make (e.g. on MacOS/FreeBSD), remove it.
include $(subst .cpp,.d,$(SRC))

# The resulting binary executable
str-imputer: $(OBJ) $(LIBDAI_LIB) $(ARGWEAVER_LIB) $(BAMTOOLS_LIB) $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Build each object file independently
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

# Auto-Generate header dependencies for each CPP file.
%.d: %.cpp
	$(CXX) -c -MP -MD $(CXXFLAGS) $(INCLUDE) $< > $@

snp_tree_test: snp_tree_test.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

vcf_test: vcf.cpp snp_tree.cpp  $(VCFLIB_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Rebuild VCFLIB if needed.
$(VCFLIB_LIB):
	git submodule update --init --recursive vcflib
	git submodule update --recursive vcflib
	cd vcflib && $(MAKE)

# Rebuild BAMTools if needed.
$(BAMTOOLS_LIB):
	git submodule update --init --recursive bamtools
	git submodule update --recursive bamtools
	( cd bamtools && mkdir build && cd build && cmake .. && $(MAKE) )
