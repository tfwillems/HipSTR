##
## Makefile for STRImputer
##

## Default compilation flags.
## Override with:
##   make CXXFLAGS=XXXXX
CXXFLAGS= -O3 -g -D_FILE_OFFSET_BITS=64 -std=c++0x -DMACOSX
#CXXFLAGS= -O0 -g -D_FILE_OFFSET_BITS=64 -std=c++0x


## Source code files, add new files to this list
SRC = main.cpp

# For each CPP file, generate an object file
OBJ := $(SRC:.cpp=.o)

LIBDAI_ROOT=/Users/tfwillems/Downloads/libDAI-0.3.1

# -lgmp -lgmpxx needed for libDAI linking
LIBS = -L./ -lz -lm -lgmp -lgmpxx
INCLUDE = -I$(LIBDAI_ROOT)/include/ -I/usr/local/opt/boost149/include
LIBDAI_LIB= $(LIBDAI_ROOT)/lib/libdai.a
ARGWEAVER_LIB= argweaver/lib/libargweaver.a

.PHONY: all
all: str-imputer

# Clean the generated files
.PHONY: clean
clean:
	rm -f *.o *.d str-imputer

# Clean all compiled files, including bamtools/vcflib
.PHONY: clean-all
clean-all: clean

# The GNU Make trick to include the ".d" (dependencies) files.
# If the files don't exist, they will be re-generated, then included.
# If this causes problems with non-gnu make (e.g. on MacOS/FreeBSD), remove it.
include $(subst .cpp,.d,$(SRC))

# The resulting binary executable
str-imputer: $(OBJ) $(LIBDAI_LIB) $(ARGWEAVER_LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

# Build each object file independently
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

# Auto-Generate header dependencies for each CPP file.
%.d: %.cpp
	$(CXX) -c -MP -MD $(CXXFLAGS) $(INCLUDE) $< > $@
