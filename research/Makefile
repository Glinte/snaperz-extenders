# Builds every analysis program under research/ into research/build/.
#
# These are one-off analysis programs, not part of the simulator; the top-level
# CMake build does not know about them.  cmake is not required here.
#
# The two C++ ones under audit/ include ../../src/blockless.h, so they check the
# shipped counter rather than a copy of it.
#
#   make            build everything
#   make build/r14a build one program
#   make clean

CC      ?= cc
CXX     ?= c++
CFLAGS  ?= -O3 -march=native
CXXFLAGS ?= -O3 -march=native -std=c++17
LDLIBS  ?= -lm

BUILD     := build
CSOURCES  := $(wildcard cumulative/*.c diagonals/*.c audit/*.c)
CXXSOURCES := $(wildcard audit/*.cpp)
BINARIES  := $(patsubst %.c,$(BUILD)/%,$(notdir $(CSOURCES))) \
             $(patsubst %.cpp,$(BUILD)/%,$(notdir $(CXXSOURCES)))

vpath %.c cumulative diagonals audit
vpath %.cpp audit

.PHONY: all clean
all: $(BINARIES)

$(BUILD)/%: %.c | $(BUILD)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

$(BUILD)/%: %.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)
