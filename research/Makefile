# Builds every C program under research/ into research/build/.
#
# These are one-off analysis programs, not part of the simulator; the top-level
# CMake build does not know about them.  cmake is not required here.
#
#   make            build everything
#   make build/r14a build one program
#   make clean

CC     ?= cc
CFLAGS ?= -O3 -march=native
LDLIBS ?= -lm

BUILD    := build
SOURCES  := $(wildcard cumulative/*.c diagonals/*.c audit/*.c)
BINARIES := $(patsubst %.c,$(BUILD)/%,$(notdir $(SOURCES)))

vpath %.c cumulative diagonals audit

.PHONY: all clean
all: $(BINARIES)

$(BUILD)/%: %.c | $(BUILD)
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)
