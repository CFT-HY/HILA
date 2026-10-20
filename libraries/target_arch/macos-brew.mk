# Platform specific makefile for macOS (Apple Silicon) with Homebrew
#
# Requires Homebrew: clang-21 (for hilapp), fftw-3, fftw-3-single, openmpi
#
# this is included from main.mk -file, which is in turn included from
# application makefile
#

# Define compiler: openmpi wrapper
CC := mpicxx
LD := mpicxx

# Define compilation flags
CPPSTD := c++17
CXXFLAGS_NOOPT := -x c++ --std=$(CPPSTD) -fno-rtti -Wno-nullability-completeness
ifndef DEBUG
	CXXFLAGS := $(CXXFLAGS_NOOPT) -O3 -mcpu=native
else
	CXXFLAGS := -g $(CXXFLAGS_NOOPT)
endif

# hilapp does not know the macOS SDK; give the std include dirs explicitly.
# Order: libc++ headers, clang builtin headers, SDK C headers
MACOS_SDK := $(shell xcrun --show-sdk-path)
HILAPP_LLVM_DIR := $(shell brew --prefix llvm@21)
CLANG_RESOURCE_DIR := $(shell $(HILAPP_LLVM_DIR)/bin/clang -print-resource-dir)
STD_INCLUDE_DIRS := -I$(MACOS_SDK)/usr/include/c++/v1 \
                    -I$(CLANG_RESOURCE_DIR)/include \
                    -I$(MACOS_SDK)/usr/include

# Homebrew fftw and openmpi
BREW_PREFIX := $(shell brew --prefix)
BREW_INCLUDE := -I$(BREW_PREFIX)/include -I$(shell brew --prefix open-mpi)/include


LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS := -L$(BREW_PREFIX)/lib

# These variables must be defined here
#
HILAPP_OPTS := $(STD_INCLUDE_DIRS) $(BREW_INCLUDE)
HILA_OPTS := -DDARWIN -I$(BREW_PREFIX)/include
