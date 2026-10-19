# Platform specific makefile for macOS (Apple Silicon) with MacPorts
#
# Requires MacPorts: clang-18 (for hilapp), fftw-3, fftw-3-single, openmpi
#
# this is included from main.mk -file, which is in turn included from
# application makefile
#

# Define compiler: MacPorts openmpi wrapper around Apple clang
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
HILAPP_LLVM_DIR := /opt/local/libexec/llvm-18
STD_INCLUDE_DIRS := -I$(MACOS_SDK)/usr/include/c++/v1 \
                    -I$(HILAPP_LLVM_DIR)/lib/clang/18/include \
                    -I$(MACOS_SDK)/usr/include

# MacPorts fftw and openmpi
MACPORTS_INCLUDE := -I/opt/local/include -I/opt/local/include/openmpi-mp

LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS := -L/opt/local/lib

# These variables must be defined here
#
HILAPP_OPTS := $(STD_INCLUDE_DIRS) $(MACPORTS_INCLUDE)
HILA_OPTS := -DDARWIN -I/opt/local/include
