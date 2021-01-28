#-------------------------------------------------------------------------------
# Makefile for building the "hilapp" hila lattice preprocessor.
#
# This version builds hilapp from "system level" distribution
#
# Modelled after the example clang tool Makefile by Eli Bendersky (eliben@gmail.com)
#    https://github.com/eliben/llvm-clang-samples
#
# Kari Rummukainen
# Jarno Rantaharju
# 
#-------------------------------------------------------------------------------
#
# The following variables will likely need to be customized, depending on where
# and how you built LLVM & Clang. They can be overridden by setting them on the
# make command line: "make VARNAME=VALUE", etc.

# If you compiled LLVM/CLANG yourself, set 
# LLVM_SRC_PATH to the path to the root of the checked out source code. This
# directory should contain the configure script, the include/ and lib/
# directories of LLVM, Clang in tools/clang/, etc.
# LLVM_SRC_PATH := /usr/lib/llvm-8

# If it is a "system library", then llvm-config should be installed:
LLVM_SRC_PATH := $(shell llvm-config --prefix)

# LLVM_BUILD_PATH is the directory in which you built LLVM - where you ran
# configure or cmake.

LLVM_BUILD_PATH := $(LLVM_SRC_PATH)

LLVM_BIN_PATH 	:= $(LLVM_BUILD_PATH)/bin

$(info -----------------------------------------------)
$(info Using LLVM_SRC_PATH = $(LLVM_SRC_PATH))
$(info Using LLVM_BUILD_PATH = $(LLVM_BUILD_PATH))
$(info Using LLVM_BIN_PATH = $(LLVM_BIN_PATH))
$(info -----------------------------------------------)

# CXX has to be a modern C++ compiler that supports C++17. gcc 8 and
# higher or Clang 8 and higher are recommended. Best of all, if you build LLVM
# from sources, use the same compiler you built LLVM with.
# Note: starting with release 3.7, llvm-config will inject flags that gcc may
# not support (for example '-Wcovered-switch-default'). If you run into this
# problem, build with CXX set to a modern clang++ binary instead of g++.

#CXX := g++-5
CXX := $(LLVM_BIN_PATH)/clang++
CXXFLAGS := -fno-rtti -std=c++17 -g

LLVM_CXXFLAGS := `$(LLVM_BIN_PATH)/llvm-config --cxxflags`
LLVM_LDFLAGS := `$(LLVM_BIN_PATH)/llvm-config --ldflags --libs --system-libs`


# These are required when compiling vs. a source distribution of Clang. For
# binary distributions llvm-config --cxxflags gives the right path.

CLANG_INCLUDES := \
	-I$(LLVM_SRC_PATH)/include \
	-I$(LLVM_SRC_PATH)/tools/clang/include \
	-I$(LLVM_BUILD_PATH)/tools/clang/include

# List of Clang libraries to link. The proper -L will be provided by the
# call to llvm-config
# Note that I'm using -Wl,--{start|end}-group around the Clang libs; this is
# because there are circular dependencies that make the correct order difficult
# to specify and maintain. The linker group options make the linking somewhat
# slower, but IMHO they're still perfectly fine for tools that link with Clang.

CLANG_LIBS := \
	-Wl,--start-group \
	-lclangAST \
	-lclangAnalysis \
	-lclangBasic \
	-lclangDriver \
	-lclangEdit \
	-lclangFrontend \
	-lclangFrontendTool \
	-lclangLex \
	-lclangParse \
	-lclangSema \
	-lclangEdit \
	-lclangASTMatchers \
	-lclangRewrite \
	-lclangRewriteFrontend \
	-lclangStaticAnalyzerFrontend \
	-lclangStaticAnalyzerCheckers \
	-lclangStaticAnalyzerCore \
	-lclangSerialization \
	-lclangToolingCore \
	-lclangTooling \
	-Wl,--end-group

# Internal paths in this project: where to find sources, and where to put
# build artifacts.

SRCDIR := src
BUILDDIR := build
BINDIR := bin

HEADERS = \
  $(SRCDIR)/hilapp.h \
  $(SRCDIR)/generalvisitor.h \
  $(SRCDIR)/toplevelvisitor.h \
  $(SRCDIR)/srcbuf.h \
  $(SRCDIR)/stringops.h \
  $(SRCDIR)/optionsparser.h \
  $(SRCDIR)/specialization_db.h


OBJECTS = \
  $(BUILDDIR)/hilapp.o \
  $(BUILDDIR)/stringops.o \
  $(BUILDDIR)/codegen.o \
  $(BUILDDIR)/codegen_cpu.o \
  $(BUILDDIR)/codegen_cuda.o \
  $(BUILDDIR)/codegen_avx_new.o \
  $(BUILDDIR)/codegen_openacc.o \
  $(BUILDDIR)/srcbuf.o \
  $(BUILDDIR)/generalvisitor.o \
  $(BUILDDIR)/toplevelvisitor.o \
  $(BUILDDIR)/loop_function_new.o \
  $(BUILDDIR)/loop_function_visitor.o \
  $(BUILDDIR)/vectorization_info.o \
  $(BUILDDIR)/depends_on_site.o \
  $(BUILDDIR)/contains_random.o \
  $(BUILDDIR)/function_contains_loop_checker.o \
  $(BUILDDIR)/addrof_and_ref_checker.o \
  $(BUILDDIR)/optionsparser.o \
  $(BUILDDIR)/specialization_db.o \
  $(BUILDDIR)/write_output.o



$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp ${HEADERS} Makefile
	@mkdir -p $(BUILDDIR)
	$(CXX) -c $(CXXFLAGS) $(LLVM_CXXFLAGS) $(CLANG_INCLUDES) $< -o $@

#$(BUILDDIR)/hilapp: ${OBJECTS}
#	$(CXX) $(CXXFLAGS) $(LLVM_CXXFLAGS) $(OBJECTS) \
#		-static $(CLANG_LIBS) $(LLVM_LDFLAGS) -o $@

$(BUILDDIR)/hilapp: ${OBJECTS}
	$(CXX) $(CXXFLAGS) $(LLVM_CXXFLAGS) $(OBJECTS) \
	$(CLANG_LIBS) $(LLVM_LDFLAGS) -o $@

install: $(BUILDDIR)/hilapp
	@mkdir -p $(BINDIR)
	mv $(BUILDDIR)/hilapp $(BINDIR)

# Insert git sha value to the code
GIT_SHA := $(shell git rev-parse --short=8 HEAD)

ifneq "$(GIT_SHA)" "" 

SHA_OPT := -DGIT_SHA_VALUE=$(GIT_SHA)
GIT_SHA_FILE := build/.git_sha_number_$(GIT_SHA)

# Force recompilation if git number has changed

$(GIT_SHA_FILE):
	-rm -f build/.git_sha_number_*
	touch $(GIT_SHA_FILE)

# make separate rule for stringops
$(BUILDDIR)/stringops.o : $(SRCDIR)/stringops.cpp $(HEADERS) Makefile $(GIT_SHA_FILE)
	@mkdir -p $(BUILDDIR)
	$(CXX) -c $(CXXFLAGS) $(LLVM_CXXFLAGS) $(CLANG_INCLUDES) $(SHA_OPT) $< -o $@

endif




.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/* *.dot test/*.pyc test/__pycache__
