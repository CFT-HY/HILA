#%-----------------------------------------------------------------------------
#% Makefile for building the "hilapp" hila lattice preprocessor.
#%
#% This version builds hilapp from "system level" distribution
#%
#% Use:
#%    make help        - print this help text
#%    make [-j N]      - build hilapp to directory hilapp/build
#%    make install     - build and move hilapp to hilapp/bin
#%    make static      - build a statically compiled hilapp.  This can
#%                       be copied to systems without llvm/clang libraries
#%                       (assuming architectures are binary compatible)
#%    make clean       - clean up
#%    make [..] LLVM_VERSION=N   - build using llvm version 
#%                       (assumes llvm-config-N exists)
#%-----------------------------------------------------------------------------
#
# The Makefile variables may be customized, depending on where
# and how you built LLVM & Clang. They can be overridden by setting them on the
# make command line: "make VARNAME=VALUE", etc.

# If you compiled LLVM/CLANG yourself, set 
# LLVM_SRC_PATH to the path to the root of the checked out source code. This
# directory should contain the configure script, the include/ and lib/
# directories of LLVM, Clang in tools/clang/, etc.
# LLVM_SRC_PATH := /usr/lib/llvm-8

# If it is a system distribution (linux), then llvm-config should be installed:

ifdef LLVM_VERSION 
	LLVM_PATH := $(shell llvm-config-$(LLVM_VERSION) --prefix)
else
	LLVM_PATH := $(shell llvm-config --prefix)
endif

LLVM_BIN_PATH := $(LLVM_PATH)/bin

$(info --------------------------------------------------------------)
$(info Type "make help" for help on options)
$(info Using LLVM_PATH = $(LLVM_PATH))
$(info Using LLVM_BIN_PATH = $(LLVM_BIN_PATH))
$(info --------------------------------------------------------------)


# Main part of the Makefile

# CXX has to be a modern C++ compiler that supports C++17. gcc 8 and
# higher or Clang 8 and higher are recommended.  In modern systems either of
# these should work

# CXX := g++
# CXX := clang++
CXX := $(LLVM_BIN_PATH)/clang++

CXXFLAGS := -fno-rtti -std=c++17 -g 

# Directories where stuff is found and built
SRCDIR := src
BUILDDIR := build
BINDIR := bin

# Set some flags for compile
# thse should be present also for own build for llvm/clang
LLVM_CXXFLAGS := `$(LLVM_BIN_PATH)/llvm-config --cxxflags`

###############################################################################
# Select non-static (standard) or static compile 
# Static is useful for portability, the binary hilapp can be copied directly to
# platforms without full llvm/clang development environment.  This unfortunately
# includes many supercomputing environments.
#
# It may be that hilapp does not find all include files needed.  

ifneq ($(findstring static,$(MAKECMDGOALS)),static)
# This branch for non-static compile

LLVM_LDFLAGS := `$(LLVM_BIN_PATH)/llvm-config --ldflags --libs --system-libs`
STATIC := 

else

# Now static target
# option --ignore-libllvm gives all libs as .a -files and libLLVM.so is not included

LLVM_LDFLAGS := `$(LLVM_BIN_PATH)/llvm-config --ldflags --link-static --libs --system-libs`
STATIC := -static

# and target "static" depends only on hilapp
static:  $(BUILDDIR)/hilapp ; @:

$(info --------------------------------------------------------------)
$(info Building statically linked program, portable )
$(info NOTE: compilation will give warnings about possible glibc incompatibility)
$(info issues when calling 'getpwnam' and 'getpwuid'. This seems to be harmless.)
$(info --------------------------------------------------------------)


endif
# End of dynamic/static linking part
###############################################################################

# These are required when compiling vs. a source distribution of Clang. For
# binary distributions llvm-config --cxxflags gives the right path.
#
# CLANG_INCLUDES := \
#	-I$(LLVM_SRC_PATH)/include \
#	-I$(LLVM_BUILD_PATH)/tools/clang/include \
#	-I$(LLVM_SRC_PATH)/tools/clang/include 



# List of Clang libraries to link. The proper -L will be provided by the
# call to llvm-config
# The whole thing (all -l libs) can be encapslated withing
# -Wl,--start-group
# -l ...
# -Wl,--end-group
# which makes correct order irrelevant.  However, this does not seem
# to work on a mac.  The order below hopefully works

CLANG_LIBS := \
	-lclangTooling \
	-lclangASTMatchers \
	-lclangFrontend \
	-lclangDriver \
	-lclangParse \
	-lclangSerialization \
	-lclangSema \
	-lclangEdit \
	-lclangAnalysis \
	-lclangToolingCore \
	-lclangAST \
	-lclangRewrite \
	-lclangLex \
	-lclangBasic \
	-lclangFrontendTool \
	-lclangRewriteFrontend \
	-lclangStaticAnalyzerFrontend \
	-lclangStaticAnalyzerCheckers \
	-lclangStaticAnalyzerCore


# Sources for this project


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
  $(BUILDDIR)/loop_function.o \
  $(BUILDDIR)/loop_function_visitor.o \
  $(BUILDDIR)/vectorization_info.o \
  $(BUILDDIR)/depends_on_site_visitor.o \
  $(BUILDDIR)/contains_random_visitor.o \
  $(BUILDDIR)/contains_loop_local_var_visitor.o \
  $(BUILDDIR)/function_contains_loop_visitor.o \
  $(BUILDDIR)/addrof_and_ref_visitor.o \
  $(BUILDDIR)/optionsparser.o \
  $(BUILDDIR)/specialization_db.o \
  $(BUILDDIR)/write_output.o




$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp ${HEADERS} Makefile
	@mkdir -p $(BUILDDIR)
	$(CXX) -c $(CXXFLAGS) $(STATIC) $(LLVM_CXXFLAGS) $(CLANG_INCLUDES) $< -o $@

$(BUILDDIR)/hilapp: ${OBJECTS} ${STATIC_INCLUDE}
	$(CXX) $(STATIC) $(CXXFLAGS) $(LLVM_CXXFLAGS) $(OBJECTS) \
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
	$(CXX) -c $(CXXFLAGS) $(STATIC) $(LLVM_CXXFLAGS) $(CLANG_INCLUDES) $(SHA_OPT) $< -o $@

endif


help:
	@sed -ne '/@sed/!s/#%//p' $(MAKEFILE_LIST)


.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/*
