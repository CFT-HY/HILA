DISTRO := $(shell cat /etc/*-release | grep -iE "^(ID=)" | cut -f2 -d"=")

$(info --------------------------------------------------------------)
$(info We are on $(DISTRO) distribution)
$(info --------------------------------------------------------------)

## Then, find out the llvm version, either from LLVM_VERSION variable or
## snooping llvm-config-NN, which must be along PATH

ifdef LLVM_VERSION 

ifeq (, $(shell which llvm-config-$(LLVM_VERSION)))
$(error llvm-config-$(LLVM_VERSION) not found along PATH!)
else
LLVM_PATH := $(shell llvm-config-$(LLVM_VERSION) --prefix)
endif

else   
# Now LLVM_VERSION not set
# Try plain llvm-config
$(info LLVM_VERSION not set, trying llvm-config)
ifneq (, $(shell which llvm-config))
LLVM_PATH := $(shell llvm-config --prefix)
else

# llvm-config not found, next step...
$(info "llvm-config" not found along PATH, trying version from clang.. )
VERSION_FROM_CLANG := $(shell clang --version | head -1 | awk 'BEGIN{FS="version *"} {print $$2}' | awk 'BEGIN{FS="."} {print $$1}' )
ifneq (, $(shell which llvm-config-$(VERSION_FROM_CLANG)))
LLVM_PATH := $(shell llvm-config-$(VERSION_FROM_CLANG) --prefix)
$(info clang version $(VERSION_FROM_CLANG))
else

$(info "clang" failed, trying to glob llvm-config-NN )
LLVM_CONFIG := $(shell ls /*/llvm-config-* | sort | tail -1)
ifeq (, $(shell which $(LLVM_CONFIG)))
$(error No llvm-config-NN found along PATH)
else
LLVM_PATH := $(shell $(LLVM_CONFIG) --prefix)
$(info using $(LLVM_CONFIG))
endif

endif
endif

# get the llvm version number if not set before
LLVM_VERSION = $(shell $(LLVM_PATH)/bin/llvm-config --version | cut -d. -f1)

endif

## Now LLVM_VERSION is set

LLVM_BIN_PATH := $(LLVM_PATH)/bin



CC= $(LLVM_BIN_PATH)/clang++
CFLAGS= -ggdb -std=c++17 -Wall -pedantic -Wextra -fsanitize=address

.DEFAULT_GOAL := build/test_idx_eval

CFLAGS += $(shell $(LLVM_BIN_PATH)/llvm-config --cflags) $(shell $(LLVM_BIN_PATH)/llvm-config --libs) $(shell $(LLVM_BIN_PATH)/llvm-config --ldflags)

CLANG_LIBS := \
	-Wl,--start-group \
	-lclangTooling \
	-lclangASTMatchers \
	-lclangFrontend \
	-lclangDriver \
	-lclangParse \
	-lclangSerialization \
	-lclangSema \
	-lclangEdit \
	-lclangAnalysis \
	-lclangTooling \
	-lclangToolingCore \
	-lclangRewrite \
	-lclangLex \
	-lclangBasic \
	-lclangFrontendTool \
	-lclangRewriteFrontend \
	-lclangStaticAnalyzerFrontend \
	-lclangStaticAnalyzerCheckers \
	-lclangStaticAnalyzerCore \
	-lclangAST \
	-lclangToolingSyntax
# Add new libs if clang version is > 14
compare = $(shell if [ $(1) -gt $(2) ] ; then echo gt ; else echo le ; fi)

ifeq ($(call compare,$(LLVM_VERSION),14),gt)
	CLANG_LIBS += \
		-lclangSupport \
		-lclangAPINotes 
endif

# And close group
CLANG_LIBS += -Wl,--end-group

ifeq ($(DISTRO), arch)
ifeq ($(STATIC),)
	CLANG_LIBS = -lclang-cpp
endif
endif


build/test_idx_eval: test_idx_eval.cpp tracked_var.h ../../src/idx_expr_evaluator.cpp
	@mkdir -p build
	$(CC) $(CFLAGS) -DTEST_IDX_EVAL test_idx_eval.cpp -o $@ $(CLANG_LIBS)

build/idx_eval_file: test_idx_eval.cpp ../../src/idx_expr_evaluator.cpp
	@mkdir -p build
	$(CC) $(CFLAGS) -DIDX_EVAL_FILE test_idx_eval.cpp -o $@ $(CLANG_LIBS)

STD_INCLUDE_DIRS := $(addprefix --extra-arg=-I, $(shell echo | g++ -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))
HILA_INCLUDES = \
  --extra-arg=-I$(MY_HILA_DIR) \
  --extra-arg=-I$(MY_HILA_DIR)/libraries \
  --extra-arg=-I$(MY_HILA_DIR)/libraries/plumbing \
  --extra-arg=-I$(MY_HILA_DIR)/libraries/datatypes \
  --extra-arg=-I$(MY_HILA_DIR)/libraries/tools 

# helper to run idx_eval_file with hila headers
.PHONY: run_test_input
run_test_input:
	./build/idx_eval_file --dot $(STD_INCLUDE_DIRS) $(HILA_INCLUDES) test_input.cpp


ifndef APPTAINER_BIN
APPTAINER_BIN:= apptainer
endif
.PHONY: run_test_sif
run_test_sif:
	$(APPTAINER_BIN) run --app idx_eval_file ./test.sif $(STD_INCLUDE_DIRS) $(HILA_INCLUDES) test_input.cpp