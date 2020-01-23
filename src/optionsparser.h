// -*- mode: c++ -*-
//===- OptionsParser.h - common options for clang tools -*- C++ -*-=====//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
//  This file implements the CommonOptionsParser class used to parse common
//  command-line options for clang tools, so that they can be run as separate
//  command-line applications with a consistent common interface for handling
//  compilation database and input files.
//
//  It provides a common subset of command-line options, common algorithm
//  for locating a compilation database and source files, and help messages
//  for the basic command-line interface.
//
//  It creates a CompilationDatabase and reads common command-line options.
//
//  This class uses the Clang Tooling infrastructure, see
//    http://clang.llvm.org/docs/HowToSetupToolingForLLVM.html
//  for details on setting it up with LLVM source tree.
//
//===----------------------------------------------------------------------===//
 
#ifndef OPTIONSPARSER_H
#define OPTIONSPARSER_H
 
#include "clang/Tooling/CompilationDatabase.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/Error.h"
 
namespace clang {
  namespace tooling {
  
    class OptionsParser {
    public:
      /// Parses command-line, initializes a compilation database.
      ///
      /// This constructor can change argc and argv contents, e.g. consume
      /// command-line options used for creating FixedCompilationDatabase.
      ///
      /// All options not belonging to \p Category become hidden.
      ///
      /// This constructor exits program in case of error.
      OptionsParser(int &argc, const char **argv,
                    llvm::cl::OptionCategory &Category,
                    const char *Overview = nullptr)
        : OptionsParser(argc, argv, Category, llvm::cl::OneOrMore,
                        Overview) {}
 
     
      OptionsParser(int &argc, const char **argv,
                    llvm::cl::OptionCategory &Category,
                    llvm::cl::NumOccurrencesFlag OccurrencesFlag,
                    const char *Overview = nullptr);
 
      /// Returns a reference to the loaded compilations database.
      CompilationDatabase &getCompilations() {
        return *Compilations;
      }
 
      /// Returns a list of source file paths to process.
      const std::vector<std::string> &getSourcePathList() const {
        return SourcePathList;
      }
  
    private:
      OptionsParser() = default;
 
      llvm::Error init(int &argc, const char **argv,
                       llvm::cl::OptionCategory &Category,
                       llvm::cl::NumOccurrencesFlag OccurrencesFlag,
                       const char *Overview);
 
      std::unique_ptr<CompilationDatabase> Compilations;
      std::vector<std::string> SourcePathList;
    };
 
 
  }  // namespace tooling
}  // namespace clang
 
#endif  // OPTIONSPARSER_H
