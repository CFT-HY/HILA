//===--- CommonOptionsParser.cpp - common options for clang tools ---------===//
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
 
#include "optionsparser.h"
#include "clang/Tooling/Tooling.h"
#include "llvm/Support/CommandLine.h"
 
using namespace clang::tooling;
using namespace llvm;
 
 
llvm::Error OptionsParser::init(
    int &argc, const char **argv, cl::OptionCategory &Category,
    llvm::cl::NumOccurrencesFlag OccurrencesFlag, const char *Overview) {
  
  static cl::opt<bool> Help("h", cl::desc("Alias for -help"), cl::Hidden,
       cl::sub(*cl::AllSubCommands));
 
  static cl::opt<std::string> BuildPath("p", cl::desc("Build path"),
       cl::Optional, cl::cat(Category),
       cl::sub(*cl::AllSubCommands));
 
  static cl::list<std::string> SourcePaths(
       cl::Positional, cl::desc("<source files>"), OccurrencesFlag,
       cl::cat(Category), cl::sub(*cl::AllSubCommands));
 
 
  cl::ResetAllOptionOccurrences();
 
  cl::HideUnrelatedOptions(Category);
 
  std::string ErrorMessage;
  Compilations =
    FixedCompilationDatabase::loadFromCommandLine(argc, argv, ErrorMessage);
  if (!ErrorMessage.empty())
    ErrorMessage.append("\n");
  llvm::raw_string_ostream OS(ErrorMessage);
  // Stop initializing if command-line option parsing failed.
  if (!cl::ParseCommandLineOptions(argc, argv, Overview, &OS)) {
    OS.flush();
    return llvm::make_error<llvm::StringError>(ErrorMessage,
                                               llvm::inconvertibleErrorCode());
  }
 
  cl::PrintOptionValues();
 
  SourcePathList = SourcePaths;
  if ((OccurrencesFlag == cl::ZeroOrMore || OccurrencesFlag == cl::Optional) &&
      SourcePathList.empty())
    return llvm::Error::success();
  if (!Compilations) {
    Compilations.reset(new FixedCompilationDatabase(".", std::vector<std::string>()));
  }
  
  return llvm::Error::success();
}
  
OptionsParser::OptionsParser(
      int &argc, const char **argv, cl::OptionCategory &Category,
      llvm::cl::NumOccurrencesFlag OccurrencesFlag, const char *Overview) {

  llvm::Error Err = init(argc, argv, Category, OccurrencesFlag, Overview);
  if (Err || SourcePathList.empty()) {
    if (Err) llvm::errs() << llvm::toString(std::move(Err));
    if (SourcePathList.empty()) llvm::errs() << argv[0] << ": no input files specified\n";
    exit(1);
  }
}

