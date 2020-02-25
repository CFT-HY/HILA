//------------------------------------------------------------------------------
// Transformer tools to convert "lattice loops" into
// hardware-dependent "kernels".
//
// Uses Clang RecursiveASTVisitor and Rewriter
// interfaces
//
// Kari Rummukainen 2017-19
//
//------------------------------------------------------------------------------
#include <sstream>
#include <iostream>
#include <string>

#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"

#include "transformer.h"
#include "optionsparser.h"
#include "stringops.h"
#include "myastvisitor.h"
#include "specialization_db.h"

//definitions for global variables
ClassTemplateDecl * field_decl = nullptr; 
ClassTemplateDecl * field_storage_decl = nullptr;   
const std::string field_storage_type = "field_storage<";
const std::string field_type = "field<";
std::list<field_ref> field_ref_list = {};
std::list<field_info> field_info_list = {};
std::list<var_info> var_info_list = {};
std::list<var_decl> var_decl_list = {};
std::list<array_ref> array_ref_list = {};
std::list<special_function_call> special_function_call_list = {};
std::vector<Expr *> remove_expr_list = {};

unsigned state::skip_children = 0;
unsigned state::scope_level = 0;
bool state::in_loop_body = false;
bool state::accept_field_parity = false;
bool state::loop_found = false;
bool state::dump_ast_next = false;
bool state::compile_errors_occurred = false;
bool state::check_loop = false;
bool state::loop_function_next = false;

///definition of command line options
llvm::cl::OptionCategory TransformerCat(program_name);

llvm::cl::opt<bool> cmdline::dump_ast("dump-ast", llvm::cl::desc("Dump AST tree"),
          llvm::cl::cat(TransformerCat));

llvm::cl::opt<bool> cmdline::no_include("noincl",
             llvm::cl::desc("Do not insert \'#include\'-files (for debug)"),
             llvm::cl::cat(TransformerCat));

llvm::cl::opt<std::string> cmdline::dummy_def("D", 
            llvm::cl::value_desc("name"),
            llvm::cl::desc("Define name/symbol for preprocessor"),
            llvm::cl::cat(TransformerCat));

llvm::cl::opt<std::string> cmdline::dummy_incl("I", 
             llvm::cl::desc("Directory for include file search"),
             llvm::cl::value_desc("directory"),
             llvm::cl::cat(TransformerCat));

llvm::cl::opt<bool> cmdline::function_spec_no_inline("function-spec-no-inline",
                          llvm::cl::desc("Do not mark generated function specializations \"inline\""),
                          llvm::cl::cat(TransformerCat));

llvm::cl::opt<bool> cmdline::method_spec_no_inline("method-spec-no-inline",
                        llvm::cl::desc("Do not mark generated method specializations \"inline\""),
                        llvm::cl::cat(TransformerCat));
  
llvm::cl::opt<bool> cmdline::funcinfo("ident-functions",
           llvm::cl::desc("Comment function call types in output"),
           llvm::cl::cat(TransformerCat));
  
llvm::cl::opt<bool> cmdline::no_output("no-output",
            llvm::cl::desc("No output file, for syntax check"),
            llvm::cl::cat(TransformerCat));
  
llvm::cl::opt<bool> cmdline::syntax_only("syntax-only",
              llvm::cl::desc("Same as no-output"),
              llvm::cl::cat(TransformerCat));
  
llvm::cl::opt<std::string> cmdline::output_filename("o",
                  llvm::cl::desc("Output file (default: <file>.cpt, write to stdout: -o - "),
                  llvm::cl::value_desc("name"),
                  llvm::cl::cat(TransformerCat));


// List of targets that can be specified in command line arguments

llvm::cl::opt<bool> cmdline::kernel("target:vanilla-kernel",
         llvm::cl::desc("Generate kernels"),
         llvm::cl::cat(TransformerCat));
  
llvm::cl::opt<bool> cmdline::vanilla("target:vanilla",
          llvm::cl::desc("Generate loops in place"),
          llvm::cl::cat(TransformerCat));

llvm::cl::opt<bool> cmdline::CUDA("target:CUDA",
          llvm::cl::desc("Generate CUDA kernels"),
          llvm::cl::cat(TransformerCat));

llvm::cl::opt<bool> cmdline::AVX512("target:AVX512",
          llvm::cl::desc("Generate AVX512 vectorized loops"),
          llvm::cl::cat(TransformerCat));

llvm::cl::opt<bool> cmdline::AVX("target:AVX",
          llvm::cl::desc("Generate AVX vectorized loops"),
          llvm::cl::cat(TransformerCat));

llvm::cl::opt<bool> cmdline::SSE("target:SSE",
          llvm::cl::desc("Generate SSE vectorized loops"),
          llvm::cl::cat(TransformerCat));

llvm::cl::opt<int> cmdline::VECTORIZE("target:VECTORIZE",
          llvm::cl::desc("Generate vectorized loops with given vector size \n"
          "For example -target:VECTORIZE=32 is equivalent to -target:AVX"),
          llvm::cl::cat(TransformerCat));

llvm::cl::opt<bool> cmdline::openacc("target:openacc",
          llvm::cl::desc("Offload to GPU using openACC"),
          llvm::cl::cat(TransformerCat));


// Debug and Utility arguments

llvm::cl::opt<bool> cmdline::func_attribute("function-attributes",
         llvm::cl::desc("write pragmas/attributes to functions called from loops"),
         llvm::cl::cat(TransformerCat));

CompilerInstance *myCompilerInstance; //this is needed somewhere in the code
global_state global;
loop_parity_struct loop_parity;
codetype target;


/// Check command line arguments and set appropriate flags in target
void get_target_struct(codetype & target) {
  if (cmdline::CUDA) {
    target.CUDA = true;
  } else if (cmdline::openacc) {
    target.openacc = true;
  } else if (cmdline::AVX) {
    target.VECTORIZE = true;
    target.vector_size = 32;
  } else if (cmdline::AVX512) {
    target.VECTORIZE = true;
    target.vector_size = 64;
  } else if (cmdline::SSE) {
    target.VECTORIZE = true;
    target.vector_size = 16;
  } else if (cmdline::VECTORIZE) {
    target.VECTORIZE = true;
    target.vector_size = cmdline::VECTORIZE;
  }
}



/// Call the backend function for handling loop functions
void MyASTVisitor::backend_handle_loop_function(FunctionDecl *fd) {
  // we should mark the function, but it is not necessarily in the
  // main file buffer
  if (target.CUDA) {
    handle_loop_function_cuda(fd);
  } else if (target.openacc) {
    handle_loop_function_openacc(fd);
  } else if (target.VECTORIZE) {
    handle_loop_function_avx(fd);
  }
}

/// Call the backend function for generating loop code
std::string MyASTVisitor::backend_generate_code(Stmt *S, bool semi_at_end, srcBuf & loopBuf) {
  std::stringstream code;
  if( target.CUDA ){
    code << generate_code_cuda(S,semi_at_end,loopBuf);
  } else if( target.openacc ){
    code << generate_code_openacc(S,semi_at_end,loopBuf);
  } else if(target.VECTORIZE) {
    code << generate_code_avx(S,semi_at_end,loopBuf);
  } else {
    code << generate_code_cpu(S,semi_at_end,loopBuf);
  }
  return code.str();
}

/// Call the backend function for generating the field storage type
void MyASTVisitor::backend_generate_field_storage_type(std::string typestr){
  if(target.VECTORIZE){
    generate_field_storage_type_AVX(typestr);
  }
}




#if 0
// It is very hard to anything with pragmas in clang.  Requires modification of
// clang itself, libtooling is not sufficient
// Define a pragma handler for #pragma heLpp
// NOTE: This is executed before AST analysis
class heLppPragmaHandler : public PragmaHandler {
  public:
    heLppPragmaHandler() : PragmaHandler("transformer") { }
    void HandlePragma(Preprocessor &PP, PragmaIntroducerKind Introducer,
                      Token &Tok) {
      // Handle the pragma

      llvm::errs() << "Got the pragma! name " << getName() << " Token " << Tok.getName() << '\n';

      static Token tok_dump_ast;

      Token pragma_name = Tok;
      PP.Lex(Tok);  // lex next token
      if (Tok.is(tok::identifier)) {

        if (PP.getSpelling(Tok) == "_transformer_cmd_dump_ast_") {

          llvm::errs() << "Got " << PP.getSpelling(Tok) << '\n';
          tok_dump_ast.setIdentifierInfo(Tok.getIdentifierInfo());

        } else if (PP.getSpelling(Tok) == "dump_ast") {

          llvm::errs() << "Dumping ast for next cmd\n";
          std::vector<Token> tokenlist;
          Token t;
          // synthesize "extern long _transformer_cmd_dump_ast_;"
          t.startToken();
          t.setKind(tok::kw_extern);
          t.setLocation(pragma_name.getLocation());
          t.setLength(pragma_name.getLength());
          tokenlist.push_back(t);

          t.setKind(tok::kw_long);  //
          t.setLocation(pragma_name.getLocation());
          t.setLength(pragma_name.getLength());
          tokenlist.push_back(t);

          t.startToken();
          t.setIdentifierInfo(tok_dump_ast.getIdentifierInfo());  // _transformer_cmd_
          t.setKind(tok::identifier);
          t.setLocation(pragma_name.getLocation());
          tokenlist.push_back(t);

          t.startToken();
          t.setKind(tok::semi);
          t.setLocation(pragma_name.getLocation());
          tokenlist.push_back(t);

          auto TokenArray = llvm::make_unique<Token[]>(tokenlist.size());
          std::copy(tokenlist.begin(), tokenlist.end(), TokenArray.get());
          PP.EnterTokenStream(std::move(TokenArray), tokenlist.size(),
                              /*DisableMacroExpansion=*/false);
        }
      }

    return;
  }

};

static PragmaHandlerRegistry::Add<heLppPragmaHandler> Y("heLpp","heL pragma description");
#endif

reduction get_reduction_type(bool is_assign,
                             std::string & assignop,
                             var_info & vi) {
  if (is_assign && (!vi.is_loop_local)) {
    if (assignop == "+=") return reduction::SUM;
    if (assignop == "*=") return reduction::PRODUCT;
  }
  return reduction::NONE;
}

/////////////////////////////////////////////////////////////////////////////


// file_buffer_list stores the edited source of all files

struct file_buffer {
  srcBuf sbuf;
  FileID fid;
};

std::vector<file_buffer> file_buffer_list = {};

srcBuf * get_file_buffer(Rewriter & R, const FileID fid) {
  for (file_buffer & fb : file_buffer_list) {
    if (fb.fid == fid) return( &fb.sbuf );
  }
  // Now allocate and return new buffer
    
  file_buffer fb;
  fb.fid = fid;
  file_buffer_list.push_back(fb);
  SourceManager &SM = R.getSourceMgr();
  SourceRange r(SM.getLocForStartOfFile(fid),SM.getLocForEndOfFile(fid));

  llvm::errs() << "Create buf for file "
               << SM.getFilename(SM.getLocForStartOfFile(fid)) << '\n';
  
  file_buffer_list.back().sbuf.create( &R, r );
  return( &file_buffer_list.back().sbuf );
}

// Implementation of the ASTConsumer interface for reading an AST produced
// by the Clang parser.
class MyASTConsumer : public ASTConsumer {
public:
  MyASTConsumer(Rewriter &R, ASTContext *C) : Visitor(R,C) { Rewriterp = &R; }


  // HandleTranslationUnit is called after the AST for the whole TU is completed
  // Need to use this interface to ensure that specializations are present
  virtual void HandleTranslationUnit(ASTContext & ctx) override {

    SourceManager &SM = ctx.getSourceManager();
    TranslationUnitDecl *tud = ctx.getTranslationUnitDecl();
    //tud->dump();

    // Traverse each declaration in the translation unit
    for (Decl* d : tud->decls() ) {

      // We only want to traverse user files. We achieve this by
      // first checking if the location of the declaration is not
      // in a system file or a virtual file.
      
      // Find the beginning location of the Declaration
      SourceLocation beginloc = d->getBeginLoc();
      // Macro source range doesn't point to it's actual location in the file
      // Find where the macro is called
      if (beginloc.isMacroID()) {
        Preprocessor &pp = myCompilerInstance->getPreprocessor();
        // Is there an easier way?
        CharSourceRange CSR = Visitor.getRewriter().getSourceMgr().getImmediateExpansionRange( beginloc );
        beginloc = CSR.getBegin();
      }

      // analyze only user files (these should be named)
      if (!SM.isInSystemHeader(beginloc) && SM.getFilename(beginloc) != "") {
      //if (!SM.isInSystemHeader(beginloc)) {
        // llvm::errs() << "Processing file " << SM.getFilename(beginloc) << "\n";
        // TODO: ensure that we go only through files which are needed!

        state::loop_found = false;

        // get our own file edit buffer (unless it exists)
        Visitor.set_writeBuf(SM.getFileID(beginloc));

        // save this for source location
        global.location.top = d->getSourceRange().getBegin();  
        global.location.bot = Visitor.getSourceLocationAtEndOfRange(d->getSourceRange());

        // Traverse the declaration using our AST visitor.
        Visitor.TraverseDecl(d);
        // llvm::errs() << "Dumping level " << i++ << "\n";
        if (cmdline::dump_ast) {
          if (!cmdline::no_include || SM.isInMainFile(beginloc))
            d->dump();
        }
      }
    }

    // check compile errors, as long as we have context -- use diagnostics engine
    auto & DE = ctx.getDiagnostics();
    state::compile_errors_occurred = DE.hasErrorOccurred();

  }


private:
  MyASTVisitor Visitor;
  Rewriter * Rewriterp;
};


//#define NEED_PP_CALLBACKS
#ifdef NEED_PP_CALLBACKS

class MyPPCallbacks : public PPCallbacks {
public:
  SourceLocation This_hashloc;
  std::string This_name;

  // This hook is called when #include (or #import) is processed
  void InclusionDirective(SourceLocation HashLoc,
                          const Token & IncludeTok,
                          StringRef FileName,
                          bool IsAngled,
                          CharSourceRange FilenameRange,
                          const FileEntry * File,
                          StringRef SearchPath,
                          StringRef RelativePath,
                          const Module * Imported,
                          SrcMgr::CharacteristicKind FileType) { }

  // This triggers when the preprocessor changes file (#include, exit from it)
  // Use this to track the chain of non-system include files
  void FileChanged(SourceLocation Loc, FileChangeReason Reason, SrcMgr::CharacteristicKind FileType,
                   FileID PrevFID) {
    SourceManager &SM = myCompilerInstance->getSourceManager();
    if (Reason == PPCallbacks::EnterFile &&
        FileType == SrcMgr::CharacteristicKind::C_User &&
        Loc.isValid() &&
        !SM.isInSystemHeader(Loc) &&
        !SM.isInMainFile(Loc) ) {

      llvm::errs() << "FILE CHANGED to " << SM.getFilename(Loc) << '\n';

    }
  }

  // This triggers when range is skipped due to #if (0) .. #endif
  //   void SourceRangeSkipped(SourceRange Range, SourceLocation endLoc) {
  //     // llvm::errs() << "RANGE skipped\n";
  //   }
};

#endif

// This struct will be used to keep track of #include-chains.

std::vector<FileID> file_id_list = {};

// Tiny utility to search for the list

bool search_fid(const FileID FID) {
  for (const FileID f : file_id_list) {
    if (f == FID) return true;
  }
  return false;
}

void set_fid_modified(const FileID FID) {
  if (search_fid(FID) == false) {
    // new file to be added
    file_id_list.push_back(FID);
    // llvm::errs() << "New file changed " << SM.getFileEntryForID(FID)->getName() << '\n';
  }
}


// For each source file provided to the tool, a new FrontendAction is created.
class MyFrontendAction : public ASTFrontendAction {
public:
  MyFrontendAction() {}

  virtual bool BeginSourceFileAction(CompilerInstance &CI) override {
    llvm::errs() << "** Starting operation on source file "+getCurrentFile()+"\n";

#ifdef NEED_PP_CALLBACKS
    // Insert preprocessor callback functions to the stream.  This enables
    // tracking included files, ranges etc.
    Preprocessor &pp = CI.getPreprocessor();
    std::unique_ptr<MyPPCallbacks> callbacks(new MyPPCallbacks());
    pp.addPPCallbacks(std::move(callbacks));
#endif

    global.main_file_name = getCurrentFile();

    file_id_list.clear();
    file_buffer_list.clear();
    field_decl = field_storage_decl = nullptr;

    return (true);
  }

  void insert_includes_to_file_buffer(FileID myFID) {
    // this is where to write
    srcBuf * buf = get_file_buffer(TheRewriter, myFID);
    // find files to be included
    SourceManager &SM = TheRewriter.getSourceMgr();
    for (FileID f : file_id_list) {
      SourceLocation IL = SM.getIncludeLoc(f);
      if (IL.isValid() && myFID == SM.getFileID(IL)) {
        // file f is included, but do #include there first
        insert_includes_to_file_buffer(f);

        // Find now '#include "file.h"' -stmt (no obv way!)
        SourceRange SR = SM.getExpansionRange(IL).getAsRange();
        SourceLocation e = SR.getEnd();
        SourceLocation b = SR.getBegin();

        // Find the end of the include statement, which is an end of line
        // TODO: do this on "buf" instead of original file data
        for (int i=1; i<100; i++) {
          const char * p = SM.getCharacterData(e.getLocWithOffset(i));
          if (p && *p == '\n') {
            SR = SourceRange(b,e.getLocWithOffset(i));
            break;
          }
        }

        // Get the filename (the part after the includestr)
        std::string includestr = TheRewriter.getRewrittenText(SR);

        // Find the start of the include statement
        e = SR.getEnd();
        for (int i=1; i<100; i++) {
          const char * p = SM.getCharacterData(b.getLocWithOffset(-i));
          if (p && *p == '#' && strncmp(p,"#include",8) == 0) {
            SR = SourceRange(b.getLocWithOffset(-i),e);
            break;
          }
        }

        // Remove "#include"
        buf->remove(SR);
        // TheRewriter.RemoveText(SR);

        // and finally insert
        // SourceRange r(SM.getLocForStartOfFile(f),SM.getLocForEndOfFile(f));
        srcBuf * buf_from = get_file_buffer(TheRewriter, f);
        // TheRewriter.InsertText(SR.getBegin(),
        buf->insert(SR.getBegin(),
                    "// start include "+includestr
                    + "---------------------------------\n"
                    + buf_from->dump() +
                    "// end include "+includestr
                    + "---------------------------------\n",
                    false);

      }
    }
  }

  // check and add FileID's for files in #include chains if needed
  void check_include_path(const FileID FID) {
    SourceManager &SM = TheRewriter.getSourceMgr();
    SourceLocation IL = SM.getIncludeLoc(FID);
    if (IL.isValid()) {
      FileID FID_up = SM.getFileID(IL);
      if (!search_fid(FID_up)) {
        file_id_list.push_back(FID_up);
        if (FID_up != SM.getMainFileID()) check_include_path(FID_up);
      }
    }
  }


  void EndSourceFileAction() override {
    SourceManager &SM = TheRewriter.getSourceMgr();
    llvm::errs() << "** EndSourceFileAction for: " << getCurrentFile() << '\n';
    // << SM.getFileEntryForID(SM.getMainFileID())->getName() << "\n";

    // Now emit rewritten buffers.

    if (!cmdline::no_output) {
      if (!cmdline::no_include) {

        // Modified files should be substituted on top of #include -directives
        // First, find buffers which are modified
        
        for ( file_buffer & fb : file_buffer_list ) {
          if (fb.sbuf.is_modified()) set_fid_modified(fb.fid);
        }

        // then, ensure that the full include chain is present in file_id_list
        // Use iterator here, because the list can grow!

        for ( FileID f : file_id_list ) {
          check_include_path(f);
        }

        insert_includes_to_file_buffer(SM.getMainFileID());
      }

      if (!state::compile_errors_occurred) {
        write_output_file( cmdline::output_filename,
                           get_file_buffer(TheRewriter,SM.getMainFileID())->dump() );

        if (cmdline::function_spec_no_inline || cmdline::method_spec_no_inline)
          write_specialization_db();
      } else {
        llvm::errs() << program_name << ": not writing output due to compile errors\n";
      }
    }

    file_buffer_list.clear();
    file_id_list.clear();

    // EndSourceFile();

  }

  std::unique_ptr<ASTConsumer> CreateASTConsumer(CompilerInstance &CI,
                                                 StringRef file) override {
    // llvm::errs() << "** Creating AST consumer for: " << file << "\n";
    TheRewriter.setSourceMgr(CI.getSourceManager(), CI.getLangOpts());
    myCompilerInstance = &CI;
    return llvm::make_unique<MyASTConsumer>(TheRewriter, &CI.getASTContext());
  }

private:
  Rewriter TheRewriter;
  // ASTContext  TheContext;
};







int main(int argc, const char **argv) {

  // TODO: clang CommandLine.cpp/.h has strange category and help
  // msg handling, should we get rid of it?

  // av takes over from argv
  const char **av = new const char *[argc+3];
  argc = rearrange_cmdline(argc, argv, av);
  av[argc++] = "-DTRANSFORMER";  // add global defn


  OptionsParser op(argc, av, TransformerCat);
  ClangTool Tool(op.getCompilations(), op.getSourcePathList());

  // We have command line args, possibly do something with them
  get_target_struct(target);
  if (cmdline::syntax_only) cmdline::no_output = true;

  // ClangTool::run accepts a FrontendActionFactory, which is then used to
  // create new objects implementing the FrontendAction interface. Here we use
  // the helper newFrontendActionFactory to create a default factory that will
  // return a new MyFrontendAction object every time.
  // To further customize this, we could create our own factory class.
  return Tool.run(newFrontendActionFactory<MyFrontendAction>().get());
}
