// -*- mode: c++ -*-

#include <sstream>
#include <string>

#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Frontend/CompilerInstance.h"
// #include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"
#include "clang/Rewrite/Core/Rewriter.h"
//#include "llvm/Support/raw_ostream.h"

using namespace clang;
//using namespace clang::driver;
using namespace clang::tooling;


#include "srcbuf.h"
#include "myastvisitor.h"


int srcBuf::get_offset( SourceLocation s ) {
  SourceManager &SM = myRewriter->getSourceMgr();
  return( SM.getFileOffset(s) );
}

void srcBuf::create( Rewriter *R, const SourceRange &sr) {
  myRewriter = R;
  int rsize = myRewriter->getRangeSize(sr);
  assert( rsize >= 0 && "srcbuf: RangeSize should be > 0");

  full_length  = rsize;
  first_offset = get_offset(sr.getBegin());
  // TODO: check that this works as intended!
  // get_offset(s->getSourceRange().getEnd()) does not necessarily
  // give end location
    
  buf = myRewriter->getRewrittenText(sr) + " ";  // 1 extra char

  // llvm::errs() << "Got buf:  " << buf << '\n';
    
  ext_ind.clear();
  ext_ind.resize(buf.size(),1);  // default=1
  ext_ind.back() = 0;
  free_ext.clear();
    
  // tokens.clear();
}
  
void srcBuf::create( Rewriter * R, Expr *e ) {
  create( R, e->getSourceRange() );
}

void srcBuf::create( Rewriter * R, Stmt *s ) {
  create( R, s->getSourceRange() );
}

void srcBuf::create( Rewriter * R, Decl *d ) {
  create( R, d->getSourceRange() );
}
    
void srcBuf::clear() {
  buf.clear();
  ext_ind.clear();
  extents.clear();
  free_ext.clear();
  // tokens.clear();
}


int srcBuf::get_index_range_size(int i1, int i2) {
  assert( i2 >= i1 && i2 < buf.size() );
  
  int size = 0;
  for (int i=i1; i<=i2; i++) {
    if (ext_ind[i] > 0) size++;
    if (abs(ext_ind[i]) > 1) size += extents[abs(ext_ind[i])-2].size();
  }
  return size;
}
  
// the mapped size of the range
int srcBuf::get_sourcerange_size(const SourceRange & s) {
  int ind = get_index(s.getBegin());
  int len = myRewriter->getRangeSize( s );

  return get_index_range_size(ind,ind+len-1);
}

  // get buffer content from index to length len
std::string srcBuf::get(int index, int len) {
  assert(index >= 0 && len <= full_length - index
         && "srcbuf: offset error");

  std::string out;
  out.clear();
  out.resize(len,' ');

  int j = index;
  for (int i=0; i<len; i++,j++) {
    // llvm::errs() << " ext_ind at "<<j<<"  is "<<ext_ind[j] << '\n';
    while (ext_ind[j] == 0) j++;
    if (ext_ind[j] == 1) out[i] = buf[j];
    else {
      std::string &e = extents[abs(ext_ind[j])-2];
      for (int k=0; k<e.size() && i<len; k++,i++) out[i] = e[k];
      if (i<len && ext_ind[j] > 0) out[i] = buf[j];
      else i--;  // for-loop advances i too far
    }
  }
  return out;
}

// buffer content from mapped sourcelocation
std::string srcBuf::get(SourceLocation s, int len) {
  return get(get_index(s),len);
}
  
// get edited string originally from range
std::string srcBuf::get(const SourceRange & s) {
  return get(s.getBegin(),get_sourcerange_size(s));
}

// get edited string originally from range
std::string srcBuf::get_range(int i1, int i2) {
  return get(i1,get_index_range_size(i1,i2));
}

std::string srcBuf::dump() {
  return get(0,full_length);
}

bool srcBuf::isOn() { return buf.size() > 0; }

char srcBuf:: get_original(int i) {
  return buf.at(i);
}
  
int srcBuf::find_original(int idx, const char c) {
  std::string::size_type i = buf.find(c, idx);
  if (i == std::string::npos) return -1;
  else return i;
}
  
int srcBuf::find_original(SourceLocation start, const char c) {
  return find_original(get_index(start),c);
}

int srcBuf::find_original(int idx, const std::string &s) {
  std::string::size_type i = buf.find(s, idx);
  if (i == std::string::npos) return -1;
  else return i;
}

int srcBuf::find_original(SourceLocation start, const std::string &c) {
  return find_original(get_index(start),c);
}
            
bool srcBuf::is_extent(int i) {
  int e = abs(ext_ind[i]);
  if (e > 1) return true;
  return false;
}

void srcBuf::remove_extent(int i) {
  int e = abs(ext_ind[i]);
  if (e > 1) {
    full_length -= extents[e-2].size();
    extents[e-2] = "";
    free_ext.push_back(e-2);
    if (ext_ind[i] > 0) ext_ind[i] = 1;
    else ext_ind[i] = 0;
  }
}

std::string * srcBuf::get_extent_ptr(int i) {
  if (!is_extent(i)) return nullptr;
  return &extents[abs(ext_ind[i])-2];
}
    
// erase text between index values (note: not length!)
void srcBuf::remove(int index1, int index2) {
  for (int i=index1; i<=index2; i++) {
    remove_extent(i);
    if (ext_ind[i] > 0) {
      full_length--;
      ext_ind[i] = 0;
    }
  }
}

void srcBuf::remove(const SourceRange &s) {
  int i=get_index(s.getBegin());
  remove(i, i+myRewriter->getRangeSize(s)-1);
}

void srcBuf::remove(Expr *E) {
  remove(E->getSourceRange());
}

void srcBuf::insert(int i, const std::string & s, bool incl_before) {
  std::string * p = get_extent_ptr(i);
  if (p != nullptr) {
    if (incl_before) *p = s + *p;
    else p->append(s);
  } else {
    int j;
    if (free_ext.empty()) {
      extents.push_back(s);
      j = extents.size()-1;
    } else {
      j = free_ext.back();
      free_ext.pop_back();
      extents[j] = s;
    }
    if (ext_ind[i] > 0) ext_ind[i] = j+2;
    else ext_ind[i] = -(j+2);
  }
  full_length += s.length();
}
        
void srcBuf::insert(SourceLocation sl, const std::string & s,
                    bool incl_before) {
  insert(get_index(sl), s, incl_before);
}

void srcBuf::insert(Expr *e, const std::string & s, bool incl_before) {
  insert(e->getSourceRange().getBegin(), s, incl_before);
}
  
// replace is a remove + insert pair, should write with a single operation
void srcBuf::replace( int i1, int i2, const std::string &s ) {
  remove(i1,i2);
  insert(i1,s,true);
}

void srcBuf::replace( const SourceRange & r, const std::string &s ) {
  remove(r);
  insert(r.getBegin(),s,true);
}

void srcBuf::replace( Expr *e, const std::string &s ) {
  replace( e->getSourceRange(), s );
}


// replace legal name tokens with this routine
// note- search done on "unedited" string
void srcBuf::replace_tokens(SourceRange r,
                            const std::vector<std::string> &a,
                            const std::vector<std::string> &b) {
  
  int start = get_index(r.getBegin());
  int end   = start + myRewriter->getRangeSize(r) - 1;

  assert(start >= 0 && end >= start && end < buf.size());
  
  std::string tok;
  // walk through the string, form tokens
  int i=start;
  bool only_whitespace_line = true;
  while (i<=end) {
    tok = "";
    // letters, _
    if (std::isalpha(buf[i]) || '_' == buf[i]) {
      int st = i;
      only_whitespace_line = false;
      for ( ; i<=end && (std::isalnum(buf[i]) || '_' == buf[i]); i++)
        tok.push_back(buf[i]);

      // now got the token, search...
      for (int j=0; j<a.size(); j++) {
        if (tok.compare(a[j]) == 0) {
          // substutute
          replace(st,st+tok.size()-1,b[j]);
          break;
        }
      }
    } else {
      if (i > 0 && buf[i] == '/' && buf[i-1] == '/') {
        // skip commented lines
        for (i++ ; i<=end && buf[i] != '\n'; i++) ;
        only_whitespace_line = true;
      } else if (buf[i] == '#' && only_whitespace_line) {
        // skip preproc lines #
        for (i++ ; i<=end && buf[i] != '\n'; i++) ;
        only_whitespace_line = true;
      } else if (buf[i] == '"') {
        // don't mess with strings
        for (i++ ; i<=end && buf[i] != '"'; i++) if (buf[i] == '\\') i++;
        only_whitespace_line = false;
      } else if (buf[i] == '\'') {
        // char constant
        if (buf[i+1] == '\\') i++;
        i+=3;
        only_whitespace_line = false;
      } else if (buf[i] == '\n') {
        i++;
        only_whitespace_line = true;
      } else {
        if (!std::isspace(buf[i])) only_whitespace_line = false;
        i++;
      }
    }   
  }
}



