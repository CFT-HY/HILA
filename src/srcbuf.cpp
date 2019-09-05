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
  clear();
  myRewriter = R;
  int rsize = myRewriter->getRangeSize(sr);
  assert( rsize >= 0 && "srcbuf: RangeSize should be > 0");

  full_length  = rsize;
  first_offset = get_offset(sr.getBegin());
  // get_offset(s->getSourceRange().getEnd()) does not necessarily
  // give end location
    
  buf = myRewriter->getRewrittenText(sr) + "  ";  // couple of extra chars (1 needed for insert after)
  true_size = original_size = rsize;
      
  ext_ind.resize(original_size+2,1);  // default=1
  ext_ind[original_size] = ext_ind[original_size+1] = 0;  // skip these
  write_ok = true;
    
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
  prependbuf = "";
  // tokens.clear();
}


int srcBuf::get_index_range_size(int i1, int i2) {
  assert( i1 >= 0 && i1 <= i2 && i2 < true_size );
  
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
std::string srcBuf::get_mapped(int index, int len) {
  assert(index >= 0 && "srcbuf: offset error");

  std::string out(len,' ');

  int j = index;
  int i;
  for (i=0; i<len && j<true_size; i++,j++) {
    // llvm::errs() << " ext_ind at "<<j<<"  is "<<ext_ind[j] << '\n';
    while (j<true_size && ext_ind[j] == 0) j++;
    if (j < true_size) {
      if (ext_ind[j] == 1) out.at(i) = buf[j];
      else {
        const std::string &e = extents[abs(ext_ind[j])-2];
        for (int k=0; k<e.size() && i<len; k++,i++) out[i] = e[k];
        if (i<len && ext_ind[j] > 0) out[i] = buf[j];
        else i--;  // for-loop advances i too far
      }
    }
  }
  if (i < len) out.resize(i);
  return out;
}

// get edited string originally from range
std::string srcBuf::get(int i1, int i2) {
  return get_mapped(i1,get_index_range_size(i1,i2));
}

// buffer content from mapped sourcelocation: len is original "length"
std::string srcBuf::get(SourceLocation s, int len) {
  int i1 = get_index(s);
  return get(i1,i1+len-1);
}
  
// get edited string originally from range
std::string srcBuf::get(const SourceRange & s) {
  return get_mapped(get_index(s.getBegin()),get_sourcerange_size(s));
}

std::string srcBuf::dump() {
  if (prependbuf.size() > 0) return prependbuf + get_mapped(0,full_length);
  else return get_mapped(0,full_length);
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
  assert(i>=0 && i<true_size);
  int e = abs(ext_ind[i]);
  if (e > 1) return true;
  return false;
}

void srcBuf::remove_extent(int i) {
  if (!write_ok) return;
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
// return next index
int srcBuf::remove(int index1, int index2) {
  assert(index1 >= 0 && index2 >= index1 && index2 < true_size);
  if (!write_ok) return index2+1;

  for (int i=index1; i<=index2; i++) {
    remove_extent(i);
    if (ext_ind[i] > 0) {
      full_length--;
      ext_ind[i] = 0;
    }
  }
  return index2+1;
}

int srcBuf::remove(const SourceRange &s) {
  int i=get_index(s.getBegin());
  return remove(i, i+myRewriter->getRangeSize(s)-1);
}

int srcBuf::remove(const CharSourceRange &s) {
  int i=get_index(s.getBegin());
  return remove(i, i+myRewriter->getRangeSize(s)-1);
}

int srcBuf::remove(Expr *E) {
  return remove(E->getSourceRange());
}

// wipe range with included comma before or after, if found.
// return the next element 

int srcBuf::remove_with_comma(const SourceRange &s) {
  int after  = remove(s);
  int before = get_index(s.getBegin()) - 1;
  int r = before;
  while (r >= 0 && std::isspace(buf[r])) r--;
  if (r >= 0 && get_original(r) == ',') {
    remove(r,before);
    return after;
  } else {
    r = after;
    while (r < true_size && std::isspace(buf[r])) r++;
    if (r < true_size && get_original(r) == ',') {
      return remove(after,r);
    }
  }
  // here comma not found
  return after;
}

// incl_before = true -> insert include before others at this location
int srcBuf::insert(int i, const std::string & s_in, bool incl_before, bool do_indent) {

  if (!write_ok) return i+1;
  
  if (i == original_size) {
    // Now append at end, OK
    true_size = original_size + 1;
    //  ext_ind[original_size] = 0;
  }

  assert( i>=0 && i<true_size && "srcBuf insert range error");
  const std::string * s = &s_in;
  std::string indented = "";
  
  // Is there newline in s?  If so, think about indent
  // for simplicity, we use the indent of this line in original text
  // Logic not perfect, but it's all cosmetics anyway
  if (do_indent && s_in.find('\n') != std::string::npos) {
    int j;
    for (j=i-1; j>=0 && buf[j] != '\n'; j--) ;
    j++;
    size_t k;
    for (k=j; k<i && std::isspace(buf[k]); k++) ;
    std::string indent = buf.substr(j,k-j);

    k = 0;
    while (k < s_in.length()) {
      size_t nl = s_in.find('\n',k);
      if (nl != std::string::npos) {
        indented.append(s_in.substr(k,nl-k+1));
        indented.append(indent);
        k = nl + 1;
      } else {
        indented.append(s_in.substr(k,std::string::npos));
        k = s_in.length();
      }
    }
    s = &indented;
  } 

  std::string * p = get_extent_ptr(i);
  
  if (p != nullptr) {
    if (incl_before) *p = *s + *p;
    else p->append(*s);
  } else {
    int j;
    if (free_ext.empty()) {
      extents.push_back(*s);
      j = extents.size()-1;
    } else {
      j = free_ext.back();
      free_ext.pop_back();
      extents[j] = *s;
    }
    if (ext_ind[i] > 0) ext_ind[i] = j+2;
    else ext_ind[i] = -(j+2);
  }
  full_length += s->length();
  return i+1;
}
        
int srcBuf::insert(SourceLocation sl, const std::string & s,
                   bool incl_before, bool do_indent) {
  return insert(get_index(sl), s, incl_before, do_indent);
}

int srcBuf::insert(Expr *e, const std::string & s, bool incl_before, bool do_indent) {
  return insert(e->getSourceRange().getBegin(), s, incl_before, do_indent);
}
  
// replace is a remove + insert pair, should write with a single operation
// return the index after
int srcBuf::replace( int i1, int i2, const std::string &s ) {
  remove(i1,i2);
  insert(i1,s,true,false);
  return i2+1;
}

int srcBuf::replace( const SourceRange & r, const std::string &s ) {
  int i = remove(r);
  insert(r.getBegin(),s,true,false);
  return i;
}

int srcBuf::replace( Expr *e, const std::string &s ) {
  return replace( e->getSourceRange(), s );
}

void srcBuf::append( const std::string &s, bool do_indent ) {
  insert( original_size, s, false, do_indent);
}

void srcBuf::prepend( const std::string &s, bool do_indent ) {
  if (!write_ok) return;
  prependbuf = s + prependbuf;
  full_length += s.length();
}

void srcBuf::copy_from_range(srcBuf * src, SourceRange range) {
  clear();
  create(src->myRewriter, range);
  // fill in using methods, no raw copy of structs
  //llvm::errs() << "Dump: " << dump();
  int srcind = src->get_index(range.getBegin());
  for (int i=0; i<true_size; i++) {
    int j = i+srcind;
    if (src->ext_ind[j] <= 0) remove(i,i);
    if (src->is_extent(j)) insert(i,*src->get_extent_ptr(j),true,false);
  }
  //llvm::errs() << "Dump again: " << dump();

}



// replace legal name tokens with this routine
// note- search done on "unedited" string
// return value: number of replacements
int srcBuf::replace_tokens(int start, int end,
                           const std::vector<std::string> &a,
                           const std::vector<std::string> &b) {

  if (!write_ok) return 0;
  
  assert(start >= 0 && end >= start && end < true_size);

  int replacements = 0;
  
  std::string tok;
  // walk through the string, form tokens
  int i=start;
  while (i<=end) {
    tok = "";
    // letters, _
    if (std::isalpha(buf[i]) || '_' == buf[i]) {
      int st = i;
      for ( ; i<=end && (std::isalnum(buf[i]) || '_' == buf[i]); i++)
        tok.push_back(buf[i]);

      // now got the token, search... note: a[i] == "" never matches
      for (int j=0; j<a.size(); j++) {
        if (tok.compare(a[j]) == 0) {
          // substutute
          replace(st,st+tok.size()-1,b[j]);
          replacements++;
          break;
        }
      }
    } else {
      if (i > 0 && buf[i] == '/' && buf[i-1] == '/') {
        // skip commented lines
        for (i++ ; i<=end && buf[i] != '\n'; i++) ;
      } else if (i > 0 && buf[i-1] == '/' && buf[i] == '*') {
        // c-stype comment
        i = buf.find("*/",i+1) + 2;
      } else if (buf[i] == '#') {
        // skip preproc lines #
        for (i++ ; i<=end && buf[i] != '\n'; i++) ;
      } else if (buf[i] == '"') {
        // don't mess with strings
        for (i++ ; i<=end && buf[i] != '"'; i++) if (buf[i] == '\\') i++;
      } else if (buf[i] == '\'') {
        // char constant
        if (buf[i+1] == '\\') i++;
        i+=3;
      } else {
        i++;
      }
    }   
  }
  return replacements;
}

int srcBuf::replace_tokens(SourceRange r,
                           const std::vector<std::string> &a,
                           const std::vector<std::string> &b) {
  
  int start = get_index(r.getBegin());
  int end   = start + myRewriter->getRangeSize(r) - 1;

  return replace_tokens(start,end,a,b);
}
  


int srcBuf::replace_token(int start, int end, const std::string & a, const std::string & b ) {
  std::vector<std::string> va = {}, vb = {};
  va.push_back(a);
  vb.push_back(b);
  
  return replace_tokens(start,end,va,vb);
}
  

