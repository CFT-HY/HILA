// -*- mode: c++ -*-
#ifndef TRANSFORMER_SRCBUF_H
#define TRANSFORMER_SRCBUF_H

// New buffer interface.  Libtooling contains a  ReplaceText and RemoveText
// in libtooling
// NOTE: It's probably due to the funny use of SourceRange in Expr and Stmt.
// Easy to make errors

// This stores the original text in string buffer buf.  In addtion,
// "vector<int> ext_ind" of length buf is allocated, and extents are
// stored in "vector<string> extents"
// ext_ind[i] == 1: nothing special
// ext_ind[i] > 1: content of  extents[ext_ind[i]-2] is inserted before buf[i]
// ext_ind[i] == 0: content of buf[i] skipped.
// ext_ind[i] < -1: content of extents[-ext_ind[i]-2] is inserted before buf[i],
//                  and buf[i] skipped 
//
// Example:  buf="0123456789", eind="1111211111", extents[0] = "cat"
//           would read as "0123cat456789"
//       or, if eind="1111(-2)0011" otherwise as above reads "0123cat789"


// struct srcbuftoken {
//   Expr * srcExpr;
//   int begin;
//   int end;
// };


// node is Stmt, Decl or Expr - depending on the type of context

using namespace clang;

class srcBuf {
private:
  std::string buf;
  std::vector<int> ext_ind;
  std::vector<std::string> extents;
  std::vector<int> free_ext;
  std::string prependbuf;
  bool write_ok;
  // std::vector<srcbuftoken> tokens;
  Rewriter * myRewriter;
  unsigned first_offset, full_length, original_size, true_size;

public:
  srcBuf() { buf.clear(); }

  ~srcBuf() { clear(); }

  srcBuf( Rewriter *R, Expr *E ) { create(R,E); }
  srcBuf( Rewriter *R, Stmt *S ) { create(R,S); }
  srcBuf( Rewriter *R, Decl *D ) { create(R,D); }

  int get_offset( SourceLocation s );

  void off() { write_ok = false; }
  void on()  { write_ok = true; }
  
  int get_index( SourceLocation s ) { 
    int l = get_offset(s) - first_offset;
    assert(l>=0 && l<true_size);
    return l;
  }
  
  void create( Rewriter *R, const SourceRange &sr);
  void create( Rewriter * R, Expr *e );
  void create( Rewriter * R, Stmt *s );
  void create( Rewriter * R, Decl *d );
    
  void clear();

  int get_index_range_size(int i1, int i2);
  
  // the mapped size of the range
  int get_sourcerange_size(const SourceRange & s);

  // get string from index to mapped length
  std::string get_mapped(int index, int len);
  // get edited string originally from range
  std::string get(int i1, int i2);
  std::string get(SourceLocation s, int len);
  std::string get(const SourceRange & s);

  std::string dump();

  bool isOn();

  char get_original(int i);
  int find_original(int idx, const char c);
  int find_original(SourceLocation start, const char c);
  int find_original(int idx, const std::string &s);
  int find_original(SourceLocation start, const std::string &c);

  bool is_extent(int i);

  void remove_extent(int i);

  std::string * get_extent_ptr(int i);
    
  // erase text between index values (note: not length!)
  // return value: index of next char
  int remove(int index1, int index2);
  int remove(const SourceRange &s);
  int remove(const CharSourceRange &s);
  int remove(Expr *E);
  // remove including possible comma before or after the range
  // useful for removing arguments
  int remove_with_comma(const SourceRange &s); 

  
  // insert text - return is just insertion point + 1
  int insert(int i, const std::string & s, bool incl_before = false, bool indent = false);
  int insert(SourceLocation sl, const std::string & s,
             bool incl_before = false, bool indent = false);
  int insert(Expr *e, const std::string & s, bool incl_before = false, bool indent = false);
  int insert_above(int i, const std::string & s, bool incl_before, bool do_indent);
  int insert_above(SourceLocation sl, const std::string & s, bool incl_before, bool do_indent);
  int insert_above(Expr *e, const std::string & s, bool incl_before, bool do_indent);
  
  int comment_line(int i);
  int comment_line(SourceLocation sl);
  int comment_line(Expr *e);


  // replace is a remove + insert pair, should write with a single operation
  // return: next element from remove
  int replace( int i1, int i2, const std::string &s );
  int replace( const SourceRange & r, const std::string &s );
  int replace( const CharSourceRange & r, const std::string &s );
  int replace( Expr *e, const std::string &s );
  
  void append( const std::string &s, bool do_indent );
  void prepend( const std::string &s, bool do_indent );
  
  void copy_from_range(srcBuf * src, SourceRange range);
  

  int replace_tokens(int start, int end,
                     const std::vector<std::string> &a,
                     const std::vector<std::string> &b);
  int replace_tokens(SourceRange r,
                     const std::vector<std::string> &a,
                     const std::vector<std::string> &b);    
  int replace_token(int start, int end, const std::string & a, const std::string & b );
  int replace_tokens(const std::vector<std::string> &a,
                     const std::vector<std::string> &b);

};


#endif // TRANSFORMER_SRCBUF_H
