#ifndef HILAPP_SRCBUF_H
#define HILAPP_SRCBUF_H
using namespace clang;

/// New buffer interface.  Libtooling contains a  ReplaceText and RemoveText
/// in libtooling, but it is much easier to leave that to be original and
/// edit own copy
///
/// This is a simple buffer, keeps track of modifications to the original
/// text which remains in unmodified form.
/// This stores the original text in string buffer buf.  In addtion,
/// "vector<int> ext_ind" of length buf is allocated, and extents are
/// stored in "vector<string> extents"
/// ext_ind[i] == 1: nothing special
/// ext_ind[i] > 1: content of  extents[ext_ind[i]-2] is inserted before buf[i]
/// ext_ind[i] == 0: content of buf[i] skipped.
/// ext_ind[i] < -1: content of extents[-ext_ind[i]-2] is inserted before buf[i],
///                  and buf[i] skipped
///
/// Example:  buf="0123456789", ext_ind="1111211111", extents[0] = "cat"
///           would read as "0123cat456789"
///       or, if ext_ind="1111(-2)0011" otherwise as above reads "0123cat789"

class srcBuf {
  private:
    std::string buf;
    std::vector<int> ext_ind;
    std::vector<std::string> extents;
    std::vector<int> free_ext;
    std::string prependbuf;
    bool write_ok;
    bool modified;
    // std::vector<srcbuftoken> tokens;
    Rewriter *myRewriter;
    unsigned first_offset, full_length, original_size, true_size;

  public:
    srcBuf() {
        clear();
    }

    ~srcBuf() {
        clear();
    }

    srcBuf(Rewriter *R, Expr *E) {
        create(R, E);
    }
    srcBuf(Rewriter *R, Stmt *S) {
        create(R, S);
    }
    srcBuf(Rewriter *R, Decl *D) {
        create(R, D);
    }
    srcBuf(Rewriter *R, SourceRange sr) {
        create(R, sr);
    }

    int get_offset(SourceLocation s);

    void off() {
        write_ok = false;
    }
    void on() {
        write_ok = true;
    }

    /// This method returns the srcbuf index corresponding to sourcelocation
    int get_index(SourceLocation s) {
        int l = get_offset(s) - first_offset;
        assert(l >= 0 && l < true_size);
        return l;
    }

    /// returns the start and begin index of the sourcerange in std::pair
    std::pair<int, int> get_sourcerange_index(const SourceRange &sr);

    unsigned size() {
        return true_size;
    }

    void create(Rewriter *R, const SourceRange &sr);
    void create(Rewriter *R, Expr *e);
    void create(Rewriter *R, Stmt *s);
    void create(Rewriter *R, Decl *d);

    void clear();

    int get_index_range_size(int i1, int i2);

    // the mapped size of the range
    int get_sourcerange_size(const SourceRange &s);

    bool is_in_range(const SourceLocation sl);
    bool is_in_range(const SourceRange &r);

    // get string from index to mapped length
    std::string get_mapped(int index, int len);
    // get edited string originally from range
    std::string get(int i1, int i2);
    std::string get(SourceLocation s, int len);
    std::string get(const SourceRange &s);

    std::string dump();

    bool isOn();
    bool is_modified() {
        return modified;
    }

    char get_original(int i);
    int find_original(int idx, const char c);
    int find_original(SourceLocation start, const char c);
    int find_original(int idx, const std::string &s);
    int find_original(SourceLocation start, const std::string &c);

    int find_original_reverse(int idx, const char c);
    int find_original_reverse(SourceLocation start, const char c);

    int find_original_word(int idx, const std::string &s,
                           bool reverse = false); // finds full word
    int find_original_word(SourceLocation start, const std::string &s, bool reverse = false);

    /// give next or previous word or special char from buffer.  if idxp != nullptr, it
    /// contains the index
    std::string get_next_original_word(int idx, int *idxp = nullptr);
    std::string get_next_original_word(SourceLocation s, int *idxp = nullptr);
    std::string get_previous_original_word(int idx, int *idxp = nullptr);
    std::string get_previous_original_word(SourceLocation s, int *idxp = nullptr);

    bool is_edited(SourceLocation sl); // true if sl is edited (inserted, replaced or deleted)
    bool is_edited(int i); 
    
    bool is_extent(int i);

    void remove_extent(int i);

    std::string *get_extent_ptr(int i);

    // erase text between index values (note: not length!)
    // return value: index of next char
    int remove(int index1, int index2);
    int remove(const SourceRange &s);
    int remove(const CharSourceRange &s);
    int remove(Expr *E);
    // remove including possible comma before or after the range
    // useful for removing arguments
    int remove_with_comma(const SourceRange &s);
    int remove_semicolon_after(const SourceRange &s);
    int remove_semicolon_after(const Expr *E);

    // insert text - return is just insertion point + 1
    int insert(int i, const std::string &s, bool incl_before = false, bool indent = false);
    int insert(SourceLocation sl, const std::string &s, bool incl_before = false,
               bool indent = false);
    int insert(Expr *e, const std::string &s, bool incl_before = false, bool indent = false);
    int insert_after(SourceLocation sl, const std::string &s, bool incl_before = false,
                     bool indent = false);
    int insert_above(int i, const std::string &s, bool incl_before, bool do_indent);
    int insert_above(SourceLocation sl, const std::string &s, bool incl_before, bool do_indent);
    int insert_above(Expr *e, const std::string &s, bool incl_before, bool do_indent);

    int comment_line(int i);
    int comment_line(SourceLocation sl);
    int comment_line(Expr *e);
    int comment_range(int a, int b);

    // Find the end of previous statement and add after it
    int insert_before_stmt(int i, const std::string &s, bool incl_before, bool do_indent);
    int insert_before_stmt(SourceLocation sl, const std::string &s, bool incl_before,
                           bool do_indent);
    int insert_before_stmt(Expr *e, const std::string &s, bool incl_before, bool do_indent);

    // replace is a remove + insert pair, should write with a single operation
    // return: next element from remove
    int replace(int i1, int i2, const std::string &s);
    int replace(const SourceRange &r, const std::string &s);
    int replace(const CharSourceRange &r, const std::string &s);
    int replace(Expr *e, const std::string &s);

    void append(const std::string &s, bool do_indent);
    void prepend(const std::string &s, bool do_indent);

    void copy_from_range(srcBuf *src, SourceRange range);

    int replace_tokens(int start, int end, const std::vector<std::string> &a,
                       const std::vector<std::string> &b);
    int replace_tokens(SourceRange r, const std::vector<std::string> &a,
                       const std::vector<std::string> &b);
    int replace_token(int start, int end, const std::string &a, const std::string &b);
    int replace_tokens(const std::vector<std::string> &a, const std::vector<std::string> &b);
};

#endif // HILAPP_SRCBUF_H
