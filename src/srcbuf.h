#ifndef TRANSFORMER_SRCBUF_H
#define TRANSFORMER_SRCBUF_H

// New buffer interface, due to weirdly buggy(?) ReplaceText and RemoveText
// in libtooling
// NOTE: It's probably due to the funny use of SourceRange in Expr and Stmt.
// Easy to make errors

struct srcbuftoken {
  Expr * srcExpr;
  size_t begin, len;
};

class srcBuf {
public:
  srcBuf() {
    buffer.clear();
    bufStmt = NULL;
  }
  srcBuf( Rewriter * R, Stmt *s ) {
    create( R, s );
  }
  ~srcBuf() {
    clear();
  }

  // this routine creates the buffer
  void create( Rewriter * R, Stmt *s ) {
    bufStmt = s;
    myRewriter = R;
    buffer = myRewriter->getRewrittenText(s->getSourceRange());
    SourceManager &SM = myRewriter->getSourceMgr();
    startOffset = SM.getFileOffset( s->getSourceRange().getBegin() );
    tokens.clear();

    
    // llvm::errs() << buffer << '\n';
    // llvm::errs() << " * Start offs " << startOffset << '\n';
  }

  void clear() {
    buffer.clear();
    bufStmt = NULL;
    tokens.clear();
  }
  
  // mark the location of the expression pointed by e
  unsigned markExpr( Expr *e ) {
    srcbuftoken t;

    size_t start = myRewriter->getSourceMgr().getFileOffset( e->getSourceRange().getBegin() );
    
    // llvm::errs() << " ** SrcRange: " << startOffset <<" " << start << '\n';

    t.begin = start - startOffset;
    t.len   = myRewriter->getRangeSize( e->getSourceRange() );
    t.srcExpr = e;

    tokens.push_back(t);

    // Check just in case - the tokens should be the same in this case

    std::string tok( myRewriter->getRewrittenText(e->getSourceRange()) );
    if ( tok.compare(buffer.substr(t.begin,t.len)) != 0 ) {
      llvm::errs() << " * Buffer marking error, src: " << tok
                   << "  buffer " << buffer.substr(t.begin,t.len) << '\n';
      exit(0);
    }
    // ret index
    return ( tokens.size() - 1 );
  }

  // insert string to position pos - if incl == true, add 
  void insert( size_t pos, const std::string &s, bool incl = false ) {
    assert( pos <= buffer.length() );
    // llvm::errs() << " -- Inserting substr: " << s << " at pos " << pos << '\n';
    unsigned len = s.length();
    buffer.insert( pos, s );
    for ( size_t i=0; i<tokens.size(); i++ ) {
      if (tokens[i].begin > pos || (!incl && tokens[i].begin == pos)) tokens[i].begin += len;
      else if (incl && tokens[i].begin == pos) tokens[i].len += len;
      else if (tokens[i].begin + tokens[i].len >= pos) tokens[i].len += len;
    }
  }

  void remove( size_t b, size_t len ) {
    // llvm::errs() << " -- Buffer: " << buffer << '\n';
    size_t e = b+len-1;
    assert( e < buffer.length() );
    // llvm::errs() << " -- Removing substr: " << buffer.substr(b, len) << '\n';
    buffer.erase( b, len );
    for ( size_t i=0; i<tokens.size(); i++ ) {
      size_t end = tokens[i].begin + tokens[i].len - 1;
      if (tokens[i].begin > e ) {
        tokens[i].begin -= len;
      } else if (tokens[i].begin >= b) {
        if (end > e) tokens[i].len -= e-tokens[i].begin+1;
        else tokens[i].len = 0;
        tokens[i].begin = b;
      } else {
        if (end >= e) tokens[i].len -= len;
        else if (end >= b) tokens[i].len = end-b+1;
      }
    }
  }

  // replace is a remove + insert pair, should write with a single operation
  void replace( size_t b, size_t len, const std::string &s ) {
    remove(b,len);
    insert(b,s,true);
  }

  void replace( unsigned ind, const std::string &s) {
    assert(ind < tokens.size());
    replace( tokens[ind].begin, tokens[ind].len, s );
  }

  std::string dump() {
    return buffer;
  }

  std::string substr( size_t b, size_t len ) {
    assert(b+len-1 < buffer.size());
    return buffer.substr(b, len );
  }
  
  std::string token( unsigned ind ) {
    assert(ind < tokens.size());
    return substr( tokens[ind].begin, tokens[ind].len );
  }
  
    
private:
  std::string buffer;
  Stmt *bufStmt;
  std::vector<srcbuftoken> tokens;
  Rewriter * myRewriter;
  size_t startOffset;
};


#endif // TRANSFORMER_SRCBUF_H
