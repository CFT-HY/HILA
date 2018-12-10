#ifndef TRANSFORMER_SRCBUF_H
#define TRANSFORMER_SRCBUF_H

// New buffer interface, due to weirdly buggy(?) ReplaceText and RemoveText
// in libtooling

struct srcbuftoken {
  Expr * srcExpr;
  size_t begin, end;
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
    t.end   = t.begin + myRewriter->getRangeSize( e->getSourceRange() ) - 1;
    t.srcExpr = e;

    tokens.push_back(t);

    // Check just in case - the tokens should be the same in this case

    std::string tok( myRewriter->getRewrittenText(e->getSourceRange()) );
    if ( tok.compare(buffer.substr(t.begin,t.end-t.begin+1)) != 0 ) {
      llvm::errs() << " * Buffer marking error, src: " << tok
                   << "  buffer " << buffer.substr(t.begin,t.end-t.begin+1) << '\n';
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
      if (tokens[i].end >= pos) tokens[i].end += len;
    }
  }

  void remove( size_t b, size_t e ) {
    // llvm::errs() << " -- Buffer: " << buffer << '\n';
    assert( b <= e && e < buffer.length() );
    // llvm::errs() << " -- Removing substr: " << buffer.substr(b, e-b+1) << '\n';
    buffer.erase( b, e-b+1 );
    for ( size_t i=0; i<tokens.size(); i++ ) {
      if (tokens[i].begin > e ) {
        tokens[i].begin -= e-b+1;
        tokens[i].end   -= e-b+1;
      } else if (tokens[i].begin > b) {
        if (tokens[i].end > e ) tokens[i].end -= e-tokens[i].begin+1;
        else tokens[i].end = b;
        tokens[i].begin = b;
      } else {
        if (tokens[i].end > e) tokens[i].end -= e-b+1;
        else if (tokens[i].end > b) tokens[i].end = b;
      }
    }
  }

  // replace is a remove + insert pair, should write with a single operation
  void replace( size_t b, size_t e, const std::string &s ) {
    remove(b,e);
    insert(b,s,true);
  }

  void replace( unsigned ind, const std::string &s) {
    assert(ind < tokens.size());
    replace( tokens[ind].begin, tokens[ind].end, s );
  }

  std::string dump() {
    return buffer;
  }

  std::string substr( unsigned a, unsigned b ) {
    assert(a <= b && b < buffer.size());
    return buffer.substr(a, b-a+1 );
  }
  
  std::string token( unsigned ind ) {
    assert(ind < tokens.size());
    return substr( tokens[ind].begin, tokens[ind].end );
  }
  
    
private:
  std::string buffer;
  Stmt *bufStmt;
  std::vector<srcbuftoken> tokens;
  Rewriter * myRewriter;
  size_t startOffset;
};


#endif // TRANSFORMER_SRCBUF_H
