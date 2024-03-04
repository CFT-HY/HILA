
/////////////////////////////////////////////////////////////////////////////////////
/// Some sourceloc utilities for hilapp
/// These are all in global namespace, to make usage easier
/////////////////////////////////////////////////////////////////////////////////////

#include "hilapp.h"

/// Get next character and sourcelocation, while skipping comments.
/// On line comments return the eol char

SourceLocation getNextLoc(const SourceManager &SM, SourceLocation sl, bool forward) {
    bool invalid = false;

    int dir;
    if (forward)
        dir = 1;
    else
        dir = -1;
    SourceLocation s = sl.getLocWithOffset(dir);
    const char *c = SM.getCharacterData(s, &invalid);

    // skip comments - only c-style backwards
    while ('/' == *c) {

        if (forward && '/' == *SM.getCharacterData(s.getLocWithOffset(1), &invalid)) {
            // a comment, skip the rest of line
            while (!invalid && *SM.getCharacterData(s, &invalid) != '\n')
                s = s.getLocWithOffset(1);
            c = SM.getCharacterData(s, &invalid);

        } else if ('*' == *SM.getCharacterData(s.getLocWithOffset(dir), &invalid)) {
            // c-style comment
            s = s.getLocWithOffset(2 * dir);
            while (!invalid && *SM.getCharacterData(s, &invalid) != '*' &&
                   *SM.getCharacterData(s.getLocWithOffset(dir), &invalid) != '/')
                s = s.getLocWithOffset(dir);
            s = s.getLocWithOffset(2 * dir);
            c = SM.getCharacterData(s, &invalid);
        } else
            break; // exit from here
    }

    return s;
}

char getChar(const SourceManager &SM, SourceLocation sl) {
    bool invalid = false;
    const char *c = SM.getCharacterData(sl, &invalid);
    if (invalid)
        return 0;
    else
        return *c;
}

// Find the location of the next searched for char.
SourceLocation findChar(const SourceManager &SM, SourceLocation sloc, char ct) {
    bool invalid = false;
    while (sloc.isValid()) {
        const char *c = SM.getCharacterData(sloc, &invalid);
        if (*c == ct)
            return sloc;
        sloc = getNextLoc(SM, sloc);
    }
    return sloc;
}

/// Skip paren expression following sl, points after the paren
/// If partype == '(', just balance par expressions (may contain strings).
/// If partype == '<', balance < > -parens (may contain () -parens, which are balanced
/// in turn) This last one is useful for template scanning
/// If partype == '[', does the same as for '<'.

SourceLocation skipParens(const SourceManager &SM, SourceLocation sl, const char partype) {

    assert((partype == '(' || partype == '<' || partype == '[') && "Unknown paren type in skipParens");

    char endpar;
    if (partype == '<')
        endpar = '>';
    else if (partype == '[')
        endpar = ']';
    else
        endpar = ')';

    while (sl.isValid() && getChar(SM, sl) != partype)
        sl = getNextLoc(SM, sl);

    sl = getNextLoc(SM, sl);
    char c = getChar(SM, sl);
    while (sl.isValid() && c != endpar) {

        if (c == partype) {
            sl = skipParens(SM, sl, partype);
        } else if (c == '(') {
            sl = skipParens(SM, sl);
        } else if (c == '"' || c == '\'') {
            sl = skipString(SM, sl);
        } else {
            // default option, next char loc
            sl = getNextLoc(SM, sl);
        }

        // and check next char
        c = getChar(SM, sl);
    }

    // successful exit here
    if (c == endpar)
        return getNextLoc(SM, sl);
    // this must be !isValid()
    return sl;
}

/// Skip " " or ' ' after location, returns the next location

SourceLocation skipString(const SourceManager &SM, SourceLocation sl) {

    char c;
    do {
        c = getChar(SM, sl);
        sl = getNextLoc(SM, sl);
    } while (sl.isValid() && c != '\'' && c != '"');

    // first single char
    if (c == '\'') {
        c = getChar(SM, sl);
        if (c == '\\')
            sl = sl.getLocWithOffset(1);
        sl = sl.getLocWithOffset(2); // jumps over the closing '
        return sl;
    }

    // now string " "
    int lev = 1;
    while (lev > 0 && sl.isValid()) {
        char c = getChar(SM, sl);
        if (c == '\\')
            sl = sl.getLocWithOffset(1);
        if (c == '"')
            lev = 0;
        sl = sl.getLocWithOffset(1);
    }

    return sl;
}

/// Get next word starting from sl -- if end is non-null, return the end of the
/// word string here (points to last char)

std::string getNextWord(const SourceManager &SM, SourceLocation sl, SourceLocation *end) {
    while (std::isspace(getChar(SM, sl)))
        sl = getNextLoc(SM, sl); // skip spaces

    std::string res;
    SourceLocation endloc = sl;
    char c = getChar(SM, sl);
    if (std::isalnum(c) || c == '_') {
        while (sl.isValid() && (std::isalnum(c) || c == '_')) {
            res.push_back(c);
            endloc = sl;
            sl = getNextLoc(SM, sl);
            c = getChar(SM, sl);
        }
    } else
        res.push_back(c);
    if (end != nullptr)
        *end = endloc;
    return res;
}

/// Return the text within the range, inclusive (note - skip comments, as usual)
std::string getRangeText(const SourceManager &SM, SourceLocation begin, SourceLocation end) {

    // need to be in the same file - and end needs to be after begin
    if (SM.getFileID(begin) != SM.getFileID(end) || begin > end)
        return "";

    std::string res;
    do {
        res.push_back(getChar(SM, begin));
        begin = getNextLoc(SM, begin);
    } while (begin <= end);
    return res;
}
