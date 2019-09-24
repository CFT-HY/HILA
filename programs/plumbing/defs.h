#ifndef DEFS_H
#define DEFS_H

// Useful global definitions here -- this file should be included by (almost) all others


// TODO: default type real_t definition somewhere (makefile?)
using real_t = double;


// move these somewhere - use consts?
#define NDIM 4


// Direction and parity

#if   NDIM==4
enum direction { XUP, YUP, ZUP, TUP, TDOWN, ZDOWN, YDOWN, XDOWN, NDIRS };
#elif NDIM==3
enum direction { XUP, YUP, ZUP, ZDOWN, YDOWN, XDOWN, NDIRS };
#elif NDIM==2
enum direction { XUP, YUP, YDOWN, XDOWN, NDIRS };
#elif NDIM==1
enum direction { XUP, XDOWN, NDIRS };
#endif

static inline direction opp_dir(const direction d) { return static_cast<direction>(NDIRS - 1 - d); }

enum class parity : unsigned { none, even, odd, all, x };
// use here #define instead of const parity. Makes EVEN a protected symbol
const parity EVEN = parity::even;
const parity ODD  = parity::odd;
const parity ALL  = parity::all;
const parity X    = parity::x;

// turns EVEN <-> ODD, ALL remains.  X->none, none->none
static inline parity opp_parity(const parity p) {
  unsigned u = 0x3 & static_cast<unsigned>(p);
  return static_cast<parity>(0x3 & ((u<<1)|(u>>1)));
}


#endif
