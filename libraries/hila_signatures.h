#ifndef HILA_SIGNATURES_H_
#define HILA_SIGNATURES_H_

/// Running HILA signature number - this is used to ensure that hilapp executable and the hila
/// framework source are mutually compatible.  Incompatibilities may arise if e.g. one updates hila
/// framework without recompiling hilapp.
///
/// The signature number is incremented when hila and/or hilapp are changed in
/// a manner which makes previous versions incompatible. This file is included both in compilation
/// of hila programs and in compilation of hilapp. When hilapp is compiled, it gets the current
/// symbol values. When hilapp is used to compile hila programs, it reads the then-current symbol
/// values and checks that:
///
///  1. in hilapp compiled HILA_SIGNATURE_NUMBER >= current MINIMUM_HILAPP_SIGNATURE
///  2. in hilapp compiled MINIMUM_HILA_SIGNATURE <= current HILA_SIGNATURE_NUMBER
///
/// If 1 is not true, an error message is printed suggesting recompilation of hilapp.
/// If 2 is not true, an error message suggests updating hila framework (git pull)
///
/// TYPICAL USE CASE: Increment all numbers to a new equal value when incompatible
/// change is made.  This is the default below.

#define HILA_SIGNATURE_NUMBER    6

#define MINIMUM_HILAPP_SIGNATURE 6
#define MINIMUM_HILA_SIGNATURE   6

#endif
