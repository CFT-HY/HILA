#ifndef HILAPP_VECTOR_H_
#define HILAPP_VECTOR_H_

#include "../defs.h"

/////////////////////////////////////////////////////////
/// This header introduces vector classes vec4d etc. for hilapp
/// Needed because (statically compiled) hilapp cannot use gcc avx header files
/// (these are understandable only by gcc)
/////////////////////////////////////////////////////////

class Vec4i {
  private:
    int v[4];

  public:
    void store(int *) const;
    Vec4i &load(const int *);
    void insert(int, int);
};

class Vec4d {
  private:
    double v[4];

  public:
    void store(double *) const;
    Vec4d &load(const double *);
};

class Vec4q {
  private:
    int64_t v[4];

  public:
    void store(int64_t *) const;
    Vec4q &load(const int64_t *);
};

class Vec8f {
  private:
    float v[8];

  public:
    void store(float *) const;
    Vec8f &load(const float *);
};

class Vec8i {
  private:
    int v[8];

  public:
    void store(int *) const;
    Vec8i &load(const int *);
    void insert(int, int);
};

class Vec8d {
  private:
    double v[8];

  public:
    void store(double *) const;
    Vec8d &load(const double *);
};

class Vec8q {
  private:
    int64_t v[8];

  public:
    void store(int64_t *) const;
    Vec8q &load(const int64_t *);
};

class Vec16f {
  private:
    float v[16];

  public:
    void store(float *) const;
    Vec16f &load(const float *);
};

class Vec16i {
  private:
    int v[16];

  public:
    void store(int *) const;
    Vec16i &load(const int *);
    void insert(int, int);
};

Vec4i operator*(const Vec4i, const Vec4i);
Vec4d operator*(const Vec4d, const Vec4d);
Vec4q operator*(const Vec4q, const Vec4q);
Vec8f operator*(const Vec8f, const Vec8f);
Vec8i operator*(const Vec8i, const Vec8i);
Vec8d operator*(const Vec8d, const Vec8d);
Vec8q operator*(const Vec8q, const Vec8q);
Vec16f operator*(const Vec16f, const Vec16f);
Vec16i operator*(const Vec16i, const Vec16i);

#endif