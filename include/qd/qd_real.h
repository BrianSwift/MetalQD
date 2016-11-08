/*
 * include/qd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Quad-double precision (>= 212-bit significand) floating point arithmetic
 * package, written in ANSI C++, taking full advantage of operator overloading.
 * Uses similar techniques as that of David Bailey's double-double package 
 * and that of Jonathan Shewchuk's adaptive precision floating point 
 * arithmetic package.  See
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   http://www.cs.cmu.edu/~quake/robust.html
 *
 * for more details.
 *
 * Yozo Hida
 */
#ifndef _QD_QD_REAL_H
#define _QD_QD_REAL_H

// Configure for non-metal-ported behavior
#ifndef QDT
#define QDT double
#define QDT_double 1
#define QDT_float 0
#define QD_ASQ /* */
#define QD_GPU 0
#define QD_HAVE_SIMD 0
#endif


#if !QD_GPU
#include <iostream>
#include <string>
#include <limits>
#include <qd/qd_config.h>
#include <qd/dd_real.h>
#else
#include "qd_config.h"
// #include "qd/dd_real.h"
#endif

struct QD_API qd_real {
  QDT x[4];    /* The Components. */

  /* Eliminates any zeros in the middle component(s). */
  void zero_elim();
  void zero_elim(QD_ASQ QDT &e);

  void renorm();
  void renorm(QD_ASQ QDT &e);

  void quick_accum(QDT d, QD_ASQ QDT &e);
  void quick_prod_accum(QDT a, QDT b, QD_ASQ QDT &e);

  qd_real(QDT x0, QDT x1, QDT x2, QDT x3);
  explicit qd_real(const QD_ASQ QDT *xx);

#if !QD_GPU
  static const qd_real _2pi;
  static const qd_real _pi;
  static const qd_real _3pi4;
  static const qd_real _pi2;
  static const qd_real _pi4;
  static const qd_real _e;
  static const qd_real _log2;
  static const qd_real _log10;
  static const qd_real _nan;
  static const qd_real _inf;

  static const QDT _eps;
  static const QDT _min_normalized;
  static const qd_real _max;
  static const qd_real _safe_max;
  static const int _ndigits;
#endif

  qd_real();
#if !QD_GPU
  qd_real(const char *s);
#endif
#ifdef _QD_DD_REAL_H
  qd_real(const dd_real &dd);
#endif
  qd_real(QDT d);
  qd_real(int i);
    
// for double to quad float on CPU, !!! need to rework on metal GPU with double
#ifndef __METAL_VERSION__
#if QDT_float
    qd_real(double d);
#endif
#endif
    
#if QD_HAVE_SIMD
    qd_real(simd::float4);
#endif

  QDT operator[](int i) const;
  QD_ASQ QDT &operator[](int i);

#if !QD_GPU
  static void error(const char *msg);
#endif

  bool isnan() const;
  bool isfinite() const { return QD_ISFINITE(x[0]); }
  bool isinf() const { return QD_ISINF(x[0]); }

  static qd_real ieee_add(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
  static qd_real sloppy_add(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);

  QD_ASQ qd_real &operator+=(QDT a);
#ifdef _QD_DD_REAL_H
  qd_real &operator+=(const dd_real &a);
#endif
  QD_ASQ qd_real &operator+=(const QD_ASQ qd_real &a);

  QD_ASQ qd_real &operator-=(QDT a);
#ifdef _QD_DD_REAL_H
  qd_real &operator-=(const dd_real &a);
#endif
  QD_ASQ qd_real &operator-=(const QD_ASQ qd_real &a);

  static qd_real sloppy_mul(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
  static qd_real accurate_mul(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);

  QD_ASQ qd_real &operator*=(QDT a);
#ifdef _QD_DD_REAL_H
  qd_real &operator*=(const dd_real &a);
#endif
  QD_ASQ qd_real &operator*=(const QD_ASQ qd_real &a);

#ifdef _QD_DD_REAL_H
  static qd_real sloppy_div(const qd_real &a, const dd_real &b);
  static qd_real accurate_div(const qd_real &a, const dd_real &b);
#endif
  static qd_real sloppy_div(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
  static qd_real accurate_div(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);

  QD_ASQ qd_real &operator/=(QDT a);
#ifdef _QD_DD_REAL_H
  qd_real &operator/=(const dd_real &a);
#endif
  QD_ASQ qd_real &operator/=(const QD_ASQ qd_real &a);

  qd_real operator^(int n) const;

  qd_real operator-() const;

  QD_ASQ qd_real &operator=(QDT a);
#ifdef _QD_DD_REAL_H
  qd_real &operator=(const dd_real &a);
#endif
#if !QD_GPU
  qd_real &operator=(const char *s);
#endif
    
  bool is_zero() const;
  bool is_one() const;
  bool is_positive() const;
  bool is_negative() const;

  static qd_real rand(void);

#if !QD_GPU
  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits, 
      bool showpos = false, bool uppercase = false) const;
  std::string to_string(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  static int read(const char *s, qd_real &a);

  /* Debugging methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "", 
                 std::ostream &os = std::cerr) const;
#endif

  static qd_real debug_rand();

};

#if !QD_GPU
namespace std {
  template <>
  class numeric_limits<qd_real> : public numeric_limits<QDT> {
  public:
    inline static QDT epsilon() { return qd_real::_eps; }
    inline static QDT min() { return qd_real::_min_normalized; }
    inline static qd_real max() { return qd_real::_max; }
    inline static qd_real safe_max() { return qd_real::_safe_max; }
    static const int digits = 209;
    static const int digits10 = 62;
  };
}
#endif

#if !QD_GPU
QD_API qd_real polyeval(const qd_real *c, int n, const qd_real &x);
QD_API qd_real polyroot(const qd_real *c, int n, 
    const qd_real &x0, int max_iter = 64, QDT thresh = 0.0);
#endif

QD_API qd_real qdrand(void);
QD_ASQ QD_API qd_real sqrt(const qd_real QD_ASQ &a);

QD_API inline bool isnan(const QD_ASQ qd_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const QD_ASQ qd_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const QD_ASQ qd_real &a) { return a.isinf(); }

/* Computes  qd * d  where d is known to be a power of 2.
   This can be done component wise.                      */
QD_API qd_real mul_pwr2(const QD_ASQ qd_real &qd, QDT d);

QD_API qd_real operator+(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API qd_real operator+(const dd_real &a, const qd_real &b);
QD_API qd_real operator+(const qd_real &a, const dd_real &b);
#endif
QD_API qd_real operator+(const QD_ASQ qd_real &a, QDT b);
QD_API qd_real operator+(QDT a, const QD_ASQ qd_real &b);

QD_API qd_real operator-(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API qd_real operator-(const dd_real &a, const qd_real &b);
QD_API qd_real operator-(const qd_real &a, const dd_real &b);
#endif
QD_API qd_real operator-(const QD_ASQ qd_real &a, QDT b);
QD_API qd_real operator-(QDT a, const QD_ASQ qd_real &b);

QD_API qd_real operator*(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API qd_real operator*(const dd_real &a, const qd_real &b);
QD_API qd_real operator*(const qd_real &a, const dd_real &b);
#endif
QD_API qd_real operator*(const QD_ASQ qd_real &a, QDT b);
QD_API qd_real operator*(QDT a, const QD_ASQ qd_real &b);

QD_API qd_real operator/(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API qd_real operator/(const dd_real &a, const qd_real &b);
QD_API qd_real operator/(const qd_real &a, const dd_real &b);
#endif
QD_API qd_real operator/(const QD_ASQ qd_real &a, QDT b);
QD_API qd_real operator/(QDT a, const QD_ASQ qd_real &b);

QD_API qd_real sqr(const QD_ASQ qd_real &a);
QD_API qd_real sqrt(const QD_ASQ qd_real &a);
QD_API qd_real pow(const QD_ASQ qd_real &a, int n);
QD_API qd_real pow(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
QD_API qd_real npwr(const QD_ASQ qd_real &a, int n);

QD_API qd_real nroot(const QD_ASQ qd_real &a, int n);

QD_API qd_real rem(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
QD_API qd_real drem(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
QD_API qd_real divrem(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b, QD_ASQ qd_real &r);

#ifdef _QD_DD_REAL_H
dd_real to_dd_real(const qd_real &a);
#endif
QDT  to_double(const QD_ASQ qd_real &a);
int     to_int(const QD_ASQ qd_real &a);

#if QD_HAVE_SIMD
simd::float4 to_float4(const QD_ASQ qd_real &a);
#endif

QD_API bool operator==(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API bool operator==(const qd_real &a, const dd_real &b);
QD_API bool operator==(const dd_real &a, const qd_real &b);
#endif
QD_API bool operator==(QDT a, const QD_ASQ qd_real &b);
QD_API bool operator==(const QD_ASQ qd_real &a, QDT b);

QD_API bool operator<(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API bool operator<(const qd_real &a, const dd_real &b);
QD_API bool operator<(const dd_real &a, const qd_real &b);
#endif
QD_API bool operator<(QDT a, const QD_ASQ qd_real &b);
QD_API bool operator<(const QD_ASQ qd_real &a, QDT b);

QD_API bool operator>(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API bool operator>(const qd_real &a, const dd_real &b);
QD_API bool operator>(const dd_real &a, const qd_real &b);
#endif
QD_API bool operator>(QDT a, const QD_ASQ qd_real &b);
QD_API bool operator>(const QD_ASQ qd_real &a, QDT b);

QD_API bool operator<=(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API bool operator<=(const qd_real &a, const dd_real &b);
QD_API bool operator<=(const dd_real &a, const qd_real &b);
#endif
QD_API bool operator<=(QDT a, const QD_ASQ qd_real &b);
QD_API bool operator<=(const QD_ASQ qd_real &a, QDT b);

QD_API bool operator>=(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API bool operator>=(const qd_real &a, const dd_real &b);
QD_API bool operator>=(const dd_real &a, const qd_real &b);
#endif
QD_API bool operator>=(QDT a, const QD_ASQ qd_real &b);
QD_API bool operator>=(const QD_ASQ qd_real &a, QDT b);

QD_API bool operator!=(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
#ifdef _QD_DD_REAL_H
QD_API bool operator!=(const qd_real &a, const dd_real &b);
QD_API bool operator!=(const dd_real &a, const qd_real &b);
#endif
QD_API bool operator!=(QDT a, const QD_ASQ qd_real &b);
QD_API bool operator!=(const QD_ASQ qd_real &a, QDT b);

QD_API qd_real fabs(const QD_ASQ qd_real &a);
QD_API qd_real abs(const QD_ASQ qd_real &a);    /* same as fabs */

QD_API qd_real ldexp(const QD_ASQ qd_real &a, int n);

QD_API qd_real nint(const QD_ASQ qd_real &a);
QD_API qd_real quick_nint(const QD_ASQ qd_real &a);
QD_API qd_real floor(const QD_ASQ qd_real &a);
QD_API qd_real ceil(const QD_ASQ qd_real &a);
QD_API qd_real aint(const QD_ASQ qd_real &a);

QD_API qd_real sin(const QD_ASQ qd_real &a);
QD_API qd_real cos(const QD_ASQ qd_real &a);
QD_API qd_real tan(const QD_ASQ qd_real &a);
QD_API void sincos(const QD_ASQ qd_real &a, QD_ASQ qd_real &s, QD_ASQ qd_real &c);

QD_API qd_real asin(const QD_ASQ qd_real &a);
QD_API qd_real acos(const QD_ASQ qd_real &a);
QD_API qd_real atan(const QD_ASQ qd_real &a);
QD_API qd_real atan2(const QD_ASQ qd_real &y, const QD_ASQ qd_real &x);

QD_API qd_real exp(const QD_ASQ qd_real &a);
QD_API qd_real log(const QD_ASQ qd_real &a);
QD_API qd_real log10(const QD_ASQ QD_ASQ qd_real &a);

QD_API qd_real sinh(const QD_ASQ qd_real &a);
QD_API qd_real cosh(const QD_ASQ qd_real &a);
QD_API qd_real tanh(const QD_ASQ qd_real &a);
QD_API void sincosh(const QD_ASQ qd_real &a, QD_ASQ qd_real &sin_qd, QD_ASQ qd_real &cos_qd);

QD_API qd_real asinh(const QD_ASQ qd_real &a);
QD_API qd_real acosh(const QD_ASQ qd_real &a);
QD_API qd_real atanh(const QD_ASQ qd_real &a);

QD_API qd_real qdrand(void);

QD_API qd_real max(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
QD_API qd_real max(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b, const QD_ASQ qd_real &c);
QD_API qd_real min(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);
QD_API qd_real min(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b, const QD_ASQ qd_real &c);

QD_API qd_real fmod(const QD_ASQ qd_real &a, const QD_ASQ qd_real &b);

#if !QD_GPU
QD_API std::ostream &operator<<(std::ostream &s, const qd_real &a);
QD_API std::istream &operator>>(std::istream &s, qd_real &a);
#endif
#ifdef QD_INLINE
#if !QD_GPU
#include <qd/qd_inline.h>
#else
#include "qd_inline.h"
#endif
#endif

#endif /* _QD_QD_REAL_H */

