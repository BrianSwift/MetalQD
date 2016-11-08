/*
 * include/inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file contains the basic functions used both by double-double
 * and quad-double package.  These are declared as inline functions as
 * they are the smallest building blocks of the double-double and 
 * quad-double arithmetic.
 */
#ifndef _QD_INLINE_H
#define _QD_INLINE_H

// Configure for non-metal-ported behavior
#ifndef QDT
#define QDT double
#define QDT_double 1
#define QDT_float 0
#define QD_ASQ /* */
#define QD_GPU 0
#define QD_HAVE_SIMD 0
#endif


#if QDT_double
#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996  (996=Emax"1023"-27)
#elif QDT_float
#define _QD_SPLITTER 4097.0                 // = 2^12 + 1
#define _QD_SPLIT_THRESH 4.153837487E34     // = 2^115 (115=Emax"127"-12)
#else
#error "Must define QDT_double or QDT_float"
#endif

#ifdef QD_VACPP_BUILTINS_H
/* For VisualAge C++ __fmadd */
#include <builtins.h>
#endif

#if !QD_GPU
#include <cmath>
#include <limits>
#else
#endif

namespace qd {

#if !QD_GPU
static const QDT _d_nan = std::numeric_limits<QDT>::quiet_NaN();
static const QDT _d_inf = std::numeric_limits<QDT>::infinity();
#else
#endif
    
#ifdef __FAST_MATH__
#error Following code incompatible with -ffast-math , needs -fno-fast-math
#endif

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline QDT quick_two_sum(QDT a, QDT b, QD_ASQ QDT &err) {
  QDT s = a + b;
  err = b - (s - a);
  return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline QDT quick_two_diff(QDT a, QDT b, QD_ASQ QDT &err) {
  QDT s = a - b;
  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
inline QDT two_sum(QDT a, QDT b, QD_ASQ QDT &err) {
  QDT s = a + b;
  QDT bb = s - a;
  err = (a - (s - bb)) + (b - bb);
  return s;
}

/* Computes fl(a-b) and err(a-b).  */
inline QDT two_diff(QDT a, QDT b, QD_ASQ QDT &err) {
  QDT s = a - b;
  QDT bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

#if QD_GPU
#ifndef QD_FMS
#error QD_FMS should be defined when using metal
#endif
#endif

#ifndef QD_FMS
/* Computes high word and lo word of a */
inline void split(QDT a, QDT &hi, QDT &lo) {
  QDT temp;
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
#if QDT_double
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
    hi *= 268435456.0;          // 2^28
    lo *= 268435456.0;          // 2^28
#elif QDT_float
      a *= .00012207031250000000;  // 2^-13
      temp = _QD_SPLITTER * a;
      hi = temp - (temp - a);
      lo = a - hi;
      hi *= 8192.;          // 2^13
      lo *= 8192.;          // 2^13
#else
#error "Must define QDT_double or QDT_float"
#endif
  } else {
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
  }
}
#endif

/* Computes fl(a*b) and err(a*b). */
inline QDT two_prod(QDT a, QDT b, QD_ASQ QDT &err) {
#ifdef QD_FMS
  QDT p = a * b;
  err = QD_FMS(a, b, p);
  return p;
#else
  QDT a_hi, a_lo, b_hi, b_lo;
  QDT p = a * b;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  return p;
#endif
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline QDT two_sqr(QDT a, QD_ASQ QDT &err) {
#ifdef QD_FMS
  QDT p = a * a;
  err = QD_FMS(a, a, p);
  return p;
#else
  QDT hi, lo;
  QDT q = a * a;
  split(a, hi, lo);
  err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
  return q;
#endif
}

/* Computes the nearest integer to d. */
inline QDT nint(QDT d) {
  if (d == std::floor(d))
    return d;
  return std::floor(d + 0.5);
}

/* Computes the truncated integer. */
inline QDT aint(QDT d) {
  return (d >= 0.0) ? std::floor(d) : std::ceil(d);
}

/* These are provided to give consistent
   interface for double with double-double and quad-double. */
inline void sincosh(QDT t, QD_ASQ QDT &sinh_t, QD_ASQ QDT &cosh_t) {
  sinh_t = std::sinh(t);
  cosh_t = std::cosh(t);
}

inline QDT sqr(QDT t) {
  return t * t;
}

inline QDT to_double(QDT a) { return a; }
inline int    to_int(QDT a) { return static_cast<int>(a); }

}

#endif /* _QD_INLINE_H */
