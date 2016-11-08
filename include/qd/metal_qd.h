//
//  metal_qd.h
//  Mandelbits
//
//  Created by Brian Swift on 10/8/16.
//  Copyright © 2016 Brian Swift. All rights reserved.
//

//  Various definitions to configure QD library to be quad float based,
//  and to be usable in both Metal on GPU and C++ on CPU.

#ifndef metal_qd_h
#define metal_qd_h

#define QDT float
#define QDT_float 1
#define QDT_double 0
#define QD_GPU 1        /* Disables I/O and other definitions not ported to GPU yet */

#ifdef __METAL_VERSION__
#define QD_ASQ thread   /* address space qualifier */
#else
#define QD_ASQ /* */
#endif

#ifdef __METAL_VERSION__
namespace std = metal;       /* fixes various uses of std::math_function */
#else
#import <cmath>
#endif

#ifdef __METAL_VERSION__
#define QD_FMS(x,y,z) precise::fma(x,y,-z)
#else
#define QD_FMS(x,y,z) fma(x,y,-z)
#endif

// #define inline /* */ /* to disable inline for dubugging and analysis */

// enable float4 conversions if SIMD available
#if defined(SIMD_LIBRARY_VERSION)||defined(__METAL_VERSION__)||defined(matrix_add)
    #define QD_HAVE_SIMD 1
#else
    #define QD_HAVE_SIMD 0
#endif


#include "qd_real.h"


#ifndef __METAL_VERSION__
inline qd_real::qd_real(double d) {
    x[0]=d;
    x[1]=d-x[0];
    x[2]=d-x[0]-x[1];
    x[3]=d-x[0]-x[1]-x[2]; // this should be 0
    qd::renorm(x[0],x[1],x[2],x[3]);
}
#endif


#if QD_HAVE_SIMD
// define conversions from and to float4
inline qd_real::qd_real(simd::float4 xx) {
    x[0] = xx.x;
    x[1] = xx.y;
    x[2] = xx.z;
    x[3] = xx.w;
}

inline simd::float4 to_float4(const QD_ASQ qd_real &a) {
    simd::float4 t={a[0],a[1],a[2],a[3]};
    return t;
}
#endif


#endif /* metal_qd_h */
