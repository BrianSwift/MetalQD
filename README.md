# Metal QD: QD extended precision library for Metal

Metal QD is _(preliminary work on)_ a port of the QD extended precision library to Apple's Metal GPU compute language.

QD is a library for Quad-Double arithmetic in which extended precision numbers are represented by the unevaluated sum of four double precision machine numbers. QD comes from http://crd-legacy.lbl.gov/~dhbailey/mpdist/

Metal QD adapts QD to the Metal environment primarily by using single precision floats instead of doubles. This allows calculations with 96-bits of precision. This port also disables references to features unavailable on the GPU, such as `iostream, string, limits, cmath`.

Metal QD is based off of `qd-2.3.17`

What is supported?
* Constructors
* Arithmetic operators: `+, -, *, =, +=, -=`
* Comparison operators: `<, >, ==, !=`
* Special Functions: `sqr(), mul_pwr2()`

What hasn't been ported?
* Anything not implemented in the _inline_ library.
* Division
* Trig and other transendental functions
* The DD (Double-Double) type
* I/O

## Usage

* Copy these files to project: `metal_qd.h, qd_real.h, qd_config.h, qd_inline.h, inline.h`

* In Metal or Objective-C++ source file: `#include "metal_qd.h"`

* Disable Metal's unsafe floating-point optimizations.
Add **-fno-fast-math** to Build Phases -> Compile Sources -> Compiler Flags for any Metal files using `metal_qd.h`
* Read the `README_QD` section **C. Programming techniques**


## Example Code Snipet

    qd_real q;
    q=1.;            // native 1.0 converted to qd_real q
    q+=exp2(-95.);   // native 2^-95 added to qd_real q
    // q now represents 1.000000000000000000000000000025
    // This has more precision than can be represented in a native double
    
    // Conversion to/from float4 usefull for passing values between CPU and GPU
    simd::float4 vq;
    vq=to_float4(q); // Convert qd_real to float4
    
    qd_real qvq;     // Convert and float4 back to qd_real
    qvq=vq;

Find the above in `tests/metal_qd_test.cpp`

## Future Work

* Port non-inline and `DD` portions of library
* Create proper Metal framework
* Create demo app using Metal QD on GPU (in addition to CPU)
* Test having both QD and Metal_QD in same CPU source
* Update documentation to reflect Metal_QD
