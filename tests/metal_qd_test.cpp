/*
  Basic test of Metal_QD on CPU

  Build and run:
  c++ metal_qd_test.cpp -I../include/qd -o metal_qd_test && ./metal_qd_test

  Expected Output:
  qvq = 1.000000000000000000000000000025
  vq = 0x1p+0 0x1p-95 0x0p+0 0x0p+0
  d = 1 0x1p+0
*/

#import <simd/simd.h>

#include "metal_qd.h"

#import <iostream>

using namespace std;

#if QDT_double
#define print_precision 67
#elif QDT_float
#define print_precision 30
#else
#error "Must define QDT_double or QDT_float"
#endif

// A hack to print simple QD numbers.
// Actual code to correctly print QD value is significantly more complex.
// See to_string and to_digits in qd_real.cpp
void print_qd(const qd_real &a){
  qd_real x = a;
	int t = trunc(x[0]);
	cout << t ;
	cout << "." ;
	x -= t;
	for(int i = 0 ; i < print_precision ; i++){
	  x *= 10. ;
	  t = trunc(x[0]);
	  cout << t ;
	  x -= t;
	}
}

int main() {
	
  qd_real q;
  q=1.;            // native 1.0 converted to qd_real q
  q+=exp2(-95.);   // native 2^-95 added to qd_real q
                   // q now represents 1.000000000000000000000000000025
                   // This has more precision than can be represented in a native double

  simd::float4 vq;
  vq=to_float4(q); // Convert qd_real to float4

  qd_real qvq;     // Convert and float4 back to qd_real
  qvq=vq;


  // Demonstrate machine double lacks sufficient precision to represent this number
  double d;
  d=1.0;
  d+=exp2(-95.);   // d remains 1.0
  
  cout << "qvq = " ; print_qd(qvq);	cout << endl;  // qvq = 1.000000000000000000000000000025

  printf("vq = %a %a %a %a\n",vq.x,vq.y,vq.z,vq.w);    // vq = 0x1p+0 0x1p-95 0x0p+0 0x0p+0

  printf("d = %g %a\n",d,d);   // d = 1 0x1p+0

  return 0;
}
