//zrootsBinaryLens_wrapper.c

#include <stdlib.h>
#include <stdio.h>
#include "complex.h"
#define FLOAT double

FLOAT* zroots_5(FLOAT p0, FLOAT p1, 
		FLOAT p2, FLOAT p3, FLOAT p4, FLOAT p5, FLOAT p6, 
		FLOAT p7, FLOAT p8, FLOAT p9, FLOAT p10, FLOAT p11) {

	void zroots(fcomplex a[], int m, fcomplex roots[], int polish);
	fcomplex a[6], roots[5];
	static FLOAT sep_roots[10];
	int polish = 1;
	int i;

	a[0] = Complex(p0, p6);
	a[1] = Complex(p1, p7);
	a[2] = Complex(p2, p8);
	a[3] = Complex(p3, p9);
	a[4] = Complex(p4, p10);
	a[5] = Complex(p5, p11);

	zroots(a, 5, roots, polish);

	for (i=0; i<5; i++) {
		sep_roots[i] = roots[i+1].r;
		sep_roots[i+5] = roots[i+1].i;
	}

	return sep_roots;
}