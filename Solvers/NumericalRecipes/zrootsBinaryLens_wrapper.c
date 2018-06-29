//zrootsBinaryLens_wrapper.c

#define FLOAT double
#define float double
#include <stdlib.h>
#include <stdio.h>
#include "complex.h"

FLOAT* zroots_5(FLOAT p0, FLOAT p1, 
		FLOAT p2, FLOAT p3, FLOAT p4, FLOAT p5, FLOAT p6, 
		FLOAT p7, FLOAT p8, FLOAT p9, FLOAT p10, FLOAT p11) {

	void zroots(fcomplex a[], int m, fcomplex roots[], int polish);
	fcomplex a[6], roots[5];
	static FLOAT sep_roots[10];
	int polish = 0;
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

FLOAT* zroots_10(FLOAT p0, FLOAT p1, FLOAT p2, FLOAT p3, FLOAT p4, FLOAT p5,
		FLOAT p6, FLOAT p7, FLOAT p8, FLOAT p9, FLOAT p10, FLOAT p11,
		FLOAT p12, FLOAT p13, FLOAT p14, FLOAT p15, FLOAT p16, FLOAT p17,
		FLOAT p18, FLOAT p19, FLOAT p20, FLOAT p21) {

	void zroots(fcomplex a[], int m, fcomplex roots[], int polish);
	fcomplex a[11], roots[10];
	static FLOAT sep_roots[20];
	int polish = 0;
	int i;

	a[0] = Complex(p0, p11);
	a[1] = Complex(p1, p12);
	a[2] = Complex(p2, p13);
	a[3] = Complex(p3, p14);
	a[4] = Complex(p4, p15);
	a[5] = Complex(p5, p16);
	a[6] = Complex(p6, p17);
	a[7] = Complex(p7, p18);
	a[8] = Complex(p8, p19);
	a[9] = Complex(p9, p20);
	a[10] = Complex(p10, p21);

	zroots(a, 10, roots, polish);

	for (i=0; i<10; i++) {
		sep_roots[i] = roots[i+1].r;
		sep_roots[i+10] = roots[i+1].i;
	}

	return sep_roots;
}

