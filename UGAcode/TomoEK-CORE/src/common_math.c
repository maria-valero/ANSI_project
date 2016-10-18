/*
 * common_math.cpp
 *
 *  Created on: Feb 22, 2015
 *      Author: nishita
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include "../include/common_math.h"

/* Used for accessing the integer representation of double precision
 * floating-point numbers. */
union Double_t {
	int64_t i;
	double f;
};

	bool doubleAlmostEqualUlps(double A, double B, int maxUlpsDiff) {
		/* Adapted from AlmostEqualUlps,
		 * http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
		 *
		 * A, B 	  : two floating-point numbers (double precision) to compare.
		 * maxUlpsDiff: maximum spacing between the floating-point numbers A and B.
		 * ULP		  = unit in the last place or unit of least precision.
		 */
		union Double_t uA, uB;
		bool signedA, signedB;
		int64_t ulpsDiff;

		/* Ensure that either A or B are not close to zero. */
		assert ( (fabs(A) > 5 * DBL_EPSILON) || (fabs(B) > 5 * DBL_EPSILON) );

		/* Initialize unions with floats. */
		uA.f = A;
		uB.f = B;

		/* Extract sign bits. */
		signedA = (uA.i >> 63) != 0;
		signedB = (uB.i >> 63) != 0;

		/* Different signs means they do not match. */
		if (signedA != signedB) {
			/* Check for equality to make sure +0==-0 */
			if (A == B)
				/* Might be reached if assert is deactivated (-DNDEBUG) */
				return true;
			return false;
		}

		/* Find the difference in ULPs */
		ulpsDiff = int64_abs(uA.i - uB.i);
		if (ulpsDiff <= maxUlpsDiff)
			return true;

		return false;
		}
	bool doubleAlmostEqualUlpsAndAbs(double A, double B,double maxDiff, int maxUlpsDiff) {
		/* Adapted from AlmostEqualUlpsAndAbs,
		 * http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
		 *
		 * A, B   : two floating-point numbers (double precision) to compare.
		 * maxDiff: epsilon for floating-point absolute epsilon check (should be some
		 *          small multiple of DBL_EPSILON.
		 * maxUlps: maximum spacing between the floating-point numbers A and B.
		 * ULP    = unit in the last place or unit of least precision.
		 */
		double absDiff;
		union Double_t uA, uB;
		bool signedA, signedB;
		int64_t ulpsDiff;

		/* Check if the numbers are really close -- needed when comparing numbers near zero. */
		absDiff = fabs(A - B);
		if (absDiff <= maxDiff)
			return true;

		/* Initialize unions with floats. */
		uA.f = A;
		uB.f = B;

		/* Extract sign bits. */
		signedA = (uA.i >> 63) != 0;
		signedB = (uB.i >> 63) != 0;

		/* Different signs means they do not match. */
		if (signedA != signedB)
			return false;

		/* Find the difference in ULPs */
		ulpsDiff = int64_abs(uA.i - uB.i);
		if (ulpsDiff <= maxUlpsDiff)
			return true;

		return false;
	}


