/*
 * common_math.h
 *
 *  Created on: Feb 16, 2015
 *      Author: nishita
 */

#ifndef COMMON_MATH_H_
#define COMMON_MATH_H_

#include <stdbool.h>
#include <stdint.h>
#include <float.h>
//#include "../include/gmt_notposix.h"
#include "../include/s_rint.h"
	/* Limit casting to one place (here) for dropping lrint output to signed or unsigned ints */
#define irint(x) ((int)lrint(x))
#define urint(x) ((unsigned int)lrint(x))
#define irintf(x) ((int)lrintf(x))
#define urintf(x) ((unsigned int)lrintf(x))
#define int64_abs(x) ((int64_t)(((x) >= 0) ? (x) : -(x)))

#define doubleAlmostEqual(A, B) (doubleAlmostEqualUlps(A, B, 5))
#define doubleAlmostEqualZero(A, B) (doubleAlmostEqualUlpsAndAbs(A, B, 5*DBL_EPSILON, 5))

#define GMT_360_RANGE(w,e) (doubleAlmostEqual (fabs((e) - (w)), 360.0))
#define GMT_180_RANGE(s,n) (doubleAlmostEqual (fabs((n) - (s)), 180.0))
#define GMT_IS_POLE(y) (doubleAlmostEqual (fabs(y), 90.0))

bool doubleAlmostEqualUlps(double A, double B, int maxUlpsDiff);
bool doubleAlmostEqualUlpsAndAbs(double A, double B,double maxDiff, int maxUlpsDiff);
#endif /* COMMON_MATH_H_ */
