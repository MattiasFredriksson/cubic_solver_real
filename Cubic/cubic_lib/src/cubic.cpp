/*
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
*/
#include "cubic/cubic.h"
#include <math.h>
#include <cmath>
#include <float.h>


template<typename FP>
FP quadratic(FP a, FP b, FP c, FP x)
{
	return a * x * x + b * x + c;
}
template double quadratic(double a, double b, double c, double x);
template float quadratic(float a, float b, float c, float x);

template<typename FP>
FP cubic(FP a, FP b, FP c, FP d, FP x)
{
	FP xsq = x * x;
	return a * x * xsq + b * xsq + c * x + d;
}
template double cubic(double a, double b, double c, double d, double x);
template float cubic(float a, float b, float c, float d, float x);

/**
* Find the roots to the quadratic equation
*	f(x) = ax^2 + bx + c
*
*/
template<typename FP>
int quadratic_roots(FP a, FP b, FP c, FP* xout) {
	using namespace std;
	constexpr FP half = (FP)0.5;
	constexpr FP EPSILON = is_same<float, typename remove_cv<FP>::type>::value ? FLT_EPSILON : DBL_EPSILON;

	if (abs(a) < EPSILON)
	{
		/* Linear equation */
		if (abs(b) > EPSILON)
		{
			*xout = -c / b;
			return 1;
		}
	}
	else
	{
		/* Quadratic equation

		Implementation is based on https://people.csail.mit.edu/bkph/articles/Quadratics.pdf.
		The root combination is swapped however as the one in the paper performed far worse on MAE tests.
		*/

		/* Reduce form through division, multiplication of '1.0 / a' has a significant precision cost. */
		/*a = (FP)1.0;*/
		b = b / a;
		c = c / a;

		FP yy = b * b - (FP)4.0 * c;
		if (yy >= (FP)0.0)
		{
			FP c2 = (FP)2.0 * c;
			FP y = sqrt(yy);
			if (b < (FP)0.0) {
				xout[0] = c2 / (y - b);
				xout[1] = (y - b) * half;
			}
			else {
				xout[0] = (-b - y) * half;
				xout[1] = c2 / (-y - b);
			}
			return 2;
		}
	}
	return 0;
}
template int quadratic_roots(double a, double b, double c, double* xout);
template int quadratic_roots(float a, float b, float c, float* xout);

/**
 * Implementation of Cardano's method for solving cubic equations.
 *
 * Implementation is based on https://github.com/tatwood/solvecubic
 * @author	  Thomas Atwood (original), Mattias Fredriksson (modified)
 * @date      2011 (cloned Nov 2021)
 * @copyright unlicense / public domain
 ****************************************************************************/
template<typename FP>
int cubic_roots(FP a, FP b, FP c, FP d, FP* xout)
{
	using namespace std;
	constexpr FP cos120 = (FP)-0.5;
	constexpr FP sin120 = (FP)0.8660254037844386467637231707529361834714026269051903140279034897;
	constexpr FP PI = (FP)3.141592653589793238462643383279502884197169399375105820974944592307816406286;
	constexpr FP PIHalf = PI / 2.0;
	constexpr FP PI2 = PI * 2.0;
	constexpr FP PI2over3 = PI * 2.0 / 3.0;
	constexpr FP third = (FP)(1.0 / 3.0);
	constexpr FP zero = (FP)0.0;
	constexpr FP EPSILON = is_same<float, typename remove_cv<FP>::type>::value ? FLT_EPSILON : DBL_EPSILON;

	int n = 0;
	if (abs(d) < EPSILON)
	{
		/* First solution is x = 0 */
		*xout = zero;
		n = 1;
		++xout;
		/* Divide all terms by x, converting to quadratic equation */
		d = c;
		c = b;
		b = a;
		a = zero;
	}
	if (abs(a) < EPSILON)
	{
		return quadratic_roots<FP>(b, c, d, xout) + n;
	}
	else

		/* Cubic equation */ {
		/* Reduce form through division, multiplication of '1.0 / a' has a (small) precision cost. */
		b = b / a;
		c = c / a;
		d = d / a;
		a = (FP)1.0;

		FP bover3 = b * third;
		FP p = c - bover3 * b;
		FP halfq = bover3 * bover3 * bover3 - (FP)0.5 * bover3 * c + (FP)0.5 * d;
		FP yy = p / (FP)27.0 * p * p + halfq * halfq;

		if (yy < (FP)0.0) /* Sqrt is negative: three real solutions */
		{
			n = 3;
			if (fabs(p) < EPSILON)
			{
				xout[0] = -bover3;
				xout[1] = xout[0];
				xout[2] = xout[0];
			}
			else
			{
				FP uu = (FP)(-4.0 / 3.0) * p;
				FP u = sqrt(uu);
				FP theta = acos((FP)-8.0 * halfq / (u * uu)) * third;
				xout[0] = u * cos(theta) - bover3;
				xout[1] = u * cos(theta - PI2over3) - bover3;
				xout[2] = u * cos(theta + PI2over3) - bover3;
			}
		}
		else
		{
			/*  Sqrt is positive: one real solution */
			FP y = sqrt(yy);
			FP uuu = y - halfq;
			FP vvv = -y - halfq;
			FP www = abs(uuu) > abs(vvv) ? uuu : vvv;
			FP w = copysign(cbrt(abs(www)), www);
			*xout = w - p / ((FP)3.0 * w) - bover3;
			n = 1;
			return n;
		}
		return n;
	}

	return n;
}
template int cubic_roots(double a, double b, double c, double d, double* xout);
template int cubic_roots(float a, float b, float c, float d, float* xout);
