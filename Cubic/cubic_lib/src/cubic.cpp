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
* Implementation is based on https://people.csail.mit.edu/bkph/articles/Quadratics.pdf.
*/
template<typename FP>
int quadratic_roots(FP a, FP b, FP c, FP* xout) {
	using namespace std;
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

		Based on https://people.csail.mit.edu/bkph/articles/Quadratics.pdf.
		The combination is swapped as the combination described in the paper performed far worse on MAE tests.
		*/

		/* Reduce form through division, multiplication of '1.0 / a' has a significant precision cost. */
		/*a = (FP)1.0;*/
		b = b / a;
		c = c / a;

		FP q = b * b - (FP)4.0 * c;
		if (q >= (FP)0.0)
		{
			FP c2 = (FP)2.0 * c;
			q = sqrt(q);
			if (b < (FP)0.0) {
				xout[0] = c2 / (q - b);
				xout[1] = (q - b) * (FP)0.5;
			}
			else {
				xout[0] = (-b - q) * (FP)0.5;
				xout[1] = c2 / (-q - b);
			}
			return 2;
		}
	}
	return 0;
}
template int quadratic_roots(double a, double b, double c, double* xout);
template int quadratic_roots(float a, float b, float c, float* xout);

/**
 * Implementation uses both the trignometric and Cardano's method method for solving cubic equations.
 *
 * Implementation is based on https://github.com/tatwood/solvecubic
 * @author	  Thomas Atwood (original author), Mattias Fredriksson (optimized)
 * @date      2011 (cloned Nov 2021)
 * @copyright unlicense / public domain
 ****************************************************************************/
template<typename FP>
int cubic_roots(FP a, FP b, FP c, FP d, FP* xout)
{
	using namespace std;
	constexpr FP PI = (FP)3.141592653589793238462643383279502884197169399375105820974944592307816406286;
	constexpr FP PIHalf = (FP)(PI / 2.0);
	constexpr FP PI2 = (FP)(PI * 2.0);
	constexpr FP PI2over3 = (FP)(PI * 2.0 / 3.0);
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
		//a = (FP)1.0;

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


/**
* Find the roots to the quadratic equation
*	f(x) = ax^2 + bx + c
*
* Implementation is based on https://people.eecs.berkeley.edu/~wkahan/Math128/Cubic.pdf.
* 'To Solve a Real Cubic Equation' authored by W. Kahan.
* 
* Note* implementation only return real roots and checks if the equation is linear.
*/
template<typename FP>
inline int qdrtc(FP A, FP B, FP C, FP* xout)
{
	constexpr FP EPSILON = std::is_same<float, typename std::remove_cv<FP>::type>::value ? FLT_EPSILON : DBL_EPSILON;

	if (std::abs(A) < EPSILON)
	{
		/* Linear equation */
		if (std::abs(B) > EPSILON)
		{
			*xout = -C / B;
			return 1;
		}
		/* Constant*/
		return 0;
	}

	FP b = -B / (FP)2.0;
	FP q = b * b - A * C;
	if (q < (FP)0.0) {
		return 0;
		/* Complex roots
		X1 = b / A;
		X2 = X1;
		Y1 = std::sqrt(-q) / A;
		Y2 = -Y1;
		*/
	}
	else {
		FP r = b + std::copysign(std::sqrt(q), b); /* sqrt(q) * sign(b) as q >= 0 */
		if (r == (FP)0.0) {

			xout[0] = C / A;
			xout[1] = -xout[0];
		}
		else {
			xout[0] = C / r;
			xout[1] = r / A;
		}
	}
	return 2;
}


template<typename FP>
inline void qbc_eval(FP X, FP A, FP B, FP C, FP D, FP& Q, FP& Q_p, FP& B1, FP& C2)
{
	FP q0 = A * X;
	B1 = q0 + B;
	C2 = B1 * X + C;
	Q_p = (q0 + B1) * X + C2;
	Q = C2 * X + D;
}

/**
 * Compute the real roots for the cubic equation
 *
 *		ax^3 + bx^2 + cx + d = 0
 *
 * Implementation is based on https://people.eecs.berkeley.edu/~wkahan/Math128/Cubic.pdf.
 * 'To Solve a Real Cubic Equation' authored by W. Kahan.
 * 
 * Note* implementation only return real roots and checks if the equation is linear.
 */
template<typename FP>
int cubic_roots_qbc(FP A, FP B, FP C, FP D, FP* xout) {
	using namespace std;
	constexpr FP EPSILON = is_same<float, typename remove_cv<FP>::type>::value ? FLT_EPSILON : DBL_EPSILON;

	int N = 0;
	FP b1, c2;
	if (abs(A) < EPSILON) {
		/* Quadratic equation */
		A = B;
		b1 = C;
		c2 = D;
		// *xout++ == INFINITY;
	}
	else if (abs(D) < EPSILON) {
		/* Convert to a quadratic equation (divide by x) */
		*xout++ = (FP)0.0;
		b1 = B;
		c2 = C;
		N = 1;
	}
	else {
		FP X = -(B / A) / (FP)3.0;
		FP q, q_p, t, r, s;
		qbc_eval(X, A, B, C, D, q, q_p, b1, c2);

		t = q / A;
		r = std::cbrt(std::abs(t));
		s = std::copysign((FP)1.0, t);

		t = -q_p / A;
		if (t > 0) {
			r = (FP)1.324717957244746025960908854478097340734404056901733365 * std::fmax(r, std::sqrt(t));
		}

		FP x0 = X - r * s;
		if (x0 != X) {
			int i = 0;
			do {
				X = x0;
				qbc_eval(X, A, B, C, D, q, q_p, b1, c2);
				if (q_p == 0) {
					x0 = X;
				}
				else {
					x0 = X - (q / q_p) / (FP)1.000000000000001; /* 1.000..001 */
				}
			} while (x0 * s > X * s);

			if (std::abs(A) * X * X > std::abs(D / X)) {
				c2 = -D / X;
				b1 = (c2 - C) / X;
			}
		}
		N = 1;
		*xout++ = X;
	}
	return N + qdrtc(A, b1, c2, xout);
}
template int cubic_roots_qbc(double a, double b, double c, double d, double* xout);
template int cubic_roots_qbc(float a, float b, float c, float d, float* xout);