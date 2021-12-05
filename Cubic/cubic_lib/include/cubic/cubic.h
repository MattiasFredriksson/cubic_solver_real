#pragma once
/* Algorithms for computing real roots of a cubic or quadratic equation (3rd or 2nd order polynomial).
*
* Provided under the Unlicense License (public domain, see license terms at http://unlicense.org/).
*/

/**
* Evaluate the quadratic function for a given x.
*
 * The quadratic function is a polynomial function on the form
*	f(x) = ax^2 + bx + c
*/
template<typename FP>
FP quadratic(FP a, FP b, FP c, FP x);

/**
* Evaluate the cubic function for a given x.
*
 * The cubic function is a polynomial function on the form
*	f(x) = ax^3 + bx^2 + cx + d
*/
template<typename FP>
FP cubic(FP a, FP b, FP c, FP d, FP x);

/**
* Compute the real roots for the quadratic equation
*
*		ax^2 + bx + c = 0
*/
template<typename FP>
int quadratic_roots(FP a, FP b, FP c, FP* xout);

/**
 * Compute the real roots for the cubic equation
 *
 *		ax^3 + bx^2 + cx + d = 0
 */
template<typename FP>
int cubic_roots(FP a, FP b, FP c, FP d, FP* xout);