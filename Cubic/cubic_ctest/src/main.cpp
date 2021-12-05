// grid_test.cpp : Defines the entry point for the application.
//


#include "cubic/cubic.h"

#include<array>
#include<stdexcept>
#include<assert.h>

template<typename FP>
void assert_zero(FP val) {
	constexpr FP EPSILON = std::is_same<int, typename std::remove_cv<FP>::type>::value ?
		0 :
		std::is_same<float, typename std::remove_cv<FP>::type>::value ? (FP)FLT_EPSILON : (FP)DBL_EPSILON;
	if (std::abs(val) > EPSILON) {
		throw std::runtime_error("Not zero.");
	}
}

/* Test cases for 'cubic_roots()' for finding roots evaluating to zero in double precision.
*/
static void test_cubic_zero_double() {
	double out[3];
	double val;
	int N;

	/* Double roots: 0, 0, 1*/
	N = cubic_roots<double>(1.0, -1.0, 0.0, 0.0, out);
	assert_zero(N - 3);
	val = cubic(1.0, -1.0, 0.0, 0.0, out[0]);
	assert_zero(val);
	val = cubic(1.0, -1.0, 0.0, 0.0, out[1]);
	assert_zero(val);

}

/* Test cases for 'quadratic_roots()' for finding roots evaluating to zero in double precision.
*/
static void test_quadratic_zero_double() {
	double out[2];
	double val;
	int N;

	/* Linear*/
	N = quadratic_roots<double>(0.0, 1.0, 3.0, out);
	assert_zero(N - 1);
	assert_zero(quadratic(0.0, 1.0, 3.0, out[0]));

	/* Roots: 0, 1*/
	N = quadratic_roots<double>(1.0, -1.0, 0.0, out);
	assert_zero(N - 2);
	assert_zero(quadratic(1.0, -1.0, 0.0, out[0]));
	assert_zero(quadratic(1.0, -1.0, 0.0, out[1]));

}


static void testcases() {
	double out[3];
	double val;
	int N;

	N = cubic_roots<double>(9933.79603159, 7658.09556667, 2059.57642168, 187.33806459, out);
	val = cubic(9933.79603159, 7658.09556667, 2059.57642168, 187.33806459, out[0]);

	N = cubic_roots<double>(2.04424068e-02, -7.05085740e+04, 2.59523337e+03, 3.18743329e+04, out);
	val = cubic(2.04424068e-02, -7.05085740e+04, 2.59523337e+03, 3.18743329e+04, out[0]);

	N = cubic_roots<double>(2.14650209e-01, 9.66306444e+04, 8.98525869e+04, 5.26231876e+04, out);
	val = cubic(2.14650209e-01, 9.66306444e+04, 8.98525869e+04, 5.26231876e+04, out[0]);
}


int main(void) {

	test_quadratic_zero_double();
	test_cubic_zero_double();

	testcases();
	return 0;
}
