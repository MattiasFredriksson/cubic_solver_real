// grid_test.cpp : Defines the entry point for the application.
//


#include "cubic/cubic.h"

#include<array>
#include<stdexcept>
#include<assert.h>
#include <random>
#include <chrono>
#include<iostream>

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
template<typename FP>
static void test_cubic_zero_double(CBRT_SOLVER<FP> cbrt_solver) {
	FP out[3];
	FP val;
	int N;

	/* Double roots: 0, 0, 1*/
	N = cbrt_solver(1.0, -1.0, 0.0, 0.0, out);
	assert_zero(N - 3);
	val = cubic<FP>(1.0, -1.0, 0.0, 0.0, out[0]);
	assert_zero(val);
	val = cubic<FP>(1.0, -1.0, 0.0, 0.0, out[1]);
	assert_zero(val);

	N = cbrt_solver(6111, -51792, 109737, 0.00623, out);
	//assert_zero(N - 3);
	val = cubic<FP>(6111, -51792, 109737, 0.00623, out[0]);
	val = cubic<FP>(6111, -51792, 109737, 0.00623, out[1]);
	val = cubic<FP>(6111, -51792, 109737, 0.00623, out[2]);

	N = cbrt_solver(6111, 25939.92, 124.1808, 0.148694, out);
	//assert_zero(N - 3);
	val = cubic<FP>(6111, 25939.92, 124.1808, 0.148694, out[0]);
	val = cubic<FP>(6111, 25939.92, 124.1808, 0.148694, out[1]);
	val = cubic<FP>(6111, 25939.92, 124.1808, 0.148694, out[2]);

	N = cbrt_solver(658, -190125, 18311811, -587898164, out);
	// assert_zero(N - 1);
	val = cubic<FP>(658, -190125, 18311811, -587898164, out[0]);
	//assert_zero(val);

	N = cbrt_solver(658, 190125, 18311811, -587898164, out);
	//assert_zero(N - 1);
	val = cubic<FP>(658, 190125, 18311811, -587898164, out[0]);
	//assert_zero(val);
}

/* Test cases for 'quadratic_roots()' for finding roots evaluating to zero in double precision.
*/
template<typename FP>
static void test_quadratic_zero_double(QDRT_SOLVER<FP> qdrt_solver) {
	FP out[2];
	int N;

	/* Linear*/
	N = qdrt_solver(0.0, 1.0, 3.0, out);
	assert_zero(N - 1);
	assert_zero(quadratic(0.0, 1.0, 3.0, out[0]));

	/* Roots: 0, 1*/
	N = qdrt_solver(1.0, -1.0, 0.0, out);
	assert_zero(N - 2);
	assert_zero(quadratic(1.0, -1.0, 0.0, out[0]));
	assert_zero(quadratic(1.0, -1.0, 0.0, out[1]));

}

template<typename FP>
static void testcases(CBRT_SOLVER<FP> cbrt_solver) {
	FP out[3];
	FP val;
	int N;

	N = cbrt_solver(9933.79603159, 7658.09556667, 2059.57642168, 187.33806459, out);
	val = cubic(9933.79603159, 7658.09556667, 2059.57642168, 187.33806459, out[0]);

	N = cbrt_solver(2.04424068e-02, -7.05085740e+04, 2.59523337e+03, 3.18743329e+04, out);
	val = cubic(2.04424068e-02, -7.05085740e+04, 2.59523337e+03, 3.18743329e+04, out[0]);

	N = cbrt_solver(2.14650209e-01, 9.66306444e+04, 8.98525869e+04, 5.26231876e+04, out);
	val = cubic(2.14650209e-01, 9.66306444e+04, 8.98525869e+04, 5.26231876e+04, out[0]);


	N = cbrt_solver(20200.51346605, -9670.18955856, -81125.76266827, -51384.96475652, out);
	val = cubic(20200.51346605, -9670.18955856, -81125.76266827, -51384.96475652, out[0]);
	val = cubic(20200.51346605, -9670.18955856, -81125.76266827, -51384.96475652, out[1]);
	val = cubic(20200.51346605, -9670.18955856, -81125.76266827, -51384.96475652, out[2]);
}


template<typename FP>
static void timing_test_instance(CBRT_SOLVER<FP> cbrt_solver, std::size_t N, int seed)
{

	// Seed with a real random value, if available
	//std::random_device r;

	// Choose a random mean between 1 and 6
	std::default_random_engine e1(seed);
	std::uniform_real_distribution<FP> uniform_dist(-1.0, 1.0);


	auto start = std::chrono::high_resolution_clock::now();
	int N_root_tot = 0;
	for (std::size_t i = 0; i < N; i++)
	{
		FP A = uniform_dist(e1);
		FP B = uniform_dist(e1);
		FP C = uniform_dist(e1);
		FP D = uniform_dist(e1);

		FP out[3];
		N_root_tot += cbrt_solver(A, B, C, D, out);
	}

	auto duration = std::chrono::high_resolution_clock::now() - start;
	auto ms_tot = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
	auto per_poly = std::chrono::duration_cast<std::chrono::nanoseconds>(duration / N);

	std::cout
		<< " | Time (s): " << (ms_tot.count() / 1000.0)
		<< " | Mean (ns): " << per_poly.count()
		<< " | Roots: " << N_root_tot
		<< " |\n";
}


template<typename FP>
static void run_timing_test(CBRT_SOLVER<FP> cbrt_solver, const char* func_name, std::size_t N = 100000000, int seed = 235201124)
{
	std::cout << " | " << func_name << " | Polynomials: " << N << " |\n";
	timing_test_instance(cbrt_solver, N, 235201124);
	timing_test_instance(cbrt_solver, N, 235201124);
	timing_test_instance(cbrt_solver, N, 235201124);
	std::cout << "----\n";
}

int main(void) {

	test_quadratic_zero_double(&qdrtc<double>);

	test_cubic_zero_double(&cubic_roots_qbc<double>);
	test_cubic_zero_double(&cubic_roots<double>);

	testcases(&cubic_roots_qbc<double>);
	testcases(&cubic_roots<double>);

	run_timing_test(&cubic_roots<double>, "cubic");
	run_timing_test(&cubic_roots_qbc<double>, "qbc");

	return 0;

	test_quadratic_zero_double(&quadratic_roots<double>);
	test_cubic_zero_double(&cubic_roots<double>);
	testcases(&cubic_roots<double>);
	return 0;
}
