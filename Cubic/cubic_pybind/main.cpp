// DTW_python.cpp : Defines the entry point for the application.
//
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

#include <cubic/cubic.h>


#ifndef PROJECT_NAME_DEF
#define PROJECT_NAME_DEF no_name_test
#endif

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

template<typename FP>
using CBRT_SOLVER_BIND = std::vector<FP>(*)(FP, FP, FP, FP);

/* Simple vector based bind function.
*/
template<typename FP, CBRT_SOLVER<FP> solver>
std::vector<FP> cubic_roots_bind(FP a, FP b, FP c, FP d) {

	FP res[3];
	int real_roots = solver(a, b, c, d, res);

	std::vector<FP> out(real_roots);
	for (int i = 0; i < real_roots; i++) {
		out[i] = res[i];
	}
	return out;
}

/* Simple vector based bind function.
*/
template<typename FP, QDRT_SOLVER<FP> solver>
std::vector<FP> quadratic_roots_bind(FP a, FP b, FP c) {

	FP res[2];
	int real_roots = solver(a, b, c, res);

	std::vector<FP> out(real_roots);
	for (int i = 0; i < real_roots; i++) {
		out[i] = res[i];
	}
	return out;
}


namespace py = pybind11;

PYBIND11_MODULE(PROJECT_NAME_DEF, m) {
	m.doc() = R"pbdoc(
        Cubic solver pybinds
        -----------------------
        .. currentmodule:: cubic
        .. autosummary::
           :toctree: _generate
           cubic_roots
		   quadratic_roots
    )pbdoc";

	m.def("cubic_roots", &cubic_roots_bind<double, &cubic_roots<double>>, R"pbdoc(
        Compute the real roots for the cubic equation.
    )pbdoc");
	m.def("quadratic_roots", &quadratic_roots_bind<double, quadratic_roots<double>>, R"pbdoc(
        Compute the real roots for the quadratic equation.
    )pbdoc");


	m.def("cubic_roots_qbc", &cubic_roots_bind<double, &cubic_roots_qbc<double>>, R"pbdoc(
        Compute the real roots for the cubic equation.
    )pbdoc");

	m.def("qdrtc", &quadratic_roots_bind<double, &qdrtc<double>>, R"pbdoc(
        Compute the real roots for the cubic equation.
    )pbdoc");

#ifdef VERSION_INFO
	m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
	m.attr("__version__") = "dev";
#endif
}