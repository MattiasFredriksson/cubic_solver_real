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

/* Simple vector based bind function.
*/
std::vector<double> cubic_roots_bind(double a, double b, double c, double d) {

	double res[3];
	int real_roots = cubic_roots(a, b, c, d, res);

	std::vector<double> out(real_roots);
	for (int i = 0; i < real_roots; i++) {
		out[i] = res[i];
	}
	return out;
}

/* Simple vector based bind function.
*/
std::vector<double> quadratic_roots_bind(double a, double b, double c) {

	double res[2];
	int real_roots = quadratic_roots(a, b, c, res);

	std::vector<double> out(real_roots);
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

	m.def("cubic_roots", &cubic_roots_bind, R"pbdoc(
        Compute the real roots for the cubic equation.
    )pbdoc");
	m.def("quadratic_roots", &quadratic_roots_bind, R"pbdoc(
        Compute the real roots for the quadratic equation.
    )pbdoc");

#ifdef VERSION_INFO
	m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
	m.attr("__version__") = "dev";
#endif
}