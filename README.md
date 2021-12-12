# Cubic and quartic real root finder

C++ implementation for computing real roots of a cubic or quadratic equation (3rd or 2nd order polynomial). [Src code](https://github.com/MattiasFredriksson/cubic_solver_real/blob/master/Cubic/cubic_lib/src/cubic.cpp) contains a "closed form" solver and the numerical root finding algorithm "QBC" (see [W. Kahan notes](https://people.eecs.berkeley.edu/~wkahan/Math128/Cubic.pdf)) based on [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method). The closed form implementation utilizes both [Cardano's method](https://en.wikipedia.org/wiki/Cubic_equation#Cardano's_method) and the trigonometric method. Implementation of the closed form solver is loosely based on https://github.com/tatwood/solvecubic but adapted to reduce MAE (mean absolute error) and runtime for computing roots in 64-bit precision. 


## Performance

Execution time averaged over 1e6 calls for generating coefficients drawn from a uniform distribution and solving the cubic equation using the closed form (cubic) and QBC algorithms.

Algo. | Average time (ns)
--- | --- 
Cubic | 198
QBC |  282

## Comparison to Numpy.roots()

Tests consist of 1E6 polynomial functions with coefficients sampled from a uniform distribution with fixed seed. Each test is repeated three times and statistics is computed by evaluating both algorithms on the same coefficients. Outcome will be hardware/software dependent.

Absolute error is computed as e = |f(x_r)| where f(x) is the polynomial function and x_r is any real root computed for f(x) = 0.

### Coefficients uniformly drawn in the range [1e-0, 1e0)

Algo. | Test No. | MAE | Std | Max 
--- | --- | --- | --- | --- 
Cubic | 0 | 0.0000000000199807 | 0.0000000105674862 | 0.0000115097497277
QBC | 0 | 0.0000000000153321 | 0.0000000103069035 | 0.0000115097497277
Numpy | 0 | 0.0000000000293169 | 0.0000000169255077 | 0.0000190077410858
Cubic |  1 | 0.0000000001029146 | 0.0000000763570113 | 0.0000895105459200
QBC |  1 | 0.0000000000333290 | 0.0000000277272149 | 0.0000325598247877
Numpy |  1 | 0.0000000003381824 | 0.0000002904316352 | 0.0003336512873353
Cubic | 2 | 0.0000000018835925 | 0.0000022526917038 | 0.0026981010554332
QBC | 2 | 0.0000000008438846 | 0.0000010087043156 | 0.0012081484789055
Numpy | 2 |0.0000000018860261 | 0.0000022526918885 | 0.002698101055433


### Coefficients uniformly drawn in the range [1e-5, 1e5)

Algo. | Test No. | MAE | Std | Max 
--- | --- | --- | --- | --- 
Cubic   | 0 |  0.0000007717178477 | 0.0002363548283299 | 0.1472778703464428
QBC   | 0 | 0.0000005225410921 | 0.0002198728449618 | 0.1472781087650219
Numpy   | 0 | 0.0000042240897500 | 0.0028116804269617 | 3.1311602018395206
Cubic   | 1 | 0.0000072913011206 | 0.0051505376408809 | 4.7461626961303409
QBC  | 1 | 0.0000043226077346 | 0.0039805146698185 | 4.7461703255248722
Numpy   | 1 | 0.0000237523910556 | 0.0185939038980810 | 20.7461779549194034
Cubic  | 2 | 0.0001678360884921 | 0.2004852254811660 | 240.1257147381838877
QBC  | 2 | 0.0001677691287439 | 0.2004852521819428 | 240.1257452557620127
Numpy  | 2 | 0.0001902332604939 | 0.2269926509587637 | 271.8742242266598623




