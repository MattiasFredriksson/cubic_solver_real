# Cubic and quartic real root finder

C++ implementation for computing real roots of a cubic or quadratic equation (3rd or 2nd order polynomial). [Src code](https://github.com/MattiasFredriksson/cubic_solver_real/blob/master/Cubic/cubic_lib/src/cubic.cpp) contains a "closed form" solver and the numerical root finding algorithm "QBC" (see (W. Kahan notes)[https://people.eecs.berkeley.edu/~wkahan/Math128/Cubic.pdf]) utilizing (Newton's method)[https://en.wikipedia.org/wiki/Newton%27s_method]. The closed form solver utilizes both (Cardano's method)[https://en.wikipedia.org/wiki/Cubic_equation#Cardano's_method] and the trigonometric method. Closed form is loosely based on https://github.com/tatwood/solvecubic but adapted/optimized to reduce MAE (mean absolute error) and runtime for computing roots in 64-bit precision. 

## Comparison to Numpy.roots()

Tests consist of 1E6 polynomial functions with coefficients sampled from a uniform distribution with fixed seed. Each test is repeated three times and statistics is computed by evaluating both algorithms on the same coefficients. Outcome will be hardware/software dependent.

Absolute error is computed as e = |f(x_r)| where f(x) is the polynomial function and x_r is any real root computed for f(x) = 0.

### Cubic equation

#### Coefficients uniformly drawn in the range [1e-0, 1e0)

Algo. | Test No. | MAE | Std | Max 
--- | --- | --- | --- | --- 
Cubic | 0 |  0.0000000000199807 | 0.0000000105674862 | 0.0000115097497277
Numpy | 0 |  0.0000000000293169 | 0.0000000169255077 | 0.0000190077410858
Cubic  | 1 |  0.0000000001029146 | 0.0000000763570113 | 0.0000895105459200
Numpy  | 1 |  0.0000000003381824 | 0.0000002904316352 | 0.0003336512873353
Cubic  | 2 |  0.0000000018835925 | 0.0000022526917038 | 0.0026981010554332
Numpy  | 2 |  0.0000000018860261 | 0.0000022526918885 | 0.0026981010554332


#### Coefficients uniformly drawn in the range [1e-5, 1e5)

Algo. | Test No. | MAE | Std | Max 
--- | --- | --- | --- | --- 
Cubic  | 0 |  0.0000007717178477 | 0.0002363548283299 | 0.1472778703464428
Numpy  | 0 | 0.0000042240897500 | 0.0028116804269617 | 3.1311602018395206
Cubic  | 1 |  0.0000072913011206 | 0.0051505376408809 | 4.7461626961303409
Numpy  |  1 | 0.0000237523910556 | 0.0185939038980810 | 20.7461779549194034
Cubic  | 2 |  0.0001678360884921 | 0.2004852254811660 | 240.1257147381838877
Numpy  | 2 |  0.0001902332604939 | 0.2269926509587637 | 271.8742242266598623


### Quadratic equation

#### Coefficients uniformly drawn in the range [1e-0, 1e0)

Algo. | Test No. | MAE | Std | Max 
--- | --- | --- | --- | --- 
Quadratic | 0 | 0.0000000000000003 | 0.0000000000000319 | 0.0000000000316736
Numpy | 0 | 0.0000000000000004 | 0.0000000000000322 | 0.0000000000316736
Quadratic | 1 | 0.0000000000000004 | 0.0000000000000597 | 0.0000000000634905
Numpy | 1 | 0.0000000000000004 | 0.0000000000000598 | 0.0000000000634905
Quadratic | 2 | 0.0000000000000004 | 0.0000000000001256 | 0.0000000001389600
Numpy | 2 | 0.0000000000000004 | 0.0000000000001255 | 0.0000000001389600


#### Coefficients uniformly drawn in the range [1e-5, 1e5)

Algo. | Test No. | MAE | Std | Max 
--- | --- | --- | --- | --- 
Quadratic | 0 | 0.0000000000618843 | 0.0000000214969506 | 0.0000229584111366
Numpy | 0 | 0.0000000000638890 | 0.0000000215317834 | 0.0000229584111366
Quadratic | 1 | 0.0000000000345022 | 0.0000000021549128 | 0.0000013075477909
Numpy | 1 | 0.0000000000368232 | 0.0000000023601247 | 0.0000013075477909
Quadratic | 2 | 0.0000000000453204 | 0.0000000124810741 | 0.0000133562862175
Numpy | 2 | 0.0000000000490860 | 0.0000000126699302 | 0.0000133562862175
