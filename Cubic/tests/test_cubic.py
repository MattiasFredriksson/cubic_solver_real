import unittest
import sys
import os
module_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../out/build/x64-Release/cubic_pybind")
sys.path.insert(0, module_path)
import cubic
import numpy as np
from numpy.polynomial.polynomial import Polynomial

def cubic_sum(A, X):
    """
    """
    return A[0] * np.power(X, 3) + A[1] * np.power(X, 2) + A[2] * X + A[3]

def cubic_solve(A):
    """
    """
    roots = cubic.cubic_roots(*A)
    roots = np.sort(roots)
    return roots, cubic_sum(A, roots)

def numpy_solve(A):
    """
    """
    roots = Polynomial(A[::-1]).roots()
    roots_real = np.real(roots[np.isreal(roots)])
    roots_real = np.sort(roots_real)
    return roots_real, cubic_sum(A, roots_real)

def fqs_solve(A):
    """
    """
    roots = np.array(fqs_single_cubic(*A))
    roots_real = np.real(roots[np.isreal(roots)])
    roots_real = np.sort(roots_real)
    return roots_real, cubic_sum(A, roots_real)

def is_root(polysum):
    return np.allclose(polysum, 0, atol=1e-7)

def run_algo_comparison(max=1e5, N=10000, N_runs=1, seed=5098359162415):
    """
    """
    N = int(N)

    rng = np.random.default_rng(seed)

    iter = 0
    while (iter < N_runs):
        polys = rng.uniform(-max, max, (N, 4))
        csums, npsums = [], []
        for i, A in enumerate(polys):
            croots, csum = cubic_solve(A)
            nproots, npsum = numpy_solve(A)
                
            # Sanity check
            if len(croots) != len(nproots):
                print("Mismatch in number of roots for polynom %s | Out: %s | Ans: %s" % (str(A), str(croots), str(nproots)))

            csums.extend(csum)
            npsums.extend(npsum)

        acsums = np.abs(csums)
        anpsums = np.abs(npsums)
        print("Algorithm comparison | Run %i | N %i | Max %.2E:" % (iter, N, max))
        # Mean absolute error
        print("Cubic solver | MAE: %0.16f | MAE Std: %0.16f | EMax: %0.16f" % (np.mean(acsums), np.std(acsums), np.max(acsums)))
        print("Numpy solver | MAE: %0.16f | MAE Std: %0.16f | EMax: %0.16f" % (np.mean(anpsums), np.std(anpsums), np.max(anpsums)))
        iter += 1

class Unittest(unittest.TestCase):


    def setUpClass():
        pass

    def test_cases_uniform(self):
        return

        N = 100000
        rng = np.random.default_rng(5098359162415)


        polys = rng.uniform(-1e5, 1e5, (N, 4))
        csums, npsums = [], []
        for i, A in enumerate(polys):
            croots, csum = cubic_solve(A)
            nproots, npsum = numpy_solve(A)
                   
            assert np.allclose(nproots, croots), \
               "%i:th failed for polynom %s | Out: %s | Ans: %s" % (i, str(A), str(croots), str(nproots))

            csums.extend(csum)
            npsums.extend(npsum)

        acsums = np.abs(csums)
        anpsums = np.abs(npsums)
        print("Uniform results:")
        print("Cubic solver | MAE: %0.16f MAE | Std: %0.16f | EMax: %0.16f" % (np.mean(acsums), np.std(acsums), np.max(acsums)))
        print("Numpy solver | MAE: %0.16f MAE | Std: %0.16f | EMax: %0.16f" % (np.mean(anpsums), np.std(anpsums), np.max(anpsums)))
        
    
    def test_cmp_algos_max_1e5(self):

        N = int(1e6)
        N_runs = 3
        run_algo_comparison(1e5, N, N_runs)
        
    def test_cmp_algos_max_1e0(self):

        N = int(1e6)
        N_runs = 3
        run_algo_comparison(1e0, N, N_runs)


def fqs_single_cubic(a0, b0, c0, d0):
    ''' Analytical closed-form solver for a single cubic equation
    (3rd order polynomial), gives all three roots.
    Parameters
    ----------
    a0, b0, c0, d0: array_like
        Input data are coefficients of the Cubic polynomial::
            a0*x^3 + b0*x^2 + c0*x + d0 = 0
    Returns
    -------
    roots: tuple
        Output data is a tuple of three roots of a given polynomial.

    https://github.com/NKrvavica/fqs/blob/master/fqs.py
    '''
    import math
    import cmath


    ''' Reduce the cubic equation to to form:
        x^3 + a*x^2 + b*x + c = 0 '''
    a, b, c = b0 / a0, c0 / a0, d0 / a0

    # Some repeating constants and variables
    third = 1. / 3.
    a13 = a * third
    a2 = a13 * a13
    sqr3 = math.sqrt(3)

    # Additional intermediate variables
    f = third * b - a2
    g = a13 * (2 * a2 - b) + c
    h = 0.25 * g * g + f * f * f

    def cubic_root(x):
        ''' Compute cubic root of a number while maintaining its sign'''
        if x.real >= 0:
            return x ** third
        else:
            return -(-x) ** third

    if f == g == h == 0:
        r1 = -cubic_root(c)
        return r1, r1, r1

    elif h <= 0:
        j = math.sqrt(-f)
        k = math.acos(-0.5 * g / (j * j * j))
        m = math.cos(third * k)
        n = sqr3 * math.sin(third * k)
        r1 = 2 * j * m - a13
        r2 = -j * (m + n) - a13
        r3 = -j * (m - n) - a13
        return r1, r2, r3

    else:
        sqrt_h = cmath.sqrt(h)
        S = cubic_root(-0.5 * g + sqrt_h)
        U = cubic_root(-0.5 * g - sqrt_h)
        S_plus_U = S + U
        S_minus_U = S - U
        r1 = S_plus_U - a13
        r2 = -0.5 * S_plus_U - a13 + S_minus_U * sqr3 * 0.5j
        r3 = -0.5 * S_plus_U - a13 - S_minus_U * sqr3 * 0.5j
        return r1, r2, r3


if __name__ == '__main__':
    unittest.main()
