import unittest
import sys
import os
module_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../out/build/x64-Release/cubic_pybind")
sys.path.insert(0, module_path)
import cubic
import numpy as np
from numpy.polynomial.polynomial import Polynomial

def quad_sum(A, X):
    """
    """
    return A[0] * np.power(X, 2) + A[1] * X + A[2]

def quad_solve(A):
    """
    """
    roots = cubic.quadratic_roots(*A)
    roots = np.sort(roots)
    return roots, quad_sum(A, roots)

def qdrtc_solve(A):
    """
    """
    roots = cubic.qdrtc(*A)
    roots = np.sort(roots)
    return roots, quad_sum(A, roots)


def numpy_solve(A):
    """
    """
    roots = Polynomial(A[::-1]).roots()
    roots_real = np.real(roots[np.isreal(roots)])
    roots_real = np.sort(roots_real)
    return roots_real, quad_sum(A, roots_real)

def is_root(polysum):
    return np.allclose(polysum, 0, atol=1e-7)


def run_algo_comparison(max=1e5, N=10000, N_runs=1, seed=5098359162415):
    """
    """
    N = int(N)
    rng = np.random.default_rng(seed)

    iter = 0
    while (iter < N_runs):
        polys = rng.uniform(-max, max, (N, 3))
        csums, qdrtsums, npsums = [], [],[]
        for i, A in enumerate(polys):
            croots, csum = quad_solve(A)
            qdrtroots, qdrtsum = qdrtc_solve(A)
            nproots, npsum = numpy_solve(A)
                                           
            # Sanity check
            if len(croots) != len(nproots):
                print("Quad | Mismatch in number of roots for polynom %s | Out: %s | Ans: %s" % (i, str(A), str(croots), str(nproots)))
            if len(qdrtroots) != len(nproots):
                print("QDRT | Mismatch in number of roots for polynom %s | Out: %s | Ans: %s" % (i, str(A), str(qdrtroots), str(nproots)))
            ##

            csums.extend(csum)
            qdrtsums.extend(qdrtsum)
            npsums.extend(npsum)

        acsums = np.abs(csums)
        aqdrtsums = np.abs(qdrtsums)
        anpsums = np.abs(npsums)
        print("Algorithm comparison | Run %i | N %i | Max %.2E:" % (iter, N, max))
        # Mean absolute error
        print("Quadratic solver | MAE: %0.16f | MAE Std: %0.16f | EMax: %0.16f" % (np.mean(acsums), np.std(acsums), np.max(acsums)))
        print("QDRT solver | MAE: %0.16f | MAE Std: %0.16f | EMax: %0.16f" % (np.mean(aqdrtsums), np.std(aqdrtsums), np.max(aqdrtsums)))
        print("Numpy solver | MAE: %0.16f | MAE Std: %0.16f | EMax: %0.16f" % (np.mean(anpsums), np.std(anpsums), np.max(anpsums)))
        iter += 1

def verif_qdrt_solver_uniform(qdrt_solver, N=100000):
    rng = np.random.default_rng(5098359162415)


    polys = rng.uniform(-1e0, 1e0, (N, 3))
    csums, npsums = [], []
    for i, A in enumerate(polys):
        croots, csum = qdrt_solver(A)
        nproots, npsum = numpy_solve(A)
                   
        assert np.allclose(nproots, croots), \
            "%i:th failed for polynom %s | Out: %s | Ans: %s" % (i, str(A), str(croots), str(nproots))

        csums.extend(csum)
        npsums.extend(npsum)

    acsums = np.abs(csums)
    anpsums = np.abs(npsums)
    print("Uniform results:")
    print("Quadratic solver | MAE: %0.16f MAE | Std: %0.16f | EMax: %0.16f" % (np.mean(acsums), np.std(acsums), np.max(acsums)))
    print("Numpy solver | MAE: %0.16f MAE | Std: %0.16f | EMax: %0.16f" % (np.mean(anpsums), np.std(anpsums), np.max(anpsums)))

class Unittest(unittest.TestCase):


    def setUpClass():
        pass

    def test_cases_uniform(self):
        return
        verif_qdrt_solver_uniform(quad_solve)
    
    def test_cases_uniform(self):
        return
        verif_qdrt_solver_uniform(qdrtc_solve)

    def test_cmp_algos_max_1e5(self):
        N = 1e5
        N_runs = 3
        run_algo_comparison(1e5, N, N_runs, 5098359162415)
        
    def test_cmp_algos_max_1e0(self):
        return

        N = 1e6
        N_runs = 3
        run_algo_comparison(1e0, N, N_runs, 6702596235672025)




if __name__ == '__main__':
    unittest.main()
