#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains specialized linear algebra routines not found in Numpy or Scipy
"""

__author__ = "Douglas Allan"

import numpy as np
from numpy.linalg import norm, qr
from numpy.random import randn
from scipy.linalg import svd
from scipy.sparse.linalg import svds, norm as spnorm, splu, spsolve_triangular, spsolve
from scipy.sparse import issparse, find, spdiags, block_array

def _symmetric_inverse_iteration(H, H_inv_func, n_vec, tol, max_iter):
    """
    Function to perform simultaneous inverse iteration on a real
    symmetric matrix in order to find the eigenvalues and eigenvectors
    smallest in absolute value (i.e., closest to zero).

    Args:
        H: Real, symmetric n by n matrix H
        H_inv_func: Function to apply inv(H) to a n x n_vec matrix B
            returning B_new = inv(H) @ B
        n_vec: Number of eigenvectors and eigenvalues to calculate
        max_iter: maximum number of iterations for inverse iteration

    Returns:
        evals: 1D array of n_vec eigenvalues
        evecs: n by n_vec 2D array of eigenvectors
        converged: Boolean flag of whether the convergence criteria
            was met after max_iter iterations
    """

    assert len(H.shape) == 2
    m, n = H.shape
    assert m == n  
    # Get some random vectors to start inverse iteration
    # Supposedly normal distributions are better than uniform distributions
    # for this purpose to avoid "corner effects" in high dimensions
    # (Gleaned from random Stack Exchange comment)
    B = randn(m, n_vec)
    # Use QR factorization to get orthonormal basis for random vectors
    B, _ = qr(B, mode="reduced")

    converged = False

    for i in range(max_iter):
        B_new = H_inv_func(B)

        # Orthonormalize new basis
        B_new, R = qr(B_new, mode="reduced")

        # Project old basis set into the subspace spanned by new basis set
        # Note that parentheses are very important. B_new.T @ B is size
        # n_vec * n_vec, whereas B_new @ B_new.T is m * m.
        B_proj = B_new @ (B_new.T @ B)

        # Subtract this projection of B into the new subspace from the original
        # version of X. If B == B_new, then this difference should be zero.

        # Alternative, if R converges to a nearly diagonal matrix, then we
        # say we've converged. np.triu(R, 1) gets the portion of R above the
        # the main diagonal, then we normalize by the largest diagonal element
        err = min(
            # Need to avoid growth of Frobenius norm with matrix dimension
            # so add 
            norm(B - B_proj, ord="fro") / np.sqrt(n_vec), 
            norm(np.triu(R,1)/max(np.diag(R)))
        )

        B = B_new

        if err < tol:
            print(f"Converged in {i} iterations")
            converged = True
            break

    # If we've converged, the eigenvalues are the diagonal entries of R
    evals = np.diag(R)
    sort_vec = np.argsort(evals)
    evals = evals[sort_vec]
    evecs = B[:, sort_vec]

    return evals, evecs, converged

    


def svd_explicit_normal(A, number_singular_values=10, p=5, num_iter=10, small_sv_tol=1e-7):
    """
    Computes the n_vec smallest singular vectors of the Scipy sparse m*n matrix A (m<=n).
    Because this method is typically called on singular or near-singular A, there can be 
    significant error in this calcuation. Nevertheless, the resulting matrix U will contain
    n_vec orthonormal vectors whose left-multiplication into A, U.T @ A, results in a matrix
    of extremely small magnitude---in other words, they are left near-null vectors of A.
    Typically such vectors are sufficient for model diagnostic tools.

    Note that if singular values are close to each other, convergence can be impaired and
    the resulting singular vector estimate may be a linear combination of the singular 
    vectors associated with the other singular values.

    For example, for a matrix with condition number of 5.7e17, this method returned a 
    vector associated with an estimated singular value of 9.31e-5. Comparison to the 
    results of a dense SVD revealed the vector to be primarily a linear combination 
    of 11 singular vectors with associated singular values ranging from 9.19e-5 to 9.82e-5.
    However, for model diagnostic purposes, the singular vector estimate is adequate.
    """
    assert len(A.shape) == 2
    assert issparse(A)
    assert p >= 0
    m, n = A.shape

    if m < n:
        H = A @ A.T
    else:
        H = A.T @ A
    H = H.tocsc()

    invH = splu(H)
    def H_inv_func(B):
        return invH.solve(B)

    # Don't care about the calculated eigenvalues
    # nor about whether the inverse iteration
    # technically converged
    _, B, _ = _symmetric_inverse_iteration(
        H,
        H_inv_func=H_inv_func,
        n_vec=number_singular_values+p,
        tol=1e-10,
        max_iter=num_iter
    )

    if m < n:
        U_tilde, svals, VT = svd(B.T @ A)
        U = B @ U_tilde
    else:
        U, svals, VT_tilde = svd(A @ B)
        V = B @ VT_tilde.T

    if svals[0] < small_sv_tol:
        print(
            f"Dimension of the near-singular subspace exceeds {number_singular_values} + {p}. "
            "The returned vectors may not contain the smallest singular values. "
            "To be sure that the smallest singular values are included, call this "
            f"method again with p>{p} or call it on a submodel with input variables fixed."
        )

    if p > 0:
        U = U[:,:p-1:-1]
        svals = svals[:p-1:-1]
        V = V[:, :p-1:-1]
    else:
        U = U[:,::-1]
        svals = svals[::-1]
        V = V[:, ::-1]

    return U, svals, V