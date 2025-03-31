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
from numpy.random import rand, randn, default_rng
from scipy.linalg import svd, eigh
from scipy.sparse.linalg import svds, norm as spnorm, splu, spsolve_triangular, spsolve, eigsh, gmres
from scipy.sparse import issparse, find, spdiags, block_array, eye as speye

from idaes.core.util.exceptions import BurntToast

def _symmetric_inverse_iteration(H_inv_func, H_shape, n_vec, tol, max_iter):
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

    assert len(H_shape) == 2
    m, n = H_shape
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
    evals = 1 / np.diag(R)
    sort_vec = np.argsort(evals)
    evals = evals[sort_vec]
    evecs = B[:, sort_vec]

    return evals, evecs, converged

def _symmetric_rayleigh_ritz_iteration(H, n_vec, tol, max_iter, seed):
    assert len(H.shape) == 2
    m, n = H.shape
    assert m == n
    rng_obj = default_rng(seed)
    mu = rng_obj.uniform(low=-1e-15, high=1e-15, size=(n_vec,))
    B = rng_obj.standard_normal(size=(m, n_vec))
    B, _ =  qr(B)

    converged = False

    for i in range(max_iter):
        # Calculate how close B[:,k] and mu[k] are to an eigenpair
        # by calculating norm(H @ B[:, k]  - mu[k] * B[:,k]) for all k.
        err = norm(H @ B - B @ np.diag(mu), axis=0)

        # TODO should this tolerance be scaled by matrix size?
        if max(err) < tol:
            converged = True
            print(f"Converged in {i} iterations")
            break

        # Find the index which is furthest from being an eigenvector in
        # order to improve the estimate for that index
        max_err_idx = np.argmax(err)

        try:
            # Calculate 
            Bgrave = spsolve(
                H - mu[max_err_idx] * speye(n),
                B
            )
        except RuntimeError:
            if "Factor is exactly singular" in str(err):
                # mu[max_err_idx] is exactly equal to an
                # eigenvalue, but B[:,max_err_idx] is not
                # a good eigenvector estimate. Therefore
                # perturb mu[max_err_idx] to allow matrix
                # inversion to succeed
                try:
                    Bgrave = spsolve(
                        H - (1e-15 + mu[max_err_idx]) * speye(n),
                        B
                    )
                except RuntimeError as exc:
                    if "Factor is exactly singular" in str(err):
                        # If perturbation fails, we'll assume
                        # something bigger is wrong in the solver
                        raise RuntimeError(
                            "Unable to invert matrix H when shifted by "
                            f"{mu[max_err_idx]}*I. Check whether something "
                            "is wrong with the matrix H's scaling."
                        ) from exc
                    else:
                        raise
            else:
                raise

        Bhat, _ = qr(Bgrave)
        mu, B_tilde = eigh(Bhat.T @ H @ Bhat)
        B =  Bhat @ B_tilde

    return mu, B, converged

def _safe_inverse_splu(H):
    q = H.shape[0]
    try:
        invH = splu(H)
    except RuntimeError as err:
        if "Factor is exactly singular" in str(err):
            H_reg = H + 1e-15 * speye(q)
            try: 
                invH = splu(H_reg)
            except RuntimeError as err2:
                if "Factor is exactly singular" in str(err2):
                    raise BurntToast(
                        "Regularization of singular Gram matrix failed. Please export "
                        "this matrix using scipy.sparse.save_npz and share it with the "
                        "developers of IDAES for troubleshooting. (If this method was "
                        "called through the SVDToolbox, the matrix is stored as the "
                        ".jacobian attribute.)"
                    ) from err2
                else:
                    raise
            print("Gram matrix is singular to machine precision, regularizing it.")
        else:
            raise


    def H_inv_func(B):
        return invH.solve(B)
    
    return H_inv_func

def _aug_eig_processing(
    A,
    evecs,
    zero_tol=1e-2,
    pair_tol=1e-2,
    idp_tol=1e-2,
):
    """
    Takes the output of eigenvalue 
    """
    assert len(A.shape) == 2
    m, n = A.shape
    l = abs(m-n)
    V = evecs[:n, :]
    V_norms = norm(V, axis=0)
    U = evecs[n:, :]
    U_norms = norm(U, axis=0)
    n_vecs = evecs.shape[1]

    # For an index k corresponding to a singular vector pair, we should have
    # V_norms[k] ~= 1/sqrt(2) and U_norms[k] ~= 1/sqrt(2). However, if A
    # is nonsquare, then we'll also have elements of the right null space if
    # m < n or the left null space if m > n. For these elements, we'll have
    # one vector having norm of approximately 1 and the other vector having
    # norm of approximately zero.

    V_zero_indices = np.nonzero(V_norms < zero_tol)[0]
    V_pair_indices = np.nonzero(V_norms > 1/np.sqrt(2) - pair_tol)[0]

    U_zero_indices = np.nonzero(U_norms < zero_tol)[0]
    U_pair_indices = np.nonzero(U_norms > 1/np.sqrt(2) - pair_tol)[0]

    # Check that we don't have any vectors of intermediate size
    # This check is probably redundant now
    assert not np.any(np.logical_and(V_norms > zero_tol, V_norms < 1/np.sqrt(2) - pair_tol))
    assert not np.any(np.logical_and(U_norms > zero_tol, U_norms < 1/np.sqrt(2) - pair_tol))

    if m > n:
        assert len(U_zero_indices) == 0
        assert len(V_zero_indices) + len(V_pair_indices) == n_vecs
        null_space = U[:, V_zero_indices]
        null_norms = U_norms[V_zero_indices]
        null_space = null_space / null_norms
        assert np.all(null_norms > 1 - pair_tol)

        pair_indices = V_pair_indices
        import pdb; pdb.set_trace()


    elif m < n:
        assert len(V_zero_indices) == 0
        assert len(U_zero_indices) + len(U_pair_indices) == n_vecs
        null_space = V[:, U_zero_indices] # Actually the left null space
        null_norms = V_norms[U_zero_indices]
        assert np.all(null_norms > 1 - pair_tol)
        
        pair_indices = U_pair_indices

    else:
        assert len(U_zero_indices) == 0
        assert len(U_pair_indices) == n_vecs
        assert len(V_zero_indices) == 0
        assert len(V_pair_indices) == n_vecs

        pair_indices = U_pair_indices

    U = U[:, pair_indices]
    U_norms = U_norms[pair_indices]
    V = V[:, pair_indices]
    V_norms = V_norms[pair_indices]

    # Numpy broadcasting starts from the last axis and works its way forward
    # so even if U or V are square, this method will correctly normalize
    # the columns, as intended
    V = V / V_norms
    U = U / U_norms

    # We need to include the null space in the orthonormalization process in
    # order to avoid blending betwen small singular values and the null vectors
    if m > n:
        null_space = null_space / null_norms
        U = np.concatenate([U, null_space], axis=1)
    elif m < n:
        null_space = null_space / null_norms
        U = np.concatenate([V, null_space], axis=1)

    # Now we need to deduplicate the singular vectors. The augmented matrix has
    # two eigenvectors corresponding to a singular value. We know that the
    # provided eigenvectors are orthogonal to each other to machine precision,
    # but, unfortunately, that does not mean that the resulting columns of
    # U and V are orthogonal except for duplicate eigenvectors. Furthermore, 
    # we might not get matching pairs of singular vectors.

    U, R_U = qr(U)
    r_U = np.abs(np.diag(R_U))
    lin_idp_vecs = np.nonzero(r_U > idp_tol)[0]
    U = U[:, lin_idp_vecs]

    V, R_V = qr(V)
    r_V = np.abs(np.diag(R_V))
    lin_idp_vecs = np.nonzero(r_V > idp_tol)[0]
    V = V[:, lin_idp_vecs]
    
    # Don't need to worry about columns of U
    # matching up with columns of V, the final
    # SVD will take care of sorting  them.
    # However, we do want the same number of vectors.
    # assert U.shape[1] == V.shape[1]
    U_sub, svals, VT_sub = svd(U.T @ A @ V, full_matrices=True)
    U = U @ U_sub
    V = V @ VT_sub.T
    if m > n:
        null_space = U[:, -l:]
        U = U[:, :-l]
    elif m < n:
        null_space = V[:, -l:]
        V = V[:, :-l]
    # Sort singular values in ascending order
    U = U[:, ::-1]
    V = V[:, ::-1]
    svals = svals[::-1]

    if m == n:
        return U, svals, V
    else:
        return U, svals, V, null_space

def svd_explicit_grammian(
        A,
        number_singular_values: int = 10,
        p: int = 5,
        num_iter: int = 10,
        small_sv_tol: float = 1e-7,
    ):
    """
    Computes smallest singular vectors of the sparse m*n matrix A (m<=n) by explicitly
    forming the Gram matrix A @ A.T or A.T @ A, whichever is smaller, and then carrying
    out inverse iteration on it. Forming this matrix explicitly is typically not 
    recommended because it means that singular values of A that are less than the square 
    root of machine epsilon are indistinguishable, halving the number of significant 
    digits in the results. As a result, this method computes p additional eigenvectors 
    of the Gram matrix, in order to ensure that the entire subspace containing small
    eigenvalues is captured, then projects A into this subspace and conducts a dense
    svd on the much smaller matrix.

    Args:
        A: Scipy sparse array with real entries
        number_singular_values: number of small singular values and vectors to compute
        p: number of extra vectors to oversample Gram matrix with
        num_iter: Number of steps of inverse iteration to conduct on Gram matrix
        small_sv_tol: If the smallest singular value computed for the
            number_singular_values + p exceeds this tolerance, the user is warned that
            the returned singular values and vectors may not accurately represent the
            entire subspace of A associated with singular values less than this value.

    Returns:
        U: m by number_singular_values dense array of left singular vectors
        svals: 1D dense array of number_singular_values singular values
        V: n by number_singular_values dense array of left singular vectors
    """
    assert len(A.shape) == 2
    assert issparse(A)
    assert p >= 0
    m, n = A.shape
    q = min(m, n)
    H_shape = (q, q)

    if m < n:
        H = A @ A.T
    else:
        H = A.T @ A
    H = H.tocsc()

    H_inv_func = _safe_inverse_splu(H)

    # Don't care about the calculated eigenvalues
    # nor about whether the inverse iteration
    # technically converged
    evals, B, _ = _symmetric_inverse_iteration(
        H_inv_func=H_inv_func,
        H_shape=H_shape,
        n_vec=number_singular_values+p,
        tol=1e-10,
        max_iter=num_iter
    )

    if m < n:
        U_tilde, svals, VT = svd(B.T @ A, full_matrices=False)
        U = B @ U_tilde
        V = VT.T
    else:
        U, svals, VT_tilde = svd(A @ B, full_matrices=False)
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



def svd_rayleigh_ritz(
        A,
        number_singular_values: int = 10,
        max_iter: int = 100,
        tol: float = 1e-12,
        seed: int = None,
    ):
    """
    Computes smallest singular vectors of the sparse m*n matrix A (m<=n) by explicitly
    forming the Gram matrix A @ A.T or A.T @ A, whichever is smaller, and then carrying
    out inverse iteration on it. Forming this matrix explicitly is typically not 
    recommended because it means that singular values of A that are less than the square 
    root of machine epsilon are indistinguishable, halving the number of significant 
    digits in the results. As a result, this method computes p additional eigenvectors 
    of the Gram matrix, in order to ensure that the entire subspace containing small
    eigenvalues is captured, then projects A into this subspace and conducts a dense
    svd on the much smaller matrix.

    Args:
        A: Scipy sparse array with real entries
        number_singular_values: number of small singular values and vectors to compute
        p: number of extra vectors to oversample Gram matrix with
        num_iter: Number of steps of inverse iteration to conduct on Gram matrix
        small_sv_tol: If the smallest singular value computed for the
            number_singular_values + p exceeds this tolerance, the user is warned that
            the returned singular values and vectors may not accurately represent the
            entire subspace of A associated with singular values less than this value.

    Returns:
        U: m by number_singular_values dense array of left singular vectors
        svals: 1D dense array of number_singular_values singular values
        V: n by number_singular_values dense array of left singular vectors
    """
    # Should we try to squeeze extra dimensions first? It doesn't look like np.squeeze
    # works on sparse arrays, so we'd need to reshape it.

    # It appears that my version of Scipy does not support sparse arrays of dimension
    # higher than 2, but that such a feature is actively being worked on. Therefore
    # the shape of A first to facilitate testing
    if not len(A.shape) == 2:
        raise ValueError(
            "This method expects a 2D Scipy sparse array-like as input, but was passed "
            f"a {len(A.shape)}D array-like instead."
        )

    if not issparse(A):
        raise ValueError(
            "This method expects a Scipy sparse array-like as an input but was passed "
            "a dense array-like instead. Try using scipy.linalg.svd for a dense SVD method."
        )

    m, n = A.shape
    l = abs(m-n)

    # The augmented matrix has two eigenvectors for every singular triplet, in addition
    # to eigenvectors corresponding to the (left) null space if m<n (m>n). In order to
    # get the required number of singular triplets, we need to exhaust the null space
    # as well as the duplicate singular triplets. 
    n_samples = 2*number_singular_values + l
    A_aug = block_array([[None, A.T],[A, None]])

    _, B, converged = _symmetric_rayleigh_ritz_iteration(
        A_aug,
        n_samples,
        tol=tol,
        max_iter=max_iter,
        seed=seed
    )

    if not converged:
        raise RuntimeError(
            "Rayleigh-Ritz iteration did not converge! Consider increasing the tolerance or "
            "maximum number of iterations."
        )

    # try:
    if m == n:
        U, svals, V = _aug_eig_processing(A, B)
    else:
        U, svals, V, null = _aug_eig_processing(A, B)
    # except AssertionError as exc:
    #     raise RuntimeError(
    #         "Categorization of singular vectors failed despite Rayleigh-Ritz iteration "
    #         "converging. Please let the IDAES dev team know about this failure so that "
    #         "appropriate tolerances for categorization can be chosen."
    #     ) from exc

    # Singular values already in ascending order, so just take the number that we want
    U = U[:, :number_singular_values]
    V = V[:, :number_singular_values]
    svals = svals[:number_singular_values]

    if m == n:
        return U, svals, V
    else:
        return U, svals, V, null
    

def svd_inverse_aug(
        A,
        number_singular_values: int = 10,
        p: int = 5,
        max_iter: int = 100,
        tol: float = 1e-14
    ):
    """
    Computes smallest singular vectors of the sparse m*n matrix A (m<=n) by explicitly
    forming the Gram matrix A @ A.T or A.T @ A, whichever is smaller, and then carrying
    out inverse iteration on it. Forming this matrix explicitly is typically not 
    recommended because it means that singular values of A that are less than the square 
    root of machine epsilon are indistinguishable, halving the number of significant 
    digits in the results. As a result, this method computes p additional eigenvectors 
    of the Gram matrix, in order to ensure that the entire subspace containing small
    eigenvalues is captured, then projects A into this subspace and conducts a dense
    svd on the much smaller matrix.

    Args:
        A: Scipy sparse array with real entries
        number_singular_values: number of small singular values and vectors to compute
        p: number of extra vectors to oversample Gram matrix with
        num_iter: Number of steps of inverse iteration to conduct on Gram matrix
        small_sv_tol: If the smallest singular value computed for the
            number_singular_values + p exceeds this tolerance, the user is warned that
            the returned singular values and vectors may not accurately represent the
            entire subspace of A associated with singular values less than this value.

    Returns:
        U: m by number_singular_values dense array of left singular vectors
        svals: 1D dense array of number_singular_values singular values
        V: n by number_singular_values dense array of left singular vectors
    """
    assert len(A.shape) == 2
    assert issparse(A)
    assert p >= 0
    m, n = A.shape
    q = min(m, n)
    l = abs(m-n)
    n_samples = 2*number_singular_values + l
    A_aug = block_array([[None, A.T],[A, None]])

    shift = 1e-15 *rand()

    try:
        invH = splu(A_aug - shift * speye(m + n))
    except RuntimeError as err:
        if "Factor is exactly singular" in str(err):
            # Get a new shift and try again. The chance of choosing one
            # singular value is miniscule, the chance of choosing two in a row
            # is miniscule squared.
            shift = 1e-15 *rand()
            try: 
                invH = splu(A_aug - shift * speye(m + n))
            except RuntimeError as err2:
                if "Factor is exactly singular" in str(err2):
                    raise BurntToast(
                        "Failed to shift-invert augmented matrix. Either the random number "
                        "generator just exactly chose two singular values in a row, or "
                        "some other issue is causing matrix inversion to fail. Please export "
                        "this matrix using scipy.sparse.save_npz and share it with the "
                        "developers of IDAES for troubleshooting. (If this method was "
                        "called through the SVDToolbox, the matrix is stored as the "
                        ".jacobian attribute.)"
                    ) from err2
                else:
                    raise
        else:
            raise
    

    def H_inv_func(B):
        return invH.solve(B)

    evals, B, converged = _symmetric_inverse_iteration(
        H_inv_func=H_inv_func,
        H_shape=A_aug.shape,
        n_vec=n_samples,
        tol=1e-14,
        max_iter=max_iter,
    )
    if not converged:
        raise RuntimeError

    if m == n:
        U, svals, V = _aug_eig_processing(A, B)
    else:
        U, svals, V, null = _aug_eig_processing(A, B)

    # Singular values already in ascending order, so just take the number that we want
    U = U[:, :number_singular_values]
    V = V[:, :number_singular_values]
    svals = svals[:number_singular_values]

    if m == n:
        return U, svals, V
    else:
        return U, svals, V, null