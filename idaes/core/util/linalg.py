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
from numpy.linalg import norm
from numpy.random import rand, randn, default_rng
from scipy.linalg import svd, eigh, qr
from scipy.sparse.linalg import svds, norm as spnorm, splu, spsolve_triangular, spsolve, eigsh, gmres
from scipy.sparse import issparse, find, spdiags, block_array, eye as speye

from idaes.core.util.exceptions import BurntToast

def _symmetric_rayleigh_ritz_iteration(H, n_vec, tol, max_iter, seed):
    assert len(H.shape) == 2
    m, n = H.shape
    assert m == n
    rng_obj = default_rng(seed)
    mu = rng_obj.uniform(low=-1e-15, high=1e-15, size=(n_vec,))
    B = rng_obj.standard_normal(size=(m, n_vec))
    B, _ =  qr(B, mode="economic")

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

        Bhat, _ = qr(Bgrave, mode="economic")
        mu, B_tilde = eigh(Bhat.T @ H @ Bhat)
        B =  Bhat @ B_tilde

    return mu, B, converged


    
def _aug_eig_processing(
    A,
    evecs,
    zero_tol=1e-2,
):
    """
    Takes the output of subspace iteration on the augmented matrix [[0, A.T], [A, 0]]
    for some m by n matrix A, and obtains estimates for U, svals, V such that 
    U.T @ A @ V == diag(svals). If m < n (m > n), then a basis for the (left) null
    space is also computed.

    For every singular triplet (sigma, u, v), the augmented matrix has two eigentriplets:
    (sigma, (v, u)) and (-sigma, (-v, u)). If m!=n, there are |m-n| eigenvectors corresponding
    to the null space. If m > n, then these vectors are of the form (0, (0, w)), and if m < n
    they are of the form (0, (w, 0)). We take eigenvectors of the augmented matrix as an input,
    but there is no telling how many singular vectors they correspond to: sometimes the method
    converges to both (sigma, (v, u)) and (-sigma, (-v, u)) and sometimes it only converges to
    a single one. Furthermore, if some sigma is close to 0, we can end up with blending between
    the null space and the space corresponding to the small singular value. Finally, while the
    eigenvectors may be orthogonal to machine precision, the resulting singular vectors may not
    be orthogonal due to error cancellation between the u and v components: for (v_1, u_1) and
    (v_2, u_2), we may have v_1.T @ v_2 + u_1.T @ u_2 ~= 1e-14, but nevertheless have
    v_1.T @ v_2 ~= 1e-7 and u_1.T @ u_2 ~= 1e-7.

    The solution to this very messy situation is, as usual, the QR factorization. We partition
    the eigenvectors into U and V components, then perform a rank-revealing QR in order to
    obtain orthonormal bases of those sets of vectors. The diagonal elements of R are measures
    of how strongly a basis vector is represented in the U and V components of evecs. In an
    ideal world, if evecs contained q singular vectors and p null vectors, we would have
    R[k,k] == 1 for k<=p-1 (corresponding to the null vectors), R[k, k] = 1/sqrt(2) for 
    p <= k <= p+q-1 (corresponding to the singular vectors), and R[k,k] ~= 1e-15 for k > p+q.
    However, we often get a spectrum of values. 
    """
    assert len(A.shape) == 2
    m, n = A.shape
    p = abs(m-n)
    assert (evecs.shape[1] - p) % 2 == 0
    n_vecs = (evecs.shape[1] - p) // 2

    V = evecs[:n, :]
    U = evecs[n:, :]

    # We expect 

    U, R_U, _ = qr(U, mode="economic", pivoting=True)
    r_U = np.abs(np.diag(R_U))

    V, R_V, _ = qr(V, mode="economic", pivoting=True)
    r_V = np.abs(np.diag(R_V))

    n_U = n_vecs
    n_V = n_vecs
    for k in range(n_vecs, evecs.shape[1]):
        if m > n:
            if r_U[k] > zero_tol:
                if r_V[k] > zero_tol:
                    n_U += 1
                    n_V +=1
                elif n_U - n_V < p:
                    n_U +=1
        else:
            if r_V[k] > zero_tol:
                if r_U[k] > zero_tol:
                    n_U += 1
                    n_V +=1
                     
                elif n_V - n_U < p:
                    n_V +=1
    
    # import pdb; pdb.set_trace()
    U = U[:, :n_U]
    V = V[:, :n_V]
            

    # Want to make sure that the number of null vector candidates
    # is the expected size
    # This assertion should no longer be necessary
    l_hat = abs(U.shape[1] - V.shape[1])
    assert l_hat <= p
    U_sub, svals, VT_sub = svd(U.T @ A @ V, full_matrices=True)
    U = U @ U_sub
    V = V @ VT_sub.T

    if m > n:
        null_space = U[:, -l_hat:]
        # Should this be A.nnz instead of m*n?
        assert norm(null_space.T @ A) <= np.sqrt(m*n) * 1e-15
        U = U[:, :-l_hat]
    elif m < n:
        null_space = V[:, -l_hat:]
        # Should this be A.nnz instead of m*n?
        assert norm(A @ null_space) <= np.sqrt(m*n) * 1e-15
        V = V[:, :-l_hat]
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

