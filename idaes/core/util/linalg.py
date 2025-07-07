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
from numpy.random import default_rng
from scipy.linalg import svd, eigh, qr
from scipy.sparse.linalg import splu
from scipy.sparse import issparse, block_array, eye as speye

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def _symmetric_rayleigh_ritz_iteration(H, n_vec, tol, max_iter, seed=None):
    """
    Function to perform simultaneous Rayleigh quotient iteration on the
    real symmetric matrix H.

    Start out with initial eigenvalue estimates chosen randomly from
    -1e-15 to 1e-15 in order to encourage convergence towards the
    eigenvaleus of smallest magnitude, while the original eigenvector
    estimates are chosen as an orthogonalization of a random matrix
    whose values were chosen from the standard Gaussian distribution.

    Then this iterative procedure is followed:
        1) Compute the error in the eigenvalue approximations
            err[j] = abs(H @ B[:, j] - mu[j] * B[:,j]) for all j
        1a) If err[j] < tol for all j, terminate iteration
        2) Choose j* such that err[j*] is maximized.
        3) Calculate Bgrave = (H - mu[j*]) ** -1 @ B via a sparse
            LU decomposition
        4) Calculate an orthonormal basis Bhat for Bgrave using a
            QR factorization
        5) Calculate the Ritz values of H with respect to the
            subspace Bhat via a symmetric eigendecomposition of
            Bhat.T @ H @ Bhat and construct the Ritz vectors
            from Bhat and the eigenvectors of Bhat.T @ H @ Bhat
        6) Set mu and B to equal the Ritz values and vectors,
            respectively, then return to (1).

    The scalar version of this algorithm, in which B is a single
    vector and mu is a scalar, is known as Rayleigh quotient
    iteration, and is known to be cubically convergent for
    symmetric matrices. I (Doug A.) have not seen this vectorized
    extension anywhere in the literature. General attention
    has been given to the (implicitly shifted) QR algorithm for
    dense matrices and Krylov subspace methods for sparse matrices.
    However, dense methods do not scale well and the Lanczos
    algorithm (a Krylov method) provided by ARPACK through Scipy's
    sparse.linalg.eigsh does not provide high quality eigenvector
    approximations (or at least high enough quality approximations
    to use as part of a sparse SVD routine), even when used in
    shift-invert mode on the augmented matrix.

    For finding small eigenvalues, the methods provided in Section 5
    of Sleijpen and Van Der Vorst (2000) and Section 4.4 of
    Hochstenbach (2001) may provide performant alternatives, because
    they can be implemented in a vectorized form in Python. The
    Lanczos variant suggested by Kokiopoulou et al. (2004) may
    be worth attention, but may have to be implemented in a lower
    level language to be performant.

    Sleijpen, G.L.G., Van Der Vorst, H.A., 2000. A Jacobi--Davidson
    Iteration Method for Linear Eigenvalue Problems. SIAM Rev. 42,
    267–293. https://doi.org/10.1137/S0036144599363084

    Hochstenbach, M.E., 2001. A Jacobi--Davidson Type SVD Method.
    SIAM J. Sci. Comput. 23, 606–628.
    https://doi.org/10.1137/S1064827500372973

    Kokiopoulou, E., Bekas, C., Gallopoulos, E., 2004. Computing
    smallest singular triplets with implicitly restarted Lanczos
    bidiagonalization. Applied Numerical Mathematics 49, 39–61.
    https://doi.org/10.1016/j.apnum.2003.11.011



    Args:
        H: Real, symmetric n by n matrix H
        n_vec: Number of eigenvectors and eigenvalues to calculate
        tol: Tolerance used to decide to terminate iteration.
        max_iter: maximum number of iterations

    Returns:
        evals: 1D array of n_vec eigenvalues
        evecs: n by n_vec 2D array of eigenvectors
        converged: Boolean flag of whether the convergence criteria
            was met after max_iter iterations
    """
    assert len(H.shape) == 2
    m, n = H.shape
    assert m == n
    rng_obj = default_rng(seed)
    mu = rng_obj.uniform(low=-1e-15, high=1e-15, size=(n_vec,))
    B = rng_obj.standard_normal(size=(m, n_vec))
    B, _ = qr(B, mode="economic")
    # import pdb; pdb.set_trace()

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
            invH_shift = splu(H - mu[max_err_idx] * speye(n, format="csc"))
        except RuntimeError as exc1:
            if "Factor is exactly singular" in str(exc1):
                # mu[max_err_idx] is exactly equal to an
                # eigenvalue, but B[:,max_err_idx] is not
                # a good eigenvector estimate. Therefore
                # perturb mu[max_err_idx] to allow matrix
                # inversion to succeed
                try:
                    invH_shift = splu(
                        H - (1e-15 + mu[max_err_idx]) * speye(n, format="csc")
                    )
                except RuntimeError as exc2:
                    if "Factor is exactly singular" in str(exc2):
                        # If perturbation fails, we'll assume
                        # something bigger is wrong in the solver
                        raise RuntimeError(
                            "Unable to invert matrix H when shifted by "
                            f"{mu[max_err_idx]}*I. Check whether something "
                            "is wrong with the matrix H's scaling."
                        ) from exc2
                    else:
                        raise
            else:
                raise
        Bgrave = invH_shift.solve(B)
        Bhat, _ = qr(Bgrave, mode="economic")
        mu, B_tilde = eigh(Bhat.T @ H @ Bhat)
        B = Bhat @ B_tilde

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
    R[k,k] == 1 for k<=p-1 (corresponding to the null vectors), R[k,k] = 1/sqrt(2) for
    p <= k <= p+q-1 (corresponding to the singular vectors), and R[k,k] ~= 1e-15 for k > p+q.
    However, we often get a spectrum of values.

    Presently, we use a crude thresholding scheme: if both abs(R_U[k,k]) and abs(R_V[k,k]) are
    greater than zero_tol, we add the columns Q_U[:,k] and Q_V[:,k] to our basis. If only one
    is greater than zero_tol, then we check how many more basis vectors we have for U than for
    V (if m > n) or for V than U (if m < n). If the difference in number of basis vectors is
    less than the size of the null space, we add an additional basis vector, otherwise we drop
    it. This scheme is somewhat arbitrary, but has worked better than any replacement that I
    (Doug A.) have tried.

    After determining bases for U and V, we perform a SVD on U.T @ A @ V in order to
    reconstruct appropriate singular triplets as well as the null vectors.

    Args:
        A: Scipy sparse array with real entries
        evecs: Approximate eigenvectors of the augmented matrix [[0, A.T], [A, 0]]
        zero_tol: Tolerance to drop vectors for basis of singular spaces U and V

    Returns:
        U: Dense array of left singular vectors of A
        svals: 1D dense array ofsingular values of A
        V: Dense array of left singular vectors of A
        null_space: If m < n (m > n) then a basis for the (left) null space
            is returned. If m == n, then no null space is returned.

    """
    assert len(A.shape) == 2
    m, n = A.shape
    p = abs(m - n)

    V = evecs[:n, :]
    U = evecs[n:, :]

    U, R_U, _ = qr(U, mode="economic", pivoting=True)
    r_U = np.abs(np.diag(R_U))

    V, R_V, _ = qr(V, mode="economic", pivoting=True)
    r_V = np.abs(np.diag(R_V))

    n_U = 0
    n_V = 0
    for k in range(min(evecs.shape[1], max(m, n))):
        if m > n:
            if r_U[k] > zero_tol:
                if k < n and r_V[k] > zero_tol:
                    n_U += 1
                    n_V += 1
                elif n_U - n_V < p:
                    n_U += 1
        else:
            if r_V[k] > zero_tol:
                if k < m and r_U[k] > zero_tol:
                    n_U += 1
                    n_V += 1

                elif n_V - n_U < p:
                    n_V += 1

    U = U[:, :n_U]
    V = V[:, :n_V]

    # Want to make sure that the number of null vector candidates
    # is the expected size
    # This assertion should no longer be necessary
    p_hat = abs(U.shape[1] - V.shape[1])
    assert p_hat <= p
    U_sub, svals, VT_sub = svd(U.T @ A @ V, full_matrices=True)
    U = U @ U_sub
    V = V @ VT_sub.T

    if m > n:
        null_space = U[:, -p_hat:]
        assert norm(null_space.T @ A) <= np.sqrt(m * n) * 1e-15
        U = U[:, :-p_hat]
    elif m < n:
        null_space = V[:, -p_hat:]
        assert norm(A @ null_space) <= np.sqrt(m * n) * 1e-15
        V = V[:, :-p_hat]

    # Sort singular values in ascending order
    U = U[:, ::-1]
    V = V[:, ::-1]
    svals = svals[::-1]

    if m == n:
        return U, svals, V
    else:
        return U, svals, V, null_space


def svd_rayleigh_ritz(
    A,
    number_singular_values: int = 10,
    max_iter: int = 100,
    tol: float = 1e-12,
    seed: int = None,
    suppress_warning=False,
):
    """
    Computes smallest singular vectors of the sparse m by n matrix A via Rayleigh-Ritz
    iteration on the augmented matrix [[0, A.T], [A, 0]]. Working with this matrix
    avoids the roundoff error inherent working with the Gram matrices A.T @ A or
    A @ A.T, but also results in the (left) null space polluting the singular spectrum
    if m < n (m > n). Therefore, this method is not appropriate for diagnosing
    optimization problems in which m << n, but is appropriate if the (left) null space
    is desired for diagnosing degrees of freedom.

    Args:
        A: Scipy sparse array with real entries
        number_singular_values: number of small singular values and vectors to compute
        max_iter: Maximum number of iterations for Rayleigh-Ritz iteration
        tol: Tolerance used in stopping condition for Rayleigh-Ritz iteration
        seed: Seed for initializing random number generator in Rayleigh-Ritz iteration
        suppress_warning: Suppress the efficiency warning issued when |m - n| > 10

    Returns:
        U: m by number_singular_values dense array of left singular vectors
        svals: 1D dense array of number_singular_values singular values
        V: n by number_singular_values dense array of left singular vectors
        null: Basis for the (left) null space if m < n (m > n)
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
    l = abs(m - n)
    if l > 10 and not suppress_warning:
        if m > n:
            _log.warning(
                f"Matrix A has a left nullspace of dimension at least {l}, which "
                "degrades the efficiency of SVD algorithms based on the augmented "
                "matrix [[0, A.T], [A, 0]]."
            )
        else:
            _log.warning(
                f"Matrix A has a nullspace of dimension at least {l}, which "
                "degrades the efficiency of SVD algorithms based on the augmented "
                "matrix [[0, A.T], [A, 0]]."
            )
    # Can't obtain more singular values than the dimension of the matrix
    number_singular_values = min(number_singular_values, m, n)

    # The augmented matrix has two eigenvectors for every singular triplet, in addition
    # to eigenvectors corresponding to the (left) null space if m<n (m>n). In order to
    # get the required number of singular triplets, we need to exhaust the null space
    # as well as the duplicate singular triplets.
    n_samples = 2 * number_singular_values + l
    A_aug = block_array([[None, A.T], [A, None]], format="csc")

    _, B, converged = _symmetric_rayleigh_ritz_iteration(
        A_aug, n_samples, tol=tol, max_iter=max_iter, seed=seed
    )

    if not converged:
        raise RuntimeError(
            "Rayleigh-Ritz iteration did not converge! Consider increasing the tolerance or "
            "maximum number of iterations."
        )

    try:
        if m == n:
            U, svals, V = _aug_eig_processing(A, B)
        else:
            U, svals, V, null = _aug_eig_processing(A, B)
    except AssertionError as exc:
        raise RuntimeError(
            "Processing of singular vectors failed despite Rayleigh-Ritz iteration "
            "converging. Please let the IDAES dev team know about this failure so that "
            "this processing step can be refined."
        ) from exc

    # Singular values already in ascending order, so just take the number that we want
    U = U[:, :number_singular_values]
    V = V[:, :number_singular_values]
    svals = svals[:number_singular_values]

    if m == n:
        return U, svals, V
    else:
        return U, svals, V, null
