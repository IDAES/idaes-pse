#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for linear algebra utility methods.
"""

import pytest

import numpy as np
from numpy.linalg import norm, qr
from numpy.random import default_rng
from scipy.sparse import csc_array

from pyomo.environ import ConcreteModel, Param, log10, Var, value
from idaes.core.util.linalg import svd_explicit_grammian, svd_rayleigh_ritz, svd_inverse_aug


__author__ = "Douglas Allan"

def _random_svd(
        n_rows:int,
        n_columns:int,
        n_small:int,
        seed:int,
        eps_small: float = 1e-12,
        eps_min: float = 1e-14
    ):
    """
    Constructs an n_rows by n_columns Numpy 2D array A from its singular value decomposition
    A = U @ diag(sigma) @ V.T. U and V are created through orthogonalization of random 
    arrays. The 1D array sigma is generated such that min(n_rows, n_columns) - n_small 
    singluar values are greater than eps_small while n_small singular values are less than
    eps_small.

    Args:
        n_rows: Number of rows in outlet array A
        n_columns: Number of columns in outlet array A
        n_small: Number of singular values smaller than eps_small
        seed: Seed for Numpy's random number generator
        eps_small: Value defining threshold for small singular value (default 1e-12)
        eps_min: Smallest singular value possible
    Returns:
        A: Random matrix generated by function
        U: m by number_singular_values dense array of left singular vectors
        svals: 1D dense array of number_singular_values singular values, sorted from smallest to largest
        V: n by number_singular_values dense array of left singular vectors
    """
    n_svals = min(n_rows, n_columns)
    assert n_svals >= n_small
    rng_obj = default_rng(seed)

    U = rng_obj.standard_normal((n_rows, n_svals))
    V = rng_obj.standard_normal((n_columns, n_svals))
    U, _ = qr(U)
    V, _ = qr(V)

    svals_big = rng_obj.uniform(low=log10(eps_small), high=0, size=(n_svals-n_small,))
    svals_big = 10 ** svals_big
    svals_big.sort()
    svals_small = rng_obj.uniform(low=log10(eps_min), high=log10(eps_small), size=(n_small,))
    svals_small = 10 ** svals_small
    svals_small.sort()
    svals = np.concatenate([svals_small, svals_big])

    A = U @ np.diag(svals) @ V.T

    return A, U, svals, V

def _projection_residual(U, V):
    """
    Returns the Frobenius norm of the component of U that falls outside
    the orthonormal basis V
    """
    return norm(V@(V.T @ U) - U)

def _assert_subspace_containment(U, V, tol=1e-8):
    """
    Asserts that the matrix U is contained in the subspace spanned by
    the orthonormal basis V to an absolute tolerance of tol.
    """
    assert norm(V@(V.T @ U) - U) == pytest.approx(0, rel=0, abs=tol)

def _assert_subspace_orthogonality(U, V, tol=1e-8):
    """
    Assert that the matrix U is orthogonal to the subspace spanned by
    the orthonormal basis V to an absolute tolerance of tol.
    """
    assert norm(V@(V.T @ U)) == pytest.approx(0, rel=0, abs=tol)

def _test_svd_quality(A, Uhat, svals_hat, Vhat, null_hat=None):
    m, n = A.shape
    assert Uhat.T @ A @ Vhat == pytest.approx(np.diag(svals_hat), rel=1e-6, abs=1e-14)
    # U.T @ A = Sigma @ V.T
    assert norm(Uhat.T @ A, axis=1) == pytest.approx(svals_hat, rel=1e-4, abs=1e-14)
    # A @ V = U @ Sigma
    assert norm(A @ Vhat, axis=0) == pytest.approx(svals_hat, rel=1e-4, abs=1e-14)
    if null_hat is not None:
        if m > n:
            assert norm(null_hat.T @ A) == pytest.approx(0, rel=0, abs=1e-15)
        else:
            assert norm(A @ null_hat) == pytest.approx(0, rel=0, abs=1e-15)

    # Check orthogonality
    assert Uhat.T @ Uhat == pytest.approx(np.eye(Uhat.shape[1]), rel=0, abs=1e-14)
    assert Vhat.T @ Vhat == pytest.approx(np.eye(Vhat.shape[1]), rel=0, abs=1e-14)
    if null_hat is not None:
        # Null vectors are less orthogonal because I'm trying to resist doing a QR on them too
        assert null_hat.T @ null_hat == pytest.approx(np.eye(null_hat.shape[1]), rel=0, abs=1e-12)

# class SVDRayleighRitz:
@pytest.mark.unit
def test_not_sparse():
    A = np.eye(2)
    with pytest.raises(
        ValueError,
        match="This method expects a Scipy sparse array-like as an input but was passed "
        "a dense array-like instead. Try using scipy.linalg.svd for a dense SVD method.",
    ):  
        svd_rayleigh_ritz(A)

@pytest.mark.unit
def test_not_2D():
    A = np.zeros((2,3,4))
    with pytest.raises(
        ValueError,
        match="This method expects a 2D Scipy sparse array-like as input, but was passed "
        "a 3D array-like instead.",
    ):  
        svd_rayleigh_ritz(A)

@pytest.mark.unit
def test_not_converged():
    A = np.eye(2)
    with pytest.raises(
        RuntimeError,
        match="Rayleigh-Ritz iteration did not converge! Consider increasing "
        "the tolerance or maximum number of iterations.",
    ):
        svd_rayleigh_ritz(csc_array(A), max_iter=0)

@pytest.mark.unit
def test_square():
    A, U, svals, V = _random_svd(n_rows=25, n_columns=25, n_small=5, seed=8, eps_min=1e-16)
    Utrunc = U[:,:10]
    Vtrunc = V[:, :10]
    svals_trunc = svals[:10]
    Uhat, svals_hat, Vhat = svd_rayleigh_ritz(csc_array(A), seed=8)
    _assert_subspace_containment(Uhat, Utrunc)
    _assert_subspace_containment(Vhat, Vtrunc)
    assert pytest.approx(svals_trunc, rel=1e-4, abs=1e-14) == svals_hat

@pytest.mark.unit
def test_thin():
    A, U, svals, V = _random_svd(n_rows=30, n_columns=25, n_small=5, seed=8, eps_min=1e-16) #, eps_small=1e-10, eps_min=1e-12)
    Utrunc = U[:,:10]
    Vtrunc = V[:, :10]
    svals_trunc = svals[:10]
    Uhat, svals_hat, Vhat, null_hat = svd_rayleigh_ritz(csc_array(A), seed=8)
    # We've got blending between the 1e-16 singular values and the null space
    # _assert_subspace_containment(Uhat, Utrunc)
    # _assert_subspace_containment(Vhat, Vtrunc)
    # _assert_subspace_orthogonality(null_hat, Utrunc)
    _test_svd_quality(A, Uhat, svals_hat, Vhat, null_hat)
    assert pytest.approx(svals_trunc, rel=1e-4, abs=1e-14) == svals_hat

if __name__ == "__main__":
    A, U, svals, V = _random_svd(n_rows=25, n_columns=25, n_small=5, seed=8)
    Usmall = U[:,:10]
    Vsmall = V[:, :10]
    svals_small = svals[:10]
    Uhat, svals_hat, Vhat = svd_rayleigh_ritz(csc_array(A), max_iter=300)
    Ubetter, svals_better, Vbetter = svd_rayleigh_ritz(csc_array(A), max_iter=300)
    print("Inverse Aug")
    print(norm(Usmall @ (Usmall.T @ Uhat) - Uhat))
    print(norm(Vsmall @ (Vsmall.T @ Vhat) - Vhat))
    print(svals_small)
    print(svals_hat)   

    print("Rayleigh-Ritz")
    print(norm(Usmall @ (Usmall.T @ Ubetter) - Ubetter))
    print(norm(Vsmall @ (Vsmall.T @ Vbetter) - Vbetter))
    print(svals_small)
    print(svals_better)
    # Ubetter, svals_better, Vbetter = svd_explicit_grammian_pseudo(csc_array(A), num_iter=100)
    # Ubetter, svals_better, Vbetter = svd_explicit_grammian(csc_array(A), num_iter=100, use_sparse_chol=True)

# @pytest.fixture(scope="module")
# def simple_model():
#     """Build a simple model for testing."""
#     m = ConcreteModel()
#     m.a = Var(initialize=4.0)
#     m.b = Var(initialize=-4.0)
#     m.e = Param(default=1e-4)

#     return m


# @pytest.mark.unit
# def test_smooth_abs_maths():
#     # Test basic smooth_abs functionalliy
#     assert smooth_abs(4, 0) == 4.0
#     assert smooth_abs(-4, 0) == 4.0
#     assert smooth_abs(10.0, 0.0) == 10.0
#     assert smooth_abs(-10.0, 0.0) == 10.0

#     assert smooth_abs(2, 1e-4) == pytest.approx(2.0, abs=1e-4)
#     assert smooth_abs(-2, 1e-4) == pytest.approx(2.0, abs=1e-4)
#     assert smooth_abs(10) == pytest.approx(10.0, abs=1e-4)
#     assert smooth_abs(-10) == pytest.approx(10.0, abs=1e-4)


# @pytest.mark.unit
# def test_smooth_abs_expr(simple_model):
#     # Test that smooth_abs works with Pyomo components
#     assert value(smooth_abs(simple_model.a, 0)) == 4.0
#     assert value(smooth_abs(simple_model.b, 0)) == 4.0

#     assert value(smooth_abs(simple_model.a, simple_model.e)) == pytest.approx(
#         4.0, abs=1e-4
#     )
#     assert value(smooth_abs(simple_model.b, simple_model.e)) == pytest.approx(
#         4.0, abs=1e-4
#     )


# @pytest.mark.unit
# def test_smooth_abs_a_errors():
#     # Test that smooth_abs returns meaningful errors when given invalid arg
#     with pytest.raises(TypeError):
#         smooth_abs("foo")
#     with pytest.raises(TypeError):
#         smooth_abs([1, 2, 3])


# @pytest.mark.unit
# def test_smooth_abs_eps_errors():
#     # Test that smooth_abs returns meaningful errors when given invalid eps
#     with pytest.raises(TypeError):
#         smooth_abs(1.0, "a")
#     with pytest.raises(TypeError):
#         smooth_abs(1.0, [1, 2, 3])


# @pytest.mark.unit
# def test_smooth_minmax_maths():
#     # Test basic smooth_minmax functionality
#     assert smooth_minmax(1, 2, 0, sense="max") == 2
#     assert smooth_minmax(1, 2, 0, sense="min") == 1
#     assert smooth_minmax(5.0, 3, 0.0, sense="max") == 5
#     assert smooth_minmax(5.0, 3, 0.0, sense="min") == 3

#     assert smooth_minmax(2.0, 12.0, 1e-4, "max") == pytest.approx(12.0, abs=1e-4)
#     assert smooth_minmax(2.0, 12.0, 1e-4, "min") == pytest.approx(2.0, abs=1e-4)
#     assert smooth_minmax(32.0, 12.0, sense="max") == pytest.approx(32.0, abs=1e-4)
#     assert smooth_minmax(32.0, 12.0, sense="min") == pytest.approx(12.0, abs=1e-4)


# @pytest.mark.unit
# def test_smooth_minmax_default_sense():
#     # Test that smooth_minmax defaults to maximise
#     assert (
#         smooth_minmax(
#             1,
#             2,
#             0,
#         )
#         == 2
#     )


# @pytest.mark.unit
# def test_smooth_minmax_expr(simple_model):
#     # Test that smooth_minmax works with Pyomo components
#     assert value(smooth_minmax(simple_model.a, simple_model.b, 0, sense="max")) == 4.0
#     assert value(smooth_minmax(simple_model.a, simple_model.b, 0, sense="min")) == -4.0

#     assert value(
#         smooth_minmax(simple_model.a, simple_model.b, sense="max")
#     ) == pytest.approx(4.0, abs=1e-4)
#     assert value(
#         smooth_minmax(simple_model.a, simple_model.b, sense="min")
#     ) == pytest.approx(-4.0, abs=1e-4)

#     assert value(
#         smooth_minmax(simple_model.a, simple_model.b, simple_model.e, sense="max")
#     ) == pytest.approx(4.0, abs=1e-4)
#     assert value(
#         smooth_minmax(simple_model.a, simple_model.b, simple_model.e, sense="min")
#     ) == pytest.approx(-4.0, abs=1e-4)


# @pytest.mark.unit
# def test_smooth_abs_ab_errors():
#     # Test that smooth_abs returns meaningful errors when given invalid args
#     with pytest.raises(TypeError):
#         smooth_abs("foo", 1)
#     with pytest.raises(TypeError):
#         smooth_abs(3, [1, 2, 3])


# @pytest.mark.unit
# def test_smooth_minmax_eps_errors():
#     # Test that smooth_abs returns meaningful errors when given invalid eps
#     with pytest.raises(TypeError):
#         smooth_minmax(1.0, 1.0, "foo")
#     with pytest.raises(TypeError):
#         smooth_minmax(1.0, 1.0, [1, 2, 3])


# @pytest.mark.unit
# def test_smooth_minmax_sense_errors():
#     # Test that smooth_abs returns meaningful errors when given invalid sense
#     with pytest.raises(ValueError):
#         smooth_minmax(1.0, 1.0, sense="foo")
#     with pytest.raises(ValueError):
#         smooth_minmax(1.0, 1.0, sense=1.0)
#     with pytest.raises(ValueError):
#         smooth_minmax(1.0, 1.0, sense=[1.0])


# @pytest.mark.unit
# def test_smooth_max(simple_model):
#     # Test that smooth_max gives correct values
#     assert smooth_max(3.0, 12.0) == pytest.approx(12.0, abs=1e-4)
#     assert value(
#         smooth_max(simple_model.a, simple_model.b, simple_model.e)
#     ) == pytest.approx(4.0, abs=1e-4)


# @pytest.mark.unit
# def test_smooth_min(simple_model):
#     # Test that smooth_min gives correct values
#     assert smooth_min(3.0, 12.0) == pytest.approx(3.0, abs=1e-4)
#     assert value(
#         smooth_min(simple_model.a, simple_model.b, simple_model.e)
#     ) == pytest.approx(-4.0, abs=1e-4)


# @pytest.mark.unit
# def test_smooth_min(simple_model):
#     # Test that smooth_min gives correct values
#     assert safe_sqrt(4) == pytest.approx(2.0, abs=1e-4)
#     assert safe_sqrt(0, eps=1e-6) == pytest.approx(0.0, abs=1e-3)
#     assert safe_sqrt(-4) == pytest.approx(0.0, abs=1e-4)
#     assert value(safe_sqrt(simple_model.a, simple_model.e)) == pytest.approx(
#         2.0, abs=1e-4
#     )


# @pytest.mark.unit
# def test_smooth_min(simple_model):
#     # Test that smooth_min gives correct values
#     assert safe_log(4) == pytest.approx(1.386294, abs=1e-4)
#     assert safe_log(0) < -5
#     assert safe_log(-4) < -5
#     assert value(safe_log(simple_model.a, simple_model.e)) == pytest.approx(
#         1.386294, abs=1e-4
#     )
