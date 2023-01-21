#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
from pyomo.environ import *
from pyomo.opt import SolverFactory
import numpy as np
from scipy import linalg
from copy import copy


class RoundingRegression:
    def __init__(self, X, Y, complexity_penalty_factor):
        """
        A class for creating a Rounding Regression model.

        Returns:
        self function containing several attributes -

            self.X      : Input design matrix
            self.Y      : Input response vector
            self.LAP    : Pyomo models and object to handle OLS solution for active set
            self.B_ols_sum  :  Sum of magnitude of coefficients from OLS solution (including all variables)
            self.regressors_probability : Scaled regressor probabilities to remove dependence on big-M chosen
            self.regressors_sorted      : Sorted list of regressors by their probabilities

        """
        # Input/output matrix
        self.X = X
        self.Y = Y

        # Construct object to handle OLS solution for active set and build Pyomo models
        self.LAP = LinAlgandPyomo(X, Y, complexity_penalty_factor)

        # Find OLS solution (i.e. including all variables) and construct QP relaxation model
        _, B_ols, self.B_ols_sum = self.LAP.evaluate_obj(np.ones(X.shape[1]))
        QP, opt = self.LAP.construct_QP(X, Y, self.B_ols_sum)

        # Get rounding probabilities from relaxed binaries of QP relaxation
        regressors_probability, _, _ = self.LAP.optimize(opt, QP)

        # Scale and sort regressors
        self.regressors_probability = (
            regressors_probability / (abs(max(regressors_probability, key=abs))) * 0.9
        )
        self.regressors_sorted = np.argsort(regressors_probability)[::-1]

    def randomized_rounding(self):
        """
        Round randomly by stepping through each regressor
        and rounding with prbability equal to scaled binary from QP relaxation value
        """

        # Number of iterations of randomized rounding starting from null model
        number_RR = 5

        # Number of refinement steps of each randomized rounding iteration
        number_Refinement = 5
        opt_obj = 1e5
        opt_regressors = np.zeros(len(self.regressors_probability))

        # Build RR model from null
        for n in range(number_RR):
            # Initialize null model
            regressors = np.zeros(len(self.regressors_probability))
            i, j = 0, 0
            step_obj = 1e5
            step = True
            # Step-through regressors until number_Refinment refinement loops reached
            while step:
                select = np.random.choice(
                    [0, 1],
                    p=[
                        1 - self.regressors_probability[self.regressors_sorted[i]],
                        self.regressors_probability[self.regressors_sorted[i]],
                    ],
                )
                if select == 1 and regressors[self.regressors_sorted[i]] != 1:
                    regressors[self.regressors_sorted[i]] = 1
                    obj, coeff, _ = self.LAP.evaluate_obj(regressors)
                    if obj < step_obj:
                        step_obj = copy(obj)
                        step_coeffs = copy(coeff)
                        step_regressors = copy(regressors)
                    else:
                        regressors[self.regressors_sorted[i]] = 0

                if (
                    select == 0
                    and regressors[self.regressors_sorted[i]] != 0
                    and np.count_nonzero(regressors) != 1
                ):
                    regressors[self.regressors_sorted[i]] = 0
                    obj, coeff, _ = self.LAP.evaluate_obj(regressors)
                    if obj < step_obj:
                        step_obj = copy(obj)
                        step_coeffs = copy(coeff)
                        step_regressors = copy(regressors)
                    else:
                        regressors[self.regressors_sorted[i]] = 1
                i += 1
                if i == min(self.X.shape[1], self.X.shape[0]):
                    if np.count_nonzero(regressors) == 0:
                        i = 0
                    else:
                        i = 0
                        j += 1
                        if j == number_Refinement:
                            step = False

            # Keep current model if best found
            if step_obj < opt_obj:
                opt_obj = copy(step_obj)
                opt_coeffs = copy(step_coeffs)
                opt_regressors = copy(step_regressors)
            self.rr_regressors = copy(opt_regressors)
            self.rr_obj = copy(opt_obj)
            self.rr_coeffs = copy(opt_coeffs)

    def deterministic_rounding(self):
        """
        Round deterministically by stepping through each regressor in order of regressor probability
        """

        improvement = False
        objective_list_det2 = []
        opt_obj = 1e5
        regressors = copy(self.rr_regressors)
        step_regressors = copy(self.rr_regressors)
        step_coeffs = copy(self.rr_coeffs)
        step_obj = copy(self.rr_obj)
        step = True

        i = 0
        j = 1
        # Deterministic rounding loop, loop until no improvement
        while step:
            if regressors[self.regressors_sorted[i]] == 0:
                regressors[self.regressors_sorted[i]] = 1
                obj, coeff, _ = self.LAP.evaluate_obj(regressors)
                if obj < step_obj:
                    step_obj = copy(obj)
                    step_coeffs = copy(coeff)
                    step_regressors = copy(regressors)
                    improvement = True
                else:
                    regressors[self.regressors_sorted[i]] = 0

            else:
                regressors[self.regressors_sorted[i]] = 0
                if np.count_nonzero(regressors) != 0:
                    obj, coeff, _ = self.LAP.evaluate_obj(regressors)
                    if obj < step_obj:
                        step_obj = copy(obj)
                        step_coeffs = copy(coeff)
                        step_regressors = copy(regressors)
                        improvement = True
                    else:
                        regressors[self.regressors_sorted[i]] = 1
                else:
                    regressors[self.regressors_sorted[i]] = 1
            i += 1
            if i == self.X.shape[1]:
                if improvement == False:
                    step = False
                else:
                    improvement = False
                    i = 0
                    j += 1

        if step_obj < opt_obj:
            self.opt_obj = copy(step_obj)
            self.opt_coeffs = copy(step_coeffs)
            self.opt_regressors = copy(step_regressors)

    def build_model(self):
        """
        Method to conduct Randomized rounding and Deterministic rounding combo
        """

        self.randomized_rounding()
        self.deterministic_rounding()

        # Format model found and return
        self.opt_regressors = np.nonzero(self.opt_regressors)[0]
        coeffs = np.zeros(self.X.shape[1])
        for (idx, coefficient) in enumerate(self.opt_coeffs):
            coeffs[self.opt_regressors[idx]] = coefficient
        self.opt_coeffs = coeffs
        return self.opt_obj, self.opt_coeffs, self.opt_regressors


class LinAlgandPyomo:
    def __init__(self, x, y, complexity_penalty_factor):
        """
        Initialize linear algebra and matrix math object that uses pyomo models

        Returns:
        self function containing several attributes -

                self.x  : Input design matrix
                self.y  : Input response vector
                self.regressors_old_A : Active set from previous iteration
                self.regressors_old_QR : Active set form previous iteration used for QR
                self.Q : Q matrix from QR decompostion
                self.R : R matrix from R decomposition
                self.A : Current design matrix (with only columns current active set)
                self.b : Input response vector alias
                self.complexity_penalty : Penalty in objecive function for size of active set

        """
        self.x = x
        self.y = y

        self.regressors_old_A = [1 for i in range(self.x.shape[1])]
        self.regressors_old_QR = [1 for i in range(self.x.shape[1])]
        self.Q, self.R = linalg.qr(self.x)

        self.A = copy(x)
        self.b = copy(y)
        # Complexity penality is a fraction of maximum complexity penalty
        self.complexity_penalty = complexity_penalty_factor * np.linalg.norm(
            x.T @ y, ord=np.inf
        )

    def construct_QP(self, x, y, bigM):
        """
        Construct the QP relaxtion of best subset MIQP

        Args:

            x : Input design matrix
            y : Input response vector
            bigM : Maximum value of coefficient

        Returns:
            self.QP : Pyomo optimization ConcreteModel object
            self.opt : Pyomo optimization object

        """
        regressors = [r for r in range(1, self.x.shape[1] + 1)]
        datapoints = [d for d in range(1, self.x.shape[0] + 1)]

        self.QP = ConcreteModel()
        self.QP.Coeff = Var(regressors, domain=Reals)
        self.QP.z = Var(regressors, domain=UnitInterval)
        self.QP.V = Var(datapoints, domain=Reals)

        def ub_rule(model, i):
            return model.Coeff[i] <= float(bigM) * model.z[i]

        def lb_rule(model, i):
            return model.Coeff[i] >= -float(bigM) * model.z[i]

        def obj_rule(model, i):
            return model.V[i] == (
                float(self.y[i - 1])
                - sum(model.Coeff[j] * float(self.x[i - 1][j - 1]) for j in regressors)
            )

        self.QP.UB = Constraint(regressors, rule=ub_rule)
        self.QP.LB = Constraint(regressors, rule=lb_rule)
        self.QP.Vconst = Constraint(datapoints, rule=obj_rule)

        self.M = float(bigM)

        self.QP.complexity_penalty = Param(
            regressors, initialize=self.complexity_penalty, mutable=True
        )
        self.QP.OBJ = Objective(
            expr=sum((self.QP.V[i]) ** 2 for i in datapoints)
            + sum(self.QP.complexity_penalty[i] * self.QP.z[i] for i in regressors)
        )
        self.opt = SolverFactory("ipopt")

        return self.QP, self.opt

    def optimize(self, opt, model):
        """
        Solve QP model and return relaxed binaries as probabilities

        Arguments:
            opt : Pyomo optimization object
            model : Pyomo optimization ConcreteModel object

        Returns:
            regressors : Binary vector indicating regressors which are active
            coefficients : Coefficient vector indicating coefficients for each regressor
            time : The amount of time needed to solve the optimization problem
        """
        self.results_opt = opt.solve(model, tee=False, keepfiles=False)
        self.solve_time = self.results_opt.solver.time
        regressors = []
        coefficients = []
        for i in range(1, len(model.z) + 1):
            regressors.append(value(model.z[i]))
            coefficients.append(value(model.Coeff[i]))

        return (
            np.array(regressors),
            np.array(coefficients),
            self.results_opt.solver.time,
        )

    def updateA_col(self):
        """
        Update the columns of the design matrix A (i.e. the active set)
        """
        h = 0
        for i in range(self.x.shape[1]):
            if self.regressors_old_A[i] == 0 and self.regressors[i] == 1:
                # New variable inserted, inserts corresponding column into A
                self.A = np.insert(self.A.T, h, self.x.T[i], 0)
                self.A = self.A.T
                h = h + 1
            if self.regressors_old_A[i] == 1 and self.regressors[i] == 1:
                h = h + 1
            if (
                self.regressors_old_A[i] == 1 and self.regressors[i] == 0
            ):  # Variable removed, deletes corresponding column from A
                self.A = np.delete(self.A.T, h, 0)
                self.A = self.A.T

    def updateQR(self):
        """
        Update the QR factorization for the new active set
        """
        h = 0
        for i in range(self.x.shape[1]):
            if self.regressors_old_QR[i] == 0 and self.regressors[i] == 1:
                # New variable inserted, inserts corresponding column into A
                self.Q, self.R = linalg.qr_insert(
                    self.Q, self.R, self.x.T[i].T, h, "col"
                )
                h = h + 1
            if self.regressors_old_QR[i] == 1 and self.regressors[i] == 1:
                h = h + 1
            if (
                self.regressors_old_QR[i] == 1 and self.regressors[i] == 0
            ):  # Variable removed, deletes corresponding column from A
                self.Q, self.R = linalg.qr_delete(self.Q, self.R, h, 1, "col")

    def OLS_soln(self):
        """
        Find the OLS solution of the current active set using QR factorization (if overdetermined problem)
        or else numpy's inbuilt lnalg.lstsq routine for underdetermined case
        """
        self.updateA_col()
        if np.linalg.matrix_rank(self.A) == self.A.shape[1]:
            self.updateQR()
            Rp = self.R[
                : np.count_nonzero(self.regressors)
            ]  # Takes the first 'p' rows of R
            nb = np.dot(self.Q.T, self.b)
            c = nb[
                : np.count_nonzero(self.regressors)
            ]  # Takes the first 'p' rows of nb vector
            d = nb[np.count_nonzero(self.regressors) :]
            self.B_ols = linalg.solve_triangular(Rp, c)
            self.SSRols = sum(
                d[i] ** 2 for i in range(np.shape(self.A)[0] - np.shape(self.A)[1])
            )
            self.B_ols_sum = sum(abs(self.B_ols[i]) for i in range(np.shape(self.A)[1]))

            self.regressors_old_A = copy(self.regressors)
            self.regressors_old_QR = copy(self.regressors)
        else:
            self.B_ols, self.SSRols, rank, s = np.linalg.lstsq(self.A, self.b, rcond=-1)
            self.B_ols_sum = sum(abs(self.B_ols[i]) for i in range(self.A.shape[1]))
            if len(self.SSRols) == 0:
                self.SSRols = 0
            else:
                self.SSRols = self.SSRols[0]
            self.regressors_old_A = copy(self.regressors)

    def evaluate_obj(self, regressors):
        """
        Evaluate objective to MIQP using OLS solution to calculate squared error term plus the complexity penalty

        Arguments:
            regressors : Binary vector indicating regressors which are active

        Returns:
            self.obj    : Approximate objective to MIQP
            self.B_ols  : Coefficient vector for OLS coefficients on active set
            self.B_ols_sum : Sum of magnitude of OLS coefficients for active set

        """
        self.regressors = regressors
        self.OLS_soln()
        self.obj = self.SSRols + self.complexity_penalty * np.count_nonzero(regressors)

        return self.obj, self.B_ols, self.B_ols_sum
