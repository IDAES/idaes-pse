##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

from pyomo.environ import *
from pyomo.common.fileutils import this_file_dir
from pyomo.core.base.external import AMPLExternalFunction
from pyomo.opt import SolverFactory
from idaes.generic_models.properties import iapws95

import unittest
import pytest

import csv
import os

prop_available = iapws95.iapws95_available()

def between(y, x1, x2):
    return 0 > (y-x1)*(y-x2)


@pytest.mark.iapws
class TestIAPWS95(unittest.TestCase):
    def read_data(self, fname, col):
        dfile = os.path.join(this_file_dir(), fname)
        cond = [] # Tuple (T [K],P [Pa], data) pressure in file is MPa
        with open(dfile, 'r') as csvfile:
           dat = csv.reader(csvfile, delimiter='\t', quotechar='"')
           for i in range(7):
               next(dat) # skip header
           for row in dat:
               try:
                   x = float(row[col])
               except:
                   x = row[col]
               cond.append((float(row[0]), float(row[1])*1e6, x))
        return cond

    def make_model(self):
        model = ConcreteModel()
        model.prop_param = iapws95.Iapws95ParameterBlock()
        model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
        return model

    def un_derivs_fd_test(self, f, x1, d=1e-6, tol=0.01):
        assert(isinstance(f, AMPLExternalFunction))
        y, g, h = f.evaluate_fgh(args=(x1,))
        yf, gf, hf = f.evaluate_fgh(args=(x1 + d,))
        yb, gb, hb = f.evaluate_fgh(args=(x1 - d,))
        gfdf = (yf - y)/d
        gfdb = -(yb - y)/d
        hfdf = (gf[0] - g[0])/d
        hfdb = -(gb[0] - g[0])/d
        #print("{:+.5e} {:+.5e} {:+.5e} | {:+.5e} {:+.5e} {:+.5e}".format(
        #      g[0], gfdb, gfdf, h[0], hfdb, hfdf))

        zero_cut = 1e-9 # how close to zero before maybe it is zero?

        # check that the forward and backward FD approximations are close enough
        # that the accuracy is good enough for the test and that the detivative
        # is not ~ zero.  I know this rough but what can you do?
        if abs(g[0]) > zero_cut:
            assert(abs((gfdf - g[0])/g[0]) < tol or
                   abs((gfdb - g[0])/g[0]) < tol or
                   between(g[0], gfdf, gfdb))

        if abs(h[0]) > zero_cut:
            assert(abs((hfdf - h[0])/h[0]) < tol or
                   abs((hfdb - h[0])/h[0]) < tol or
                   between(h[0], hfdf, hfdb))

    def bin_derivs_fd_test(self, f, x1, x2, d0=1e-6, d1=1e-6, tol=0.01):
        assert(isinstance(f, AMPLExternalFunction))

        y, g, h = f.evaluate_fgh(args=(x1, x2))
        yf0, gf0, hf0 = f.evaluate_fgh(args=(x1 + d0, x2))
        yb0, gb0, hb0 = f.evaluate_fgh(args=(x1 - d0, x2))
        yf1, gf1, hf1 = f.evaluate_fgh(args=(x1, x2 + d1))
        yb1, gb1, hb1 = f.evaluate_fgh(args=(x1, x2 - d1))
        gf = [(yf0 - y)/d0, (yf1 - y)/d1]
        gb = [-(yb0 - y)/d0, -(yb1 - y)/d1]
        hf = [(gf0[0] - g[0])/d0, (gf0[1] - g[1])/d0, (gf1[1] - g[1])/d1]
        hb = [-(gb0[0] - g[0])/d0, -(gb0[1] - g[1])/d0, -(gb1[1] - g[1])/d1]

        #print("{:+.5e} {:+.5e} {:+.5e} | "
        #      "{:+.5e} {:+.5e} {:+.5e} | "
        #      "{:+.5e} {:+.5e} {:+.5e} | "
        #      "{:+.5e} {:+.5e} {:+.5e} | "
        #      "{:+.5e} {:+.5e} {:+.5e}".format(
        #        g[0], gb[0], gf[0],
        #        g[1], gb[1], gf[1],
        #        h[0], hb[0], hf[0],
        #        h[1], hb[1], hf[1],
        #        h[2], hb[2], hf[2]))

        zero_cut = 1e-9 # how close to zero before maybe it is zero?
        # check that the forward and backward FD approximations are close enough
        # that the accuracy is good enough for the test and that the detivative
        # is not ~ zero.  I know this rough but what can you do?
        if abs(g[0]) > zero_cut:
            assert(abs((gf[0] - g[0])/g[0]) < tol or
                   abs((gb[0] - g[0])/g[0]) < tol or
                   between(g[0], gf[0], gb[0]))

        if abs(g[1]) > zero_cut and abs((gf[1] - gb[1])/g[1]) < tol:
            assert(abs((gf[1] - g[1])/g[1]) < tol or
                   abs((gb[1] - g[1])/g[1]) < tol or
                   between(g[1], gf[1], gb[1]))

        if abs(h[0]) > zero_cut and abs((hf[0] - hb[0])/h[0]) < tol:
            assert(abs((hf[0] - h[0])/h[0]) < tol or
                   abs((hb[0] - h[0])/h[0]) < tol or
                   between(h[0], hf[0], hb[0]))

        if abs(h[1]) > zero_cut and abs((hf[1] - hb[1])/h[1]) < tol:
            assert(abs((hf[1] - h[1])/h[1]) < tol or
                   abs((hb[1] - h[1])/h[1]) < tol or
                   between(h[1], hf[1], hb[1]))

        if abs(h[2]) > zero_cut and abs((hf[2] - hb[2])/h[2]) < tol:
            assert(abs((hf[2] - h[2])/h[2]) < tol or
                   abs((hb[2] - h[2])/h[2]) < tol or
                   between(h[2], hf[2], hb[2]))

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    #@pytest.mark.skip(reason="temporary to save time")
    def test_derivs_sat_deltal(self):
        model = self.make_model()
        cond = self.read_data("sat_prop_iapws95_nist_webbook.txt", col=2)
        for i, c in enumerate(cond):
            f = model.prop_in.func_delta_sat_l
            print(c[0])
            self.un_derivs_fd_test(f, 647.096/c[0], d=1e-6, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    #@pytest.mark.skip(reason="temporary to save time")
    def test_derivs_sat_deltav(self):
        model = self.make_model()
        cond = self.read_data("sat_prop_iapws95_nist_webbook.txt", col=2)
        for i, c in enumerate(cond):
            f = model.prop_in.func_delta_sat_v
            print(c[0])
            self.un_derivs_fd_test(f, 647.096/c[0], d=1e-6, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_psat(self):
        model = self.make_model()
        cond = self.read_data("sat_prop_iapws95_nist_webbook.txt", col=2)
        for i, c in enumerate(cond):
            f = model.prop_in.func_p_sat
            print(c[0])
            self.un_derivs_fd_test(f, 647.096/c[0], d=1e-6, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_tau_sat(self):
        model = self.make_model()
        cond = self.read_data("sat_prop_iapws95_nist_webbook.txt", col=2)
        for i, c in enumerate(cond):
            f = model.prop_in.func_tau_sat
            print(c[1])
            if c[1] > 22040000.0: break # dont want to do critical point here.
            self.un_derivs_fd_test(f, c[1]/1000.0, d=1e-3, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_dens(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=2)
        phase = self.read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
                f_dens = model.prop_in.func_delta_liq
            else:
                p = "Vap"
                f_dens = model.prop_in.func_delta_vap
            if c[2] < 100: # less than 100 Pa is very low
                continue
            print("{} {}".format(c[0], c[1]))
            self.bin_derivs_fd_test(
                f_dens, c[1]/1000, 647.096/c[0], d0=1e-3, d1=1e-6, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_hxpt(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=2)
        phase = self.read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
                f = model.prop_in.func_hlpt
            else:
                p = "Vap"
                f = model.prop_in.func_hvpt
            print("{} {}".format(c[0], c[1]))
            self.bin_derivs_fd_test(f, c[1]/1000, 647.096/c[0], d0=1e-3, d1=1e-6, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    #@pytest.mark.skip(reason="temporary to save time")
    def test_derivs_sxpt(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=2)
        phase = self.read_data("prop_iapws95_nist_webbook.txt", col=13)
        for i, c in enumerate(cond):
            if phase[i][2] in ["liquid", "supercritical"]:
                p = "Liq"
                f = model.prop_in.func_slpt
            else:
                p = "Vap"
                f = model.prop_in.func_svpt
            print("{} {} {}".format(c[0], c[1], phase[i][2]))
            self.bin_derivs_fd_test(f, c[1]/1000, 647.096/c[0], d0=1e-3, d1=1e-6, tol=0.01)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    #@pytest.mark.skip(reason="temporary to save time")
    def test_derivs_tau(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=5)
        f = model.prop_in.func_tau
        j = 0
        for i, c in enumerate(cond):
            j += 1
            if j > 50: j = 0
            if j == 0:
                print("{} {} {}".format(c[0], c[1], c[2]))
                self.bin_derivs_fd_test(f, c[2], c[1]/1000, d0=1e-4, d1=1e-3, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    #@pytest.mark.skip(reason="temporary to save time")
    def test_derivs_tau(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=6)
        f = model.prop_in.func_tau_sp
        j = 0
        for i, c in enumerate(cond):
            j += 1
            if j > 50: j = 0
            if j == 0:
                print("{} {} {}".format(c[0], c[1], c[2]))
                self.bin_derivs_fd_test(f, c[2], c[1]/1000, d0=1e-4, d1=1e-3, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_vf_1phase(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=5)
        f = model.prop_in.func_vf
        j = 0
        for i, c in enumerate(cond):
            j += 1
            if j > 50: j = 0
            if j == 0:
                print("{} {} {}".format(c[0], c[1], c[2]))
                self.bin_derivs_fd_test(f, c[2], c[1]/1000, d0=1e-4, d1=1e-3, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_vf_2phase(self):
        model = self.make_model()
        hvdat = self.read_data("sat_prop_iapws95_nist_webbook.txt", col=17)
        hldat = self.read_data("sat_prop_iapws95_nist_webbook.txt", col=5)
        f = model.prop_in.func_vf
        j = 0
        for i, c in enumerate(hvdat):
            T = c[0]
            if T > 647: continue
            p = c[1]/1000
            hv = c[2]
            hl = hldat[i][2]
            ht = hl + (hv - hl)*0.75
            print("{} {} {}".format(T, p, ht))
            self.bin_derivs_fd_test(f, ht, p, d0=1e-4, d1=1e-3, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_phi0(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=2)
        for i, c in enumerate(cond):
            delta = c[2]/322.0
            tau = 647.096/c[0]
            if c[1] < 200: # less than 200 Pa is very low
                continue
            print("{} {}".format(c[0], c[1]))
            self.bin_derivs_fd_test(
                model.prop_in.func_phi0, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phi0_delta, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phi0_delta2, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phi0_tau, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phi0_tau2, delta, tau, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_phir(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=2)
        for i, c in enumerate(cond):
            delta = c[2]/322.0
            tau = 647.096/c[0]
            print("{} {}".format(c[0], c[1]))
            if c[1] < 200:
                continue #low pressure < 200 Pa
            self.bin_derivs_fd_test(
                model.prop_in.func_phir, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phir_delta, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phir_tau, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phir_delta2, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phir_tau2, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_phir_delta_tau, delta, tau, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    #@pytest.mark.skip(reason="temporary to save time")
    def test_derivs_pushfg(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=2)
        for i, c in enumerate(cond):
            delta = c[2]/322.0
            tau = 647.096/c[0]
            if c[1] < 200:
                continue #low pressure < 200 Pa
            print("{} {}".format(c[0], c[1]))
            self.bin_derivs_fd_test(
                model.prop_in.func_p, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_u, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_s, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_h, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_f, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_g, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_w, delta, tau, tol=0.001)

    @pytest.mark.slow
    @pytest.mark.skipif(not prop_available, reason="IAPWS not available")
    @pytest.mark.skip(reason="temporary to save time")
    def test_derivs_cp_cv(self):
        model = self.make_model()
        cond = self.read_data("prop_iapws95_nist_webbook.txt", col=2)
        for i, c in enumerate(cond):
            delta = c[2]/322.0
            tau = 647.096/c[0]
            if c[1] < 200:
                continue #low pressure < 200 Pa
            print("{} {}".format(c[0], c[1]))
            self.bin_derivs_fd_test(
                model.prop_in.func_cv, delta, tau, tol=0.001)
            self.bin_derivs_fd_test(
                model.prop_in.func_cp, delta, tau, tol=0.001)
