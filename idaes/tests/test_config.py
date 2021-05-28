##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

import os
import pytest
import idaes

class TestIdaesConfigure(object):

    @pytest.mark.unit
    def test_write_config(self):
        idaes.write_config(
            os.path.join(idaes.testing_directory, "config_testing_orig.json"))
        idaes.write_config(
            os.path.join(idaes.testing_directory, "config_testing_default.json"))


    @pytest.mark.unit
    def test_tmp_config_ctx(self):
        assert idaes.cfg["logging"]["disable_existing_loggers"] == False
        with idaes.temporary_config_ctx():
            idaes.cfg["logging"]["disable_existing_loggers"] = True
            assert idaes.cfg["logging"]["disable_existing_loggers"] == True
        assert idaes.cfg["logging"]["disable_existing_loggers"] == False


    @pytest.mark.unit
    def test_read_default_config(self):
        # set some stuff to not the default
        idaes.cfg.logging.disable_existing_loggers = True
        idaes.cfg.use_idaes_solvers = False
        assert idaes.cfg.logging.disable_existing_loggers == True
        assert idaes.cfg.use_idaes_solvers == False
        #read the default
        idaes.read_config(
            os.path.join(idaes.testing_directory, "config_testing_default.json"))
        assert idaes.cfg.ipopt.options.nlp_scaling_method == "gradient-based"
        assert idaes.cfg.use_idaes_solvers == True
        assert idaes.cfg.logging.disable_existing_loggers == False


    @pytest.mark.unit
    def test_change_env(self):
        with idaes.temporary_config_ctx():
            # default config uses idaes solvers
            idaes.cfg.use_idaes_solvers = False
            idaes.reconfig()
            pf = os.environ["PATH"]
        pt = os.environ["PATH"]
        assert pf != pt
