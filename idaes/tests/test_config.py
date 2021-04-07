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
import idaes.config as ic

class TestIdaesConfigure(object):

    @pytest.mark.unit
    def test_write_config(self):
        idaes._create_testing_dir()
        ic.write_config(
            os.path.join(idaes.testing_directory, "config_testing_orig.json"))
        ic.write_config(
            os.path.join(idaes.testing_directory, "config_testing_default.json"))

    @pytest.mark.unit
    def test_read_default_config(self):
        """For testing we'll assume we have the default config, and at the end
        of the tests, we'll reset the config back how we found it.
        """
        idaes.cfg["logging"]["disable_existing_loggers"] = True
        assert idaes.cfg["logging"]["disable_existing_loggers"] == True
        ic.read_config(
            os.path.join(idaes.testing_directory, "config_testing_default.json"))
        assert (
            idaes.cfg["ipopt"]["options"]["nlp_scaling_method"]
            == "gradient-based")
        assert idaes.cfg["use_idaes_solvers"] == True
        assert idaes.cfg["logging"]["disable_existing_loggers"] == False


    @pytest.mark.unit
    def test_change_env(self):
        # Keep original isn't True, if there is a config file or it was changed
        idaes.cfg.use_idaes_solvers = False
        idaes.reconfig()
        pf = os.environ["PATH"]
        # Set back to default
        idaes.cfg.use_idaes_solvers = True
        idaes.reconfig()
        pt = os.environ["PATH"]
        assert pf != pt

    @pytest.mark.unit
    def test_revert(self):
        """ Test that the original configuration can be read back in and that
        it is read correctly"""
        ic.read_config(
            os.path.join(idaes.testing_directory, "config_testing_orig.json"))
        orig = idaes.cfg.use_idaes_solvers
        idaes.cfg.use_idaes_solvers = not idaes.cfg.use_idaes_solvers
        ic.read_config(
            os.path.join(idaes.testing_directory, "config_testing_orig.json"))
        assert idaes.cfg.use_idaes_solvers == orig
