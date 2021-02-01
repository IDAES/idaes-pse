import os
import pytest
import idaes
import idaes.config as ic

class TestIdaesConfigure(object):

    @pytest.mark.unit
    def test_write_config(self):
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
            idaes.cfg["ipopt-idaes"]["options"]["nlp_scaling_method"]
            == "gradient-based")
        assert idaes.cfg["use_idaes_solvers"] == True
        assert idaes.cfg["logging"]["disable_existing_loggers"] == False


    @pytest.mark.unit
    def test_change_env(self):
        # Keep original isn't True, if there is a config file or it was changed
        idaes.cfg.use_idaes_solvers = False
        idaes.reconfig()
        assert not os.environ["PATH"].startswith("idaes.bin_directory")
        # Set back to default
        idaes.cfg.use_idaes_solvers = True
        idaes.reconfig()
        assert not os.environ["PATH"].startswith("idaes.bin_directory")


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
