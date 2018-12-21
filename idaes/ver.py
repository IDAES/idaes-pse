##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
A standard interface to version information.

If you want to add a version to a class, e.g. a model, then
simply inherit from ``HasVersion`` and initialize it with the
appropriate version information.
"""
__author__ = 'Dan Gunter <dkgunter@lbl.gov>'


class Version(object):
    """This class attempts to be compliant with a subset of
    `PEP 440 <https://www.python.org/dev/peps/pep-0440/>`_.

    Pre- and post- releases, as well as "release epochs" are not
    supported.
    """
    _specifiers = {
        'alpha': 'a',
        'beta': 'b',
        'candidate': 'rc',
        'development': 'dev',
        'final': ''    # this is the default
    }

    def __init__(self, major, minor, micro, releaselevel='final', serial=0):
        if releaselevel not in self._specifiers:
            raise ValueError('Value "{}" for releaselevel not in ({})'
                             .format(releaselevel,
                                     ','.join(
                                         sorted(self._specifiers.keys()))))
        self.major, self.minor, self.micro = major, minor, micro
        self.releaselevel, self.serial = releaselevel, serial

    def __str__(self):
        return '{}.{}.{}{}'.format(
            self.major,
            self.minor,
            self.micro,
            ('' if self.releaselevel == 'final'
             else self._specifiers[self.releaselevel] + str(self.serial)))


class HasVersion(object):
    """Interface for a versioned class.
    """
    def __init__(self, *args):
        self.version = Version(*args)


# Modify this line to update the version number
package_version = Version(1, 0, 0, 'candidate', 1)


__version__ = str(package_version)
