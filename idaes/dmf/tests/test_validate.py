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
Test validate module.
"""
# stdlib
import time
# third-party
import pytest
# package-local
from idaes.dmf import validate

docs = {
    'tabulardata': {
        "data": [{"name": "Density Data",
                  "units": "g/cm^3",
                  "values": [1.0053,
                             1.0188,
                             1.0023],
                  "errors": [5e-05,
                             0.0001,
                             5e-05],
                  "error_type": "absolute"}],
        "meta": [{"datatype": "foo",
                  "authors": "Philip Roth",
                  "date": "1969-01-12",
                  "title": "Portnoy's Complaint",
                  "info": "Random House"}]
    }
}


def test_validate():
    validator = validate.JsonSchemaValidator()
    for schema, instance in docs.items():
        validator.validate(instance, schema)

# TODO: Validation errors


@pytest.mark.skip
def test_validate_caching():
    n = 200
    times = {'cached': 0, 'uncached': 0}
    validators = {'cached': validate.JsonSchemaValidator(),
                  'uncached': validate.JsonSchemaValidator(do_not_cache=True)}
    for key in validators:
        v = validators[key]
        i, t0 = 0, time.time()
        while i < n:
            for schema, instance in docs.items():
                v.validate(instance, schema)
            i += 1
        t1 = time.time()
        times[key] = t1 - t0

    # if caching is working at all, cached time will be measurably
    # less than uncached
    assert times['cached'] < times['uncached']
    nschemas = n * len(docs)
    for key in times:
        mps = times[key] * 1000. / nschemas
        print('{}: {:.2f}ms per doc'.format(key, mps))
    spd = (times['uncached'] - times['cached']) / times['uncached'] * 100.
    print('Speedup from caching: {:.1f}%'.format(spd))

