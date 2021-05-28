###############################################################################
# Copyright
# =========
#
# Institute for the Design of Advanced Energy Systems Process Systems Engineering
# Framework (IDAES PSE Framework) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so. Copyright
# (C) 2018-2019 IDAES - All Rights Reserved
#
###############################################################################


def almpickle(data, vargs):
    import pickle
    import sys
    # remove lambda function from results dictionary
    if (len(vargs) == 1):
        fname = vargs[0]
    elif (len(vargs) == 2):
        fname = vargs[0]
        dopick = vargs[1]
    else:
        fname = 'alamo.pick'
        dopick = True

    if (not isinstance(data, list)):
        if ('f(model)' in data.keys()):
            data['f(model)'] = ()
        else:
            sys.stdout.write('Results dictionary does not contain lambda function')

    if (dopick):
        pickle.dump(data, open(fname, "wb"))


def postpickle(data):
    # Recompile lambda function from model string and labels
    # import sympy
    from sympy.parsing.sympy_parser import parse_expr
    from sympy import symbols, lambdify

    if (isinstance(data['model'], type({}))):
        data['f(model)'] = {}
        for olab in data['zlabels']:
            model_str = data['model'].split('=')[1].replace('^', '**')
            data['f(model)'][olab] = lambdify([symbols(data['xlabels'])], 
                                              parse_expr(model_str), "numpy")
    else:
        model_str = data['model'].split('=')[1].replace('^', '**')
        data['f(model)'] = lambdify([symbols(data['xlabels'])],
                                    parse_expr(model_str), "numpy")

    return data


def almunpickle(fname):
    # unpickle and relambdify
    from idaes.surrogate import alamopy
    # from alamopy import postpickle
    import pickle

    res = pickle.load(open(fname, 'rb'))
    if ('f(model)' in res[0].keys()):
        res = alamopy.postpickle(res)
    return res
