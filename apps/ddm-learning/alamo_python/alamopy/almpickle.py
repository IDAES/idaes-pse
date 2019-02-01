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
def almpickle(data,vargs):
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

    if (type(data) != type(list())):
        if ('f(model)' in data.keys()):
            data['f(model)'] = ()
        else:
            sys.stdout.write('Results dictionary does not contain lambda function')

    if (dopick):
        pickle.dump( data , open( fname, "wb"))


def postpickle(data):
    # Recompile lambda function from model string and labels
    import sympy
    from sympy.parsing.sympy_parser import parse_expr
    from sympy import symbols, lambdify

    if (type(data['model']) == type({})):
        data['f(model)']={}
        for olab in data['zlabels']:
            data['f(model)'][olab]=lambdify([symbols(data['xlabels'])], parse_expr(data['model'].split('=')[1].replace('^','**')), "numpy")
    else:
        data['f(model)']=lambdify([symbols(data['xlabels'])], parse_expr(data['model'].split('=')[1].replace('^','**')), "numpy")

    return data

def almunpickle(fname):
    # unpickle and relambdify
    import alamopy
    from alamopy import postpickle
    import pickle
    res = pickle.load(open(fname,'rb'))
    if ('f(model)' in res[0].keys()):
        res = alamopy.postpickle(res)
    return res
