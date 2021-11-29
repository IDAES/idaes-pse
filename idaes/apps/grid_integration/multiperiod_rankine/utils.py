def assemble_generator_data(generator_name, gen_params):

    '''
    Build a dictionary of generator parameteres given an rts-gmlc generator name and a dataframe of all generator data

    Arguments:
        generator_name: a generator name in RTS-GMLC dataset.
        gen_params: the RTS-GMLC generator data in Pandas DataFrame

    Returns:
        model_data: a dictionary which has this structure
        {data type name: value}.
    '''

    gen_params = gen_params.set_index('GEN UID',inplace = False)
    properties = ['PMin MW', 'PMax MW', 'Min Up Time Hr', 'Min Down Time Hr',\
                  'Ramp Rate MW/Min', 'Start Heat Warm MBTU', 'Fuel Price $/MMBTU',\
                  'HR_avg_0', 'HR_incr_1', 'HR_incr_2', 'HR_incr_3',\
                  'Output_pct_1','Output_pct_2','Output_pct_3']

    # to dict
    model_data = gen_params.loc[generator_name, properties].to_dict()

    model_data['RU'] = model_data['Ramp Rate MW/Min'] * 60
    model_data['RD'] = model_data['RU']
    model_data['SU'] = min(model_data['PMin MW'], model_data['RU'])
    model_data['SD'] = min(model_data['PMin MW'], model_data['RD'])
    model_data['SU Cost'] = model_data['Start Heat Warm MBTU'] * model_data['Fuel Price $/MMBTU']

    model_data['Min Load Cost'] = model_data['HR_avg_0']/1000 * \
                                     model_data['Fuel Price $/MMBTU'] *\
                                     model_data['PMin MW']

    model_data['Power Segments'] = {}
    model_data['Marginal Costs'] = {}

    model_data['Original Marginal Cost Curve'] = {}
    model_data['Original Marginal Cost Curve'][model_data['PMin MW']] = model_data['Min Load Cost']/model_data['PMin MW']

    for l in range(1,4):
        model_data['Power Segments'][l] = model_data['Output_pct_{}'.format(l)] * model_data['PMax MW']
        model_data['Marginal Costs'][l] = model_data['HR_incr_{}'.format(l)]/1000 * model_data['Fuel Price $/MMBTU']
        model_data['Original Marginal Cost Curve'][model_data['Power Segments'][l]] = model_data['Marginal Costs'][l]

    return model_data