# imports
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from ies_model_library import *
from iessolution import *
import scipy.stats
import diptest
import matplotlib.ticker as mtick
from linear_regression import read_raw_data

def read_all_files(directory):
    '''
    Loops through a given directory and converts all the files in the directory to solution objects 

    Arguments: 
        directory: string, where the files are located 

    Returns:
        solution_objects: list of IESSolution classes, classes represent each solution
    '''
    # make a list of signal names
    label_list = os.listdir(directory)

    # remove 00time item
    label_list.remove('00_time.csv')
    label_list.remove('.ipynb_checkpoints')

    # remove repeats
    remove_list = []
    for s in label_list:
        if 'ERCOT_0' in s or 'ERCOT_100' in s or 'CAISO_100' in s:
            remove_list.append(s)

    for r in remove_list:
        label_list.remove(r)


    # read the lmp statistics info
    formatted_data = read_raw_data('formatted_raw_data.csv')


    # make a list of ies_model objects
    iesmodels = {'_model0':Case0(), '_model1':Case1(), '_model3':Case3(), '_model4':Case4(), '_model5':Case5(), '_model6':Case6()}

    solution_objects = []

    # define label in solution object 
    sub = '.csv'

    # loop through the label list
    for i in range(len(label_list)):

        # read csv 
        results_csv = pd.read_csv(directory + '/' + label_list[i])

        case_name = label_list[i][:-len(sub)]

        # make solution object 
        soln = IESSolution(csv_file = directory + '/' + label_list[i], case_object = iesmodels[case_name[len(case_name)-7:]])

        # add a gas price to the solution object
        soln.gas_price = formatted_data.at[(label_list[i][:-(len(sub)+7)] + directory[-3:]), 'Natural Gas Price ($/MMBtu)']

        # make a list of soln objects
        solution_objects.append(soln)

    
    # save a signal label and natural gas price for each
    for i in range(len(solution_objects)):

        solution_objects[i].signal_label = label_list[i][:-len(sub)]

    return solution_objects

def separate_cases(solution_objects):
    '''
    Separates a list of solution objects into 6 individual case lists 
    This function assumes file contains cases 0,1,3,4,5,6 and that the files are sorted by lmp_signal_name_caseX

    Arguments:
        solution_objects: list of IESSolution classes, contains all cases

    Returns:
        case0_objects: list of IESSolution classes, only case 0 
        case1_objects: list of IESSolution classes, only case 1  
        case3_objects: list of IESSolution classes, only case 3 
        case4_objects: list of IESSolution classes, only case 4
        case5_objects: list of IESSolution classes, only case 5 
        case6_objects: list of IESSolution classes, only case 6 
    '''

    # initialize empty lists
    case0_objects = []
    case1_objects = []
    case3_objects = []
    case4_objects = []
    case5_objects = []
    case6_objects = []
        
    for i in range(len(solution_objects)):

        # return the model using the end of the string
        length = len(solution_objects[i].signal_label)

        # trim the last 7 characters '_modelX'

        model_num = solution_objects[i].signal_label[length - 7:]
        # print(model_num)

        
        # append lists based on model number 
        if model_num == '_model0':
            case0_objects.append(solution_objects[i])
        elif model_num == '_model1':
            case1_objects.append(solution_objects[i])
        elif model_num == '_model3':
            case3_objects.append(solution_objects[i])
        elif model_num == '_model4':
            case4_objects.append(solution_objects[i])
        elif model_num == '_model5':
            case5_objects.append(solution_objects[i])
        elif model_num == '_model6':
            case6_objects.append(solution_objects[i])
        else:
            raise NotImplementedError()
        
    # check if model 0 and model 1 are empty 
    if not case0_objects and not case1_objects:
        # read the market results 20 file 
        objects2 = read_all_files('../market_results_20')

        for i in range(len(objects2)):

            # return the model using the end of the string
            length = len(objects2[i].signal_label)

            model_num = objects2[i].signal_label[length - 7:]

            # append 0 and 1 lists only
            if model_num == '_model0':
                case0_objects.append(objects2[i])
            elif model_num == '_model1':
                case1_objects.append(objects2[i])


    
    return case0_objects, case1_objects, case3_objects, case4_objects, case5_objects, case6_objects

def plt_signal_comparison(directory, PROFIT = True, CAPACITY = True, OUTPUT = True, sort = 'median', relative_case = None, save = False):
    '''
    Reads all the result .csv files in given directory and plots each's profit, capacity factor and/or plant output. 
    Saves files in the directory

    Arguments:
        directory: string, directory where the csv files are 
        PROFIT: boolean, plots profit 
        CAPACITY: boolean, plots capacity factor
        OUTPUT: boolean, plots plant output
        sort: string, what to sort the lmp signals by 
            options: median (default): sorted by lmp signal median value
            bimodal: sorted by value of bimodality coefficient 
            gasprice: sorted by natural gas price of the signal
            source: sorted by signal source, chronologically
        relative_case: integer, case number to reference other values to in plot 
            default: None, no relative case, all numbers are plotted in their original scale
        save: boolean, whether or not to save copies of the plots 
            default = False, do not save

    Returns:
        none 
    '''

    solution_objects = read_all_files(directory)

    # plotting data 
    plot_data = {0 : {}, 1 : {}, 3 : {}, 4 : {}, 5: {}, 6: {}}

    # generate keys based on options selected 
    for key in plot_data.keys():
        # case_objects category 
        plot_data[key]['case objects'] = []

        # profit
        if PROFIT:
            plot_data[key]['profit'] = []

        # capacity factor
        if CAPACITY:
            plot_data[key]['p capacity'] = []
            plot_data[key]['h capacity'] = []

        # output
        if OUTPUT:
            plot_data[key]['p output'] = []
            plot_data[key]['h output'] = []


    # sort the objects chosen sorting method 
    if sort == 'median':
        sort_func = lambda soln_obj: np.median(soln_obj.lmp)
    elif sort == 'bimodal':
        sort_func = lambda soln_obj: ((scipy.stats.skew(soln_obj.lmp)**2)+1)/(scipy.stats.kurtosis(soln_obj.lmp) + ((3*(len(soln_obj.lmp)-1)**2)/((len(soln_obj.lmp)-2)*(len(soln_obj.lmp)-3))))
    elif sort == 'gasprice':
        sort_func = lambda soln_obj: soln_obj.gas_price
    elif sort == 'source':
        sort_func = None
    
    else:
        raise NotImplementedError('Please choose to sort lmp signals by median of bimodality coefficient by setting sort to median or bimodal')

    # use python sort function
    if sort != 'source':
        solution_objects.sort(key = sort_func)
    else:
        # perform sort by source 

        # separate list into years and different projections 
        list_2019 = []
        list_2022 = []
        list_nrel = []
        list_princeton = []
        list_netl = []

        for obj in solution_objects:
            # print(obj.signal_label[-11:-8])
            if obj.signal_label[-11:-7] == '2019':
                list_2019.append(obj)
            elif obj.signal_label[-11:-7] == '2022':
                list_2022.append(obj)

            elif obj.signal_label[-11:-7] == '2030':
                list_princeton.append(obj)
            elif obj.signal_label[-11:-7] == '2035':
                # check nrel 
                if obj.signal_label[0:4] == 'MiNg':
                    list_nrel.append(obj)
            else:
                list_netl.append(obj)

        # sort each of these lists by median lmp 
        sort_func = lambda soln_obj: np.median(soln_obj.lmp)

        list_2019.sort(key = sort_func)
        list_2022.sort(key = sort_func)
        list_nrel.sort(key = sort_func)
        list_princeton.sort(key = sort_func)
        list_netl.sort(key = sort_func)

        # append them all together into the new solution objects list
        solution_objects = list_2019 + list_2022 + list_princeton + list_netl + list_nrel

    # separate into cases
    plot_data[0]['case objects'], plot_data[1]['case objects'], plot_data[3]['case objects'], \
        plot_data[4]['case objects'], plot_data[5]['case objects'], plot_data[6]['case objects'] = separate_cases(solution_objects)

    if sort != 'source':
        # use python sort function to sort each of the case objects lists 
            for d in plot_data.keys():

                plot_data[d]['case objects'].sort(key = sort_func)

    else:
        for d in plot_data.keys():
            
            # perform sort by source 

            # separate list into years and different projections 
            list_2019 = []
            list_2022 = []
            list_nrel = []
            list_princeton = []
            list_netl = []

            for obj in plot_data[d]['case objects']:
                if obj.signal_label[-11:-7] == '2019':
                    list_2019.append(obj)
                elif obj.signal_label[-11:-7] == '2022':
                        list_2022.append(obj)

                elif obj.signal_label[-11:-7] == '2030':
                    list_princeton.append(obj)
                elif obj.signal_label[-11:-7] == '2035':
                    # check nrel 
                    if obj.signal_label[0:4] == 'MiNg':
                        list_nrel.append(obj)
                else:
                    list_netl.append(obj)

            # sort each of these lists by median lmp 
            sort_func = lambda soln_obj: np.median(soln_obj.lmp)

            list_2019.sort(key = sort_func)
            list_2022.sort(key = sort_func)
            list_nrel.sort(key = sort_func)
            list_princeton.sort(key = sort_func)
            list_netl.sort(key = sort_func)

            # append them all together into the new solution objects list
            plot_data[d]['case objects'] = list_2019 + list_2022 + list_princeton + list_netl + list_nrel

    # save the profits, capacity factors and outputs
    for key in plot_data.keys():
        for i in range(len(plot_data[key]['case objects'])):

            if PROFIT:
                if relative_case is not None:
                    # profit list
                    plot_data[key]['profit'].append(((plot_data[relative_case]['case objects'][i].profit - plot_data[key]['case objects'][i].profit)/plot_data[relative_case]['case objects'][i].profit)*100)
                else:
                    # profit_list 
                    plot_data[key]['profit'].append(plot_data[key]['case objects'][i].profit)

            if CAPACITY:
                # capacity_lists 
                if key == 6:
                    plot_data[key]['p capacity'].append(0)
                else:
                    plot_data[key]['p capacity'].append(sum(plot_data[key]['case objects'][i].p_dispatch)/(len(plot_data[key]['case objects'][i].horizon)*plot_data[key]['case objects'][i].plant_nameplate))

                if key == 0 or key == 1:
                    plot_data[key]['h capacity'].append(0)
                else:
                    plot_data[key]['h capacity'].append(sum(plot_data[key]['case objects'][i].h_dispatch)/(len(plot_data[key]['case objects'][i].horizon)*plot_data[key]['case objects'][i].h_nameplate))

            if OUTPUT:
                # output lists 
                plot_data[key]['p output'].append(sum(plot_data[key]['case objects'][i].p_dispatch))

                if key == 0 or key == 1:
                    plot_data[key]['h output'].append(0)
                else:
                    plot_data[key]['h output'].append(sum(plot_data[key]['case objects'][i].h_dispatch)*3600) 



    # in one list, remove the final 7 characters to create signal only labels (excluding _modelX)
    for i in range(len(plot_data[0]['case objects'])):
        # remove the end of the strings 
        plot_data[0]['case objects'][i].signal_label = plot_data[0]['case objects'][i].signal_label[:-7]

    final_labels = [plot_data[0]['case objects'][i].signal_label for i in range(len(plot_data[0]['case objects']))]

    # define formatting strings for each case
    fmt_str = {0: '.--', 1: 'o--', 3: 's-', 4: '^-', 5: 'P-', 6: 'D:'}
    colors = {0: 'r', 1 : 'b', 3 : 'g', 4 : 'm', 5 : 'darkorange', 6 : 'deeppink'}
    labels = {0: 'NGCC', 1 : 'SOFC', 3 : 'NGCC + SOEC', 4 : 'rSOC', 5 : 'SOFC + SOEC', 6 : 'SOEC'}
    if relative_case is not None:
        y_labels = {'profit': 'Optimal Annual Profit (M$/yr)',\
            'relative profit': 'Optimal Percent Change in Profit \nrelative to {0}'.format(labels[relative_case]),\
            'p capacity':'Annual Power Capacity Factor\n(MW Produced/MW Capacity)' ,\
            'h capacity': 'Annual Hydrogen Capacity Factor\n(kg Produced/kg Capacity)', \
            'p output': 'Annual Power Output (MW/yr)', \
            'h output': 'Annual Hydrogen Output (kg/yr)'}
    else:
        y_labels = {'profit': 'Optimal Annual Profit (M$/yr)',\
            'p capacity':'Annual Power Capacity Factor\n(MW Produced/MW Capacity)' ,\
            'h capacity': 'Annual Hydrogen Capacity Factor\n(kg Produced/kg Capacity)', \
            'p output': 'Annual Power Output (MW/yr)', \
            'h output': 'Annual Hydrogen Output (kg/yr)'}

    for key in plot_data[0]:

        # skip case objects
        if key == 'case objects':
            continue
        
        # make profit plot 
        fig, ax = plt.subplots(figsize = (20, 9))

        # create color map 
        cmap = plt.cm.get_cmap('gist_rainbow')

        if sort == 'median':
            # normalize it to be between the min and max median LMP 
            vmin = min([sort_func(c) for c in plot_data[0]['case objects']])
            vmax = max([sort_func(c) for c in plot_data[0]['case objects']])

            norm = Normalize(vmin, vmax)

            # set background color sections based on median lmp
            for c in range(len(plot_data[0]['case objects'])):
                ax.axvspan(c-0.4999, c+0.4999, ymin = 0.9, ymax = 0.99, alpha = 1.0, color = cmap(norm(sort_func(plot_data[0]['case objects'][c]))))
            

            # surround the outer box with a black border
            plt.axvline(x = 0-0.5, ymin = 0.9, ymax = 0.99, color = 'k', linewidth = 4)
            plt.axvline(x = len(final_labels)-1 + 0.5, ymin = 0.9, ymax = 0.99, color = 'k', linewidth = 4)


        # plot 0 line
        plt.axhline(y=0, color = 'black', linewidth = 3, alpha = 0.5)


        # plot data
        # print(final_labels)
        for c in plot_data.keys():
            ax.plot(range(len(final_labels)), plot_data[c][key], fmt_str[c], color = colors[c], label = labels[c], markersize = 8, linewidth = 3, zorder = 3)

        # add horizontal lines on colorbar    
        axis_len = ax.get_ylim()[1] - ax.get_ylim()[0]

        if sort == 'median':
            plt.axhline(y=ax.get_ylim()[1]-((1-0.9)*axis_len),xmin = 0.0095, xmax = 0.9901,  color = 'black', linewidth = 4, zorder = 1)
            plt.axhline(y=ax.get_ylim()[1]-((1-0.99)*axis_len),xmin = 0.0095, xmax = 0.9901,  color = 'black', linewidth = 4, zorder = 1)


            sortfuncvals = []
            multiples_of_10 = {10: 1000, 20:1000, 30:1000, 40:1000, 50:1000, 65: 1000}

            # label multiples of 10 on color bar 
            for c in range(len(plot_data[0]['case objects'])):
                # make a list of the values that sort func returns 
                sortfuncvals.append(sort_func(plot_data[0]['case objects'][c]))


            for k in multiples_of_10.keys():
                difference = 1000
                for j in range(len(sortfuncvals)):
                # find the closest median LMP to each multiple of 10 from 10 - 70  and save it in the dict 

                    if abs(sortfuncvals[j] - k) <= difference: 
                        multiples_of_10[k] = j
                        difference = abs(sortfuncvals[j] - k)

            for k in multiples_of_10.keys():
                # print the number over the bar at that value
                plt.text(x = multiples_of_10[k], y = ax.get_ylim()[1] - ((1-0.94)*axis_len), s = ('$' + str(k) + '/MWh'), color = 'k', fontweight = 'bold', fontsize = 14, zorder = 3)


        # remove margins
        plt.margins(x=0.01)

        # add hydrogen price label 
        plt.title('Hydrogen Price = ${0}/kg'.format(solution_objects[2].h2_price[0]), pad = 47, fontsize = 30, fontweight = 'bold')

        # set axes labels 
        if sort == 'median':
            plt.xlabel('Price Scenarios (sorted from lowest to highest median LMP)', fontsize = 30, fontweight = 'bold')
        elif sort == 'bimodal':
            plt.xlabel('Price Scenarios (sorted from lowest to highest bimodality)', fontsize = 30, fontweight = 'bold')
        elif sort == 'gasprice':
            plt.xlabel('Price Scenarios (sorted from lowest to highest natural gas price')
        elif sort == 'source':
            plt.xlabel('Price Scenarios (sorted chronologically and by projection source)', fontsize = 30, fontweight = 'bold')
        if relative_case is not None and key == 'profit':
            plt.ylabel(y_labels['relative profit'], fontsize = 30, fontweight = 'bold')
        else:
            plt.ylabel(y_labels[key], fontsize = 30, fontweight = 'bold')

        # major and minor ticks
        plt.tick_params(direction='in', top = True, right = True)
        plt.minorticks_on()
        plt.tick_params(which = 'minor', direction = 'in', top = False, bottom = False, right = True)


        # set y limits
        if relative_case is not None and key == 'profit':
            plt.ylim(ymin = -250, ymax = 250)

        # set labels
        plt.xticks(range(len(final_labels)), final_labels, size = 12,  fontweight = 'bold', rotation = 90)

        # add tick labels every 5 lines
        # calculate the length of the y axis 
        for i, label in enumerate(final_labels):
            if i % 5 == 0:
                if i == 0:
                    # label 1 not 0
                    plt.text(i, ax.get_ylim()[0]+(0.01*axis_len), i+1, ha = 'center', fontsize = 12, fontweight = 'bold')
                else:
                    plt.text(i-1, ax.get_ylim()[0]+(0.01*axis_len), i, ha='center', fontsize=12, fontweight = 'bold')

        plt.yticks(size = 20, fontweight = 'bold')

        # legend
        legend_properties = {'weight':'bold', 'size':20}
        plt.legend(loc = 'upper center', ncol = 6, bbox_to_anchor = (0.5, 1.1), prop = legend_properties)
        plt.grid(zorder = 2)

        

        # generate plot name 
        h2_price_string = str(solution_objects[3].h2_price[0])
        if relative_case is not None and key == 'profit':
            plot_name = 'all_{0}_{1}{2}_{3}_relativeto{4}'.format(key, h2_price_string.split('.')[0], h2_price_string.split('.')[1], sort, labels[relative_case])
        else:
            plot_name = 'all_{0}_{1}{2}_{3}'.format(key, h2_price_string.split('.')[0], h2_price_string.split('.')[1], sort)
        # print(plot_name)

        if save:
            plt.savefig((plot_name + '.png'), dpi = 300, bbox_inches = 'tight')
            plt.savefig((plot_name + '.pdf'), dpi = 300, bbox_inches = 'tight')

        plt.show()


    return

def plot_profit_scatter(fdest_data = 'formatted_raw_data.csv', x = 'log10p', y = 'medianlmp', h2_price = 20, side_by_side = False, save = False):
    '''
    Generates a scatter plot for a single case where marker size is determined by NPV.

    Arguments:
        fdest_data: string, destination of data file 
            default = 'formatted_raw_data.csv' data file in this directory
        x: string, what to plot on the x axes
            options = log10p, gasprice
        y: string, what to plot on the y axes
            options = medianlmp, gasprice
        h2_price: int, hydrogen price of data to be plotted. options = 0 (all prices), 10 ($1/kg), 15 ($1.5/kg), 20 ($2/kg), 25 ($2.5/kg), 30 ($3/kg)
            default = 20, $2/kg hydrogen data plotted
        side_by_side: boolean, if True, creates two subplots side by side that separate the negative profit and positive profit bubbles
            defualt = False, produces a single plot per case with all points regardless of NPV
        save: boolean, whether or not to save the plots
            default = False, do not save

    Returns:
        None
    '''

    # read raw data 
    formatted_data = read_raw_data(fdest_data)

    # establish labels and plot colors 
    colors = {'NGCC': 'r', 'SOFC' : 'b', 'NGCC+SOEC' : 'g', 'rSOC' : 'm', 'SOFC+SOEC' : 'darkorange', 'SOEC' : 'deeppink'}
    labels = {'NGCC': 'NGCC', 'SOFC' : 'SOFC', 'NGCC+SOEC' : 'NGCC + SOEC', 'rSOC' : 'rSOC', 'SOFC+SOEC' : 'SOFC + SOEC', 'SOEC' : 'SOEC'}
    xy_labels = {'medianlmp': '50th Percentile ($/MWh)', 'log10p': 'Dip Test p-value', 'gasprice': 'Natural Gas Price ($/MMBtu)'}

    for c in colors.keys():

        # make lists of data 
        xlist = []
        ylist = []
        profit = []

        for s in formatted_data.index:
            if h2_price == 0:
                xlist.append(formatted_data.at[s, xy_labels[x]])
                ylist.append(formatted_data.at[s, xy_labels[y]])
                profit.append(formatted_data.at[s, (c + ' profit (M$)')])
            else:
                if s[len(s)-2:] == str(h2_price):
                    xlist.append(formatted_data.at[s, xy_labels[x]])
                    ylist.append(formatted_data.at[s, xy_labels[y]])
                    profit.append(formatted_data.at[s, (c + ' profit (M$)')])

        # create two subplots with shared y axis 
        fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (16,8), sharey = True)

        # make a scatter plot
        for i in range(len(profit)):
            # if points are positive make them the corresponding technology color, if not make them grey 
            if profit[i] >= 0:
                axs[1].scatter(xlist[i], ylist[i], facecolor = colors[c], edgecolor = 'k', s = abs(profit[i])*5, zorder = 2)
            else:
                axs[0].scatter(xlist[i], ylist[i], edgecolor = 'k', facecolor = 'gray', s = abs(profit[i])*5, zorder = 2)

        # titles and axes labels 
        if h2_price == 0:
            fig.suptitle('{0} profits (M$)'.format(labels[c]), fontsize = 30, fontweight = 'bold')
        else:
            fig.suptitle(r'{0} profits (M\$) - \${1}.{2}/kg H$_2$'.format(labels[c], str(h2_price)[0], str(h2_price)[1]), \
                        fontsize = 30, fontweight = 'bold')
        axs[0].set_title('Negative Profit', fontsize = 20, fontweight = 'bold')
        axs[1].set_title('Positive Profit', fontsize = 20, fontweight = 'bold')
        fig.supxlabel(xy_labels[x], fontsize = 25, fontweight = 'bold')
        fig.supylabel(xy_labels[y], fontsize = 25, fontweight = 'bold', x = 0.05)

        # increase size of axes ticks 
        for label in axs[0].get_xticklabels() + axs[0].get_yticklabels() + axs[1].get_xticklabels() + axs[1].get_yticklabels():
            label.set_fontsize(15)
            label.set_fontweight('bold')

        gray_patch = mpatches.Patch(facecolor = 'gray', edgecolor = 'k', label = 'Negative\nProfit')
        color_patch = mpatches.Patch(facecolor = colors[c], edgecolor = 'k', label = 'Positive\nProfit')

         # legend
        legend_properties = {'weight':'bold', 'size':30}
        
        # individual legends for each plot
        axs[0].legend(handles = [gray_patch], loc = 'lower right', prop = legend_properties)
        axs[1].legend(handles = [color_patch], loc = 'lower right', prop = legend_properties)

        # grid 
        axs[0].grid()
        axs[1].grid()

        # save figure in this directory if option is selected 
        if save:
            plt.savefig('{0}_{1}_{2}_{3}.png'.format(c, x, y, str(h2_price)), dpi = 300, bbox_inches = 'tight')
            plt.savefig('{0}_{1}_{2}_{3}.pdf'.format(c, x, y, str(h2_price)), dpi = 300, bbox_inches = 'tight')


        plt.show()


    return 

def plot_positiveprofit_barplot(fdest_data = 'formatted_raw_data.csv', h2_price = 0, save = False):
    '''
    Plots a bar chart depicting the percentage of LMP scenarios that result in positive profit. Displays plots and optionally saves in this directory.

    Arguments:
        fdest_data: string, destination of raw data file 
            default = 'formatted_raw_data.csv', LMP statistics and optimization results from analysis 
        h2_price: int, hydrogen price of scenarios. options = 0 (all prices), 10 ($1/kg), 15 ($1.5/kg), 20 ($2/kg), 25 ($2.5/kg), 30 ($3/kg)
            default = 0, uses all 5 hydrogen price scenarios. 
        save: boolean, whether or not to save the plot as a .png and .pdf in this directory
            default = False, does not save

    Returns:
        None 
    '''

    # make a dictionary of all the cases 
    winners = {'NGCC': 0, 'SOFC': 0, 'NGCC+SOEC': 0, 'rSOC': 0, 'SOFC+SOEC': 0, 'SOEC':0}

    # read data from formatted data file 
    formatted_data = pd.read_csv('formatted_raw_data.csv')
    formatted_data.set_index('Unnamed: 0', inplace = True)

    # loop through every line in the formatted data frame 
    for s in formatted_data.index:

        # plot results across all the hydrogen prices 
        if h2_price == 0:
            # save all profit values 
            profits = {'NGCC': formatted_data.at[s,'NGCC profit (M$)'], 'SOFC': formatted_data.at[s,'SOFC profit (M$)'], \
                       'NGCC+SOEC': formatted_data.at[s,'NGCC+SOEC profit (M$)'], 'rSOC': formatted_data.at[s,'rSOC profit (M$)'], \
                        'SOFC+SOEC': formatted_data.at[s,'SOFC+SOEC profit (M$)'], 'SOEC': formatted_data.at[s,'SOEC profit (M$)']}

            # increment counter if the profit value is positive 
            for key in winners.keys():
                if profits[key] >= 0:
                    winners[key] += 1

        # plot results for a single hydrogen price        
        else:
            if s[len(s)-2:] == str(h2_price):
                 # save all profit values 
                profits = {'NGCC': formatted_data.at[s,'NGCC profit (M$)'], 'SOFC': formatted_data.at[s,'SOFC profit (M$)'], \
                        'NGCC+SOEC': formatted_data.at[s,'NGCC+SOEC profit (M$)'], 'rSOC': formatted_data.at[s,'rSOC profit (M$)'], \
                            'SOFC+SOEC': formatted_data.at[s,'SOFC+SOEC profit (M$)'], 'SOEC': formatted_data.at[s,'SOEC profit (M$)']}

                # increment counter if the profit value is positive 
                for key in winners.keys():
                    if profits[key] >= 0:
                        winners[key] += 1

    colors = {'NGCC': 'r', 'SOFC' : 'b', 'NGCC+SOEC' : 'g', 'rSOC' : 'm', 'SOFC+SOEC' : 'darkorange', 'SOEC' : 'deeppink'}
    labels = {'NGCC': 'NGCC', 'SOFC' : 'SOFC', 'NGCC+SOEC' : 'NGCC + SOEC', 'rSOC' : 'rSOC', 'SOFC+SOEC' : 'SOFC + SOEC', 'SOEC' : 'SOEC'}

    fig, ax = plt.subplots(figsize = (12,6))

    # loop through the items in the dictionary 
    for c in colors.keys():
        plt.barh(labels[c], winners[c], color = colors[c], edgecolor = 'k', label = labels[c], zorder = 2)

    if h2_price == 0:
        plt.xlabel('Scenarios with Positive Profits', fontsize = 25, fontweight = 'bold')
    else:
        plt.xlabel(r'Scenarios with Positive Profits - \${0}.{1}/kg H$_2$'.format(str(h2_price)[0], str(h2_price)[1]), \
                   fontsize = 25, fontweight = 'bold')

    # set the xlabels as the names not the numbers
    plt.yticks(ha = 'right', fontsize = 15, fontweight = 'bold')
    plt.xticks(fontsize = 15, fontweight = 'bold')

    # format the x axis as percentages
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(61))

    # place grid behind the bars
    plt.grid(zorder = 1)

    # save if option was selected
    if save:
        plt.savefig('positiveprofits_bar_{0}.png'.format(h2_price), dpi = 300, bbox_inches = 'tight')
        plt.savefig('positiveprofits_bar_{0}.pdf'.format(h2_price), dpi = 300, bbox_inches = 'tight')

    # display plot
    plt.show()

    return

def plot_bimodalprofit_scatter(fdest_data = 'formatted_raw_data.csv', technology = 'SOFC', h2_price = 0, save = False):
    '''
    Plot a scatter plot of annual profit results comparing integrated systems with standalone systems. 
    x axis represents bimodality coefficient and y axis represents dip test statistic.
    Bubble size corresponds to absolute value of profit difference, with color representing positive 
    difference and grey representing negative difference.

    Arguments:
        fdest_data: string, destination of raw data file
            default = 'formatted_raw_data.csv', csv file saved in this directory 
        technology: string, technology type to run comparison, options = 'SOFC' and 'NGCC'
            default = 'SOFC', compares SOFC+SOEC with standalone SOFC and SOEC
        h2_price: int, hydrogen price of scenarios. options = 0 (all prices), 10 ($1/kg), 15 ($1.5/kg), 20 ($2/kg), 25 ($2.5/kg), 30 ($3/kg)
            default = 0, uses all 5 hydrogen price scenarios. 
        save: boolean, whether or not to save the plot as a .png and .pdf in this directory
            default = False, does not save

    Returns:
        None
    '''

    # read data from formatted data file 
    formatted_data = pd.read_csv('formatted_raw_data.csv')
    formatted_data.set_index('Unnamed: 0', inplace = True)

    # create lists for the data 
    BC = []
    diptest = []
    profitpow = []
    profith2 = []

    # define the systems to be compared based on technology 
    if technology == 'SOFC':
        integratedsystem = 'SOFC+SOEC profit (M$)'
        powersystem = 'SOFC profit (M$)'
        h2system = 'SOEC profit (M$)'
    elif technology == 'NGCC':
        integratedsystem = 'NGCC+SOEC profit (M$)'
        powersystem = 'NGCC profit (M$)'
        h2system = 'SOEC profit (M$)'

    # loop through the df and append the lists with all data points or only data points with selected hydrogen price 
    for i in formatted_data.index:
        if h2_price == 0:
            BC.append(formatted_data.at[i, 'Bimodality Coefficient'])
            diptest.append(formatted_data.at[i, 'Dip Test Statistic'])
            profitpow.append(formatted_data.at[i, integratedsystem] - formatted_data.at[i, powersystem])
            profith2.append(formatted_data.at[i, integratedsystem] - formatted_data.at[i, h2system])
        else:
            if i[len(i)-2:] == str(h2_price):
                BC.append(formatted_data.at[i, 'Bimodality Coefficient'])
                diptest.append(formatted_data.at[i, 'Dip Test Statistic'])
                profitpow.append(formatted_data.at[i, integratedsystem] - formatted_data.at[i, powersystem])
                profith2.append(formatted_data.at[i, integratedsystem] - formatted_data.at[i, h2system])

    # create two subplots with shared y axis 
    fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (16,8), sharey = True)

    # define technology colors 
    colors = {'SOFC': 'darkorange', 'NGCC': 'green'}

    # make a scatter plot
    for i in range(len(profitpow)):
        # if points are positive make them the corresponding technology color, if not make them grey 
        if profitpow[i] >= 0:
            axs[0].scatter(BC[i], diptest[i], facecolor = colors[technology], edgecolor = 'k', s = abs(profitpow[i])*5, zorder = 2)
        else:
            axs[0].scatter(BC[i], diptest[i], edgecolor = 'k', facecolor = 'gray', s = abs(profitpow[i])*5, zorder = 2)
            

        if profith2[i]>= 0:
            axs[1].scatter(BC[i], diptest[i], facecolor = colors[technology],edgecolor = 'k', s = abs(profith2[i])*5, zorder = 2)
        else:
            axs[1].scatter(BC[i], diptest[i], edgecolor = 'k', facecolor = 'gray', s = abs(profith2[i])*5, zorder = 2)


    # convert y axis to log scale between 0 and 1 (format at scalars)
    plt.yscale('log')
    plt.ylim(ymin = 0.0, ymax = 1.0)
    axs[0].yaxis.set_major_formatter(mtick.ScalarFormatter())

    # set x axis limits from 0 to 1
    axs[0].set_xlim(xmin = 0.0, xmax = 1.0)
    axs[1].set_xlim(xmin = 0.0, xmax = 1.0)


    # titles and axes labels 
    if h2_price == 0:
        fig.suptitle('Benefits of Integrated System ({0} + SOEC)'.format(technology), fontsize = 30, fontweight = 'bold')
    else:
        fig.suptitle(r'Benefits of Integrated System ({0} + SOEC) - \${1}.{2}\kg H$_2$'.format(technology, str(h2_price)[0], str(h2_price)[1]), \
                     fontsize = 30, fontweight = 'bold')
    axs[0].set_title('vs. Standalone {}'.format(technology), fontsize = 20, fontweight = 'bold')
    axs[1].set_title('vs. Standalone SOEC', fontsize = 20, fontweight = 'bold')
    fig.supxlabel('Bimodality Coefficient', fontsize = 25, fontweight = 'bold')
    fig.supylabel('Dip Test Statistic', fontsize = 25, fontweight = 'bold', x = 0.05)

    # increase size of axes ticks 
    for label in axs[0].get_xticklabels() + axs[0].get_yticklabels() + axs[1].get_xticklabels() + axs[1].get_yticklabels():
        label.set_fontsize(15)
        label.set_fontweight('bold')

    # add shaded region of plot where signals are bimodal 
    for j in range(len(axs)):
        axs[j].fill_between(x = np.linspace(0.555, 1.0), y1 = 0.05, y2 = 1.0, color = 'b', alpha = 0.1, zorder = 1)
        axs[j].text(0.78, 0.5,'Bimodal\nSignals', fontsize = 20, fontweight = 'bold', ha = 'center')
        axs[j].vlines(0.555, ymin = 0.05, ymax = 1.0, color = 'k', alpha = 0.5)
        axs[j].hlines(0.05, xmin = 0.555, xmax = 1.0, color = 'k', alpha = 0.5)
        axs[j].grid()


    # define legend using patches
    gray_patch = mpatches.Patch(facecolor = 'gray', edgecolor = 'k', label = 'Less\nProfitable')
    color_patch = mpatches.Patch(facecolor = colors[technology], edgecolor = 'k', label = 'More\nProfitable')

    # legend
    legend_properties = {'weight':'bold', 'size':20}

    # individual legends for each plot
    axs[0].legend(handles = [gray_patch, color_patch], loc = 'upper left', prop = legend_properties)
    axs[1].legend(handles = [gray_patch, color_patch], loc = 'upper left', prop = legend_properties)

    # save figure in this directory if option is selected 
    if save:
        plt.savefig('bimodality_scatter_{0}.png'.format(str(h2_price)), dpi = 300, bbox_inches = 'tight')
        plt.savefig('bimodality_scatter_{0}.pdf'.format(str(h2_price)), dpi = 300, bbox_inches = 'tight')


    plt.show()

    return 

def compute_computation_times(directory, save = False):
    '''
    Returns a table of minimum, average, and maximum computation times in minutes 
    as a data frame.

    Arguments:
        directory: Location of market results files 
            default = '../', folders are one directory up 
        save: boolean, whether or not to save the table as a csv
            default = False, do not save 

    Returns:
        comp_times: pd.DataFrame, minimum, average, and maximum compuitational times by system 
    '''

    # read all computational time file 
    time_10 = pd.read_csv((directory + '/market_results_10/00_time.csv'), header = None)
    time_15 = pd.read_csv((directory + '/market_results_15/00_time.csv'), header = None)
    time_20 = pd.read_csv((directory + '/market_results_20/00_time.csv'), header = None)
    time_25 = pd.read_csv((directory + '/market_results_25/00_time.csv'), header = None)
    time_30 = pd.read_csv((directory + '/market_results_30/00_time.csv'), header = None)  

    # create lists for each system
    NGCC = []
    SOFC = []
    NGCCSOEC = []
    rSOC = []
    SOFCSOEC = []
    SOEC = []

    # save time values from dataframe into lists
    for df in [time_10, time_15, time_20, time_25, time_30]:
        for i ,row in df.iterrows():
            if row[0] == 'model0':
                NGCC.append(row[4])

            elif row[0] == 'model1':
                SOFC.append(row[4])

            elif row[0] == 'model3':
                NGCCSOEC.append(row[4])

            elif row[0] == 'model4':
                rSOC.append(row[4])

            elif row[0] == 'model5':
                SOFCSOEC.append(row[4])

            elif row[0] == 'model6':
                SOEC.append(row[4])

            else:
                raise NotImplementedError()
            
    # make dataframe for the computation times 

    comp_times = pd.DataFrame(index = ['NGCC', 'SOFC', 'NGCC + SOEC', 'rSOC', 'SOFC + SOEC', 'SOEC'], columns = ['Minimum Solving Time (min)', 'Average Solving Time (min)', 'Maximum Solving Time (min)'])

    comp_times.at['NGCC', 'Minimum Solving Time (min)'] = round(np.min(NGCC),2)
    comp_times.at['NGCC', 'Average Solving Time (min)'] = round(np.mean(NGCC),2)
    comp_times.at['NGCC', 'Maximum Solving Time (min)'] = round(max(NGCC),2)

    comp_times.at['SOFC', 'Minimum Solving Time (min)'] = round(min(SOFC),2)
    comp_times.at['SOFC', 'Average Solving Time (min)'] = round(np.mean(SOFC),2)
    comp_times.at['SOFC', 'Maximum Solving Time (min)'] = round(max(SOFC),2)

    comp_times.at['NGCC + SOEC', 'Minimum Solving Time (min)'] = round(min(NGCCSOEC),2)
    comp_times.at['NGCC + SOEC', 'Average Solving Time (min)'] = round(np.mean(NGCCSOEC),2)
    comp_times.at['NGCC + SOEC', 'Maximum Solving Time (min)'] = round(max(NGCCSOEC),2)

    comp_times.at['rSOC', 'Minimum Solving Time (min)'] = round(min(rSOC),2)
    comp_times.at['rSOC', 'Average Solving Time (min)'] = round(np.mean(rSOC),2)
    comp_times.at['rSOC', 'Maximum Solving Time (min)'] = round(max(rSOC),2)

    comp_times.at['SOFC + SOEC', 'Minimum Solving Time (min)'] = round(min(SOFCSOEC),2)
    comp_times.at['SOFC + SOEC', 'Average Solving Time (min)'] = round(np.mean(SOFCSOEC),2)
    comp_times.at['SOFC + SOEC', 'Maximum Solving Time (min)'] = round(max(SOFCSOEC),2)

    comp_times.at['SOEC', 'Minimum Solving Time (min)'] = round(min(SOEC),2)
    comp_times.at['SOEC', 'Average Solving Time (min)'] = round(np.mean(SOEC),2)
    comp_times.at['SOEC', 'Maximum Solving Time (min)'] = round(max(SOEC),2)

    if save:
        comp_times.to_csv('comp_times.csv')

    return comp_times