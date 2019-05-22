# -*- coding: utf-8 -*-
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
from enum import Enum
from bokeh.models import Arrow, OpenHead, NormalHead
from bokeh.models import Label, HoverTool, CircleCross, X
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from statistics import median
import random


class HENStreamType(Enum):
    """Enum type defining hot and cold streams
    """
    hot = 1
    cold = 2
    hot_utility = 3
    cold_utility = 4



def turn_off_grid_and_axes_ticks(plot):
    """Turn off axis ticks and grid lines on a bokeh figure object.

    Args:
        plot: bokeh.plotting.plotting.figure instance. 

    Returns:
        modified bokeh.plotting.plotting.figure instance.

    Raises:
        None

    """
    plot.xaxis.major_tick_line_color = None
    plot.xaxis.minor_tick_line_color = None
    plot.yaxis.major_tick_line_color = None
    plot.yaxis.minor_tick_line_color = None
    plot.xaxis.visible = None
    plot.yaxis.visible = None
    plot.xgrid.grid_line_color = None
    plot.ygrid.grid_line_color = None
    return plot        


def get_color_dictionary(set_to_color):
    """Given a set, return a dictionary of the form: 

    .. code-block:: python

        {'set_member': valid_bokeh_color}
    ..

        Args:
            set_to_color: set of unique elements, e.g: [1,2,3] or ["1", "2", "3"]
        
        Returns:
            Dictionary of the form: 

            .. code-block:: python

                {'set_member': valid_bokeh_color}


        Raises:
            None            
    """
    color_dict = {}
    for stg in set_to_color:
        def r(): return random.randint(0, 255)
        color_dict[stg] = '#%02X%02X%02X' % (r(), r(), r())
    return color_dict


def add_module_markers_to_heat_exchanger_plot(plot,
                                              x, y, modules, line_color,
                                              fill_color,
                                              mark_modules_with_tooltips):
    """Plot module markers as tooltips to a heat exchanger network diagram.
        
        Args:
            plot : bokeh.plotting.plotting.figure instance. 
            x : x-axis coordinate of module marker tooltip.
            y : y-axis coordinate of module marker tooltip.
            modules : dict containing modules.
            line_color : color of border of the module marker.
            fill_color : color inside the module marker.
            mark_modules_with_tooltips : whether to add tooltips to plot or not (currently not utilized).
        
        Returns: 
            bokeh.plotting.plotting.figure instance with module markers added. 

        Raises:
            None
    """
    module_text_template = '{0} x {1} m^2, '
    module_text = ''.join([module_text_template.format(
        modules[m], m) for m in modules])[0:-2]
    source = ColumnDataSource(data=dict(x=[x],
                                        y=[y],
                                        modules=['{0}'.format(
                                            str(module_text))]))
    marker_module = X(x=x, y=y,
                      line_color=line_color,
                      fill_color=fill_color,
                      size=5)

    glyph_marker = plot.add_glyph(
        source_or_glyph=source, glyph=marker_module)
    hover_modules = HoverTool(renderers=[glyph_marker],
                              tooltips=[('Modules', '@modules')])
    plot.add_tools(hover_modules)
    return plot


def is_hot_or_cold_utility(exchanger):
    """Return if an exchanger is a hot or a cold utility by checking if it has the key `utility_type`.
        
        Args:
            exchanger: dict representing the exchanger. 

        Returns:        
            True if `utility_type` in the `exchanger` dict passed.

        Raises:
            None            
    """
    return 'utility_type' in exchanger


def get_stream_y_values(exchangers, hot_streams, cold_streams, y_stream_step=1):
    """Return a dict containing the layout of the heat exchanger diagram including any stage splits.

        Args:
            exchangers: List of exchangers where each exchanger is a dict of the form:

                        .. code-block:: python
                        
                            {'hot': 'H2', 'cold': 'C1', 'Q': 1400, 'A': 159, 'annual_cost': 28358, 
                            'stg': 2}
                        

                        where hot is the hot stream name, cold is the cold stream name, A is the area (in m^2), annual_cost
                        is the annual cost in $, Q is the amount of heat transferred from one stream
                        to another in a given exchanger and stg is the stage the exchanger belongs to.
                        Additionally a 'utility_type' can specify if we draw the cold stream as water (:class:`idaes.vis.plot_utils.HENStreamType.cold_utility`)
                        or the hot stream as steam (:class:`idaes.vis.plot_utils.HENStreamType.hot_utility`). 

                        Additionally, the exchanger could have the key 'modules', like this:

                        .. code-block:: python
                        
                            {'hot': 'H1', 'cold': 'C1', 'Q': 667, 'A': 50, 'annual_cost': 10979, 'stg': 3, 
                            'modules': {10: 1, 20: 2}} 


            hot_streams: List of dicts representing hot streams where each item is a dict of the form:

                         .. code-block:: python
                            
                            {'name':'H1', 'temps': [443, 435, 355, 333], 'type': HENStreamType.hot}


            cold_streams: List of dicts representing cold streams where each item is a dict of the form:

                          .. code-block:: python

                            {'name':'H1', 'temps': [443, 435, 355, 333], 'type': HENStreamType.hot}                             


            y_stream_step: how many units on the HEN diagram to leave between each stream (or sub-stream) and the one above it. Defaults to 1. 
        
        Returns: 
            Tuple containing 3 dictionaries to be used when plotting the HEN: 
                * stream_y_values_dict : a dict of each stream name as key and value being a dict of the form

                .. code-block:: python
                
                    {'default_y_value': 2, 'split_y_values': [1,3]}.
                ..

                This indicates what the default y value of this stream will be on the diagram and what values we'll 
                use when it splits. 
            
                * hot_split_streams : list of tuples of the form (a,b) where a is a hot stream name and b is the max.
                times it will split over all the stages.
            
                * cold_split_streams : list of tuples of the form (a,b) where a is a cold stream name and b is the max.
                times it will split over all the stages.

        Raises:
            None                
    """
    STREAM_NAME_INDEX = 0
    STREAM_SPLIT_COUNT = 3
    stages = set([exchanger['stg'] for exchanger in exchangers])
    stream_y_values = {}
    split_streams = []

    stream_names = [list(stream.keys())[0] for stream in hot_streams] + \
        [list(stream.keys())[0] for stream in cold_streams]
    for name in stream_names:
        stream_y_values[name] = {
            'default_y_value': -1, 'split_y_values': []}

    for stage in stages:
        hot_stage_streams = [exchanger['hot']
                             for exchanger in exchangers if exchanger['stg'] == stage
                             and not is_hot_or_cold_utility(exchanger)]
        cold_stage_streams = [exchanger['cold']
                              for exchanger in exchangers if exchanger['stg'] == stage
                              and not is_hot_or_cold_utility(exchanger)]
        hot_streams_to_split = set(
            [(stream, HENStreamType.hot, stage, hot_stage_streams.count(stream)) for stream in hot_stage_streams if hot_stage_streams.count(stream) > 1])
        cold_streams_to_split = set(
            [(stream, HENStreamType.cold, stage, cold_stage_streams.count(stream)) for stream in cold_stage_streams if cold_stage_streams.count(stream) > 1])
        split_streams.extend(list(hot_streams_to_split))
        split_streams.extend(list(cold_streams_to_split))

    split_stream_names = [split_stream_record[STREAM_NAME_INDEX]
                          for split_stream_record in split_streams]
    cold_split_streams = [
        split_stream_record for split_stream_record in split_streams if split_stream_record[1] == HENStreamType.cold]
    hot_split_streams = [
        split_stream_record for split_stream_record in split_streams if split_stream_record[1] == HENStreamType.hot]

    max_y = 0
    for stream in cold_streams:
        stream_name = list(stream.keys())[0]
        if stream_name in split_stream_names:
            max_split_count = max([split_stream_record[STREAM_SPLIT_COUNT]
                                   for split_stream_record in cold_split_streams])
            stream_split_y_values = []
            for i in range(max_y, max_y+max_split_count):
                stream_split_y_values.append(i)
            stream_y_values[stream_name]['split_y_values'] = sorted(
                stream_split_y_values)
            if len(stream_split_y_values) % 2 != 0:
                stream_y_values[stream_name]['default_y_value'] = median(
                    stream_split_y_values[1:])
            else:
                stream_y_values[stream_name]['default_y_value'] = median(
                    stream_split_y_values)
            max_y = max(stream_split_y_values) + y_stream_step
        else:
            stream_y_values[stream_name]['default_y_value'] = max_y
            max_y += y_stream_step

    for stream in hot_streams:
        stream_name = list(stream.keys())[0]
        if stream_name in split_stream_names:
            max_split_count = max([split_stream_record[STREAM_SPLIT_COUNT]
                                   for split_stream_record in hot_split_streams])
            stream_split_y_values = []
            for i in range(max_y, max_y+max_split_count):
                stream_split_y_values.append(i)
            stream_y_values[stream_name]['split_y_values'] = sorted(
                stream_split_y_values)
            if len(stream_split_y_values) % 2 != 0:
                stream_y_values[stream_name]['default_y_value'] = median(
                    stream_split_y_values[1:])
            else:
                stream_y_values[stream_name]['default_y_value'] = median(
                    stream_split_y_values)
            max_y = max(stream_split_y_values) + y_stream_step
        else:
            stream_y_values[stream_name]['default_y_value'] = max_y
            max_y += y_stream_step
    stream_y_values = sorted(stream_y_values.items(
    ), key=lambda k_v: k_v[1]['default_y_value'])
    stream_y_values_dict = {}
    for item in stream_y_values:
        stream_y_values_dict[item[0]] = item[1]
    return stream_y_values_dict, hot_split_streams, cold_split_streams


def plot_line_segment(plot, x_start, x_end, y_start, y_end, color="white", legend=None):
    """Plot a line segment on a bokeh figure. 

        Args:
            plot: bokeh.plotting.plotting.figure instance. 
            x_start: x-axis coordinate of 1st point in line.         
            x_end: x-axis coordinate of 2nd point in line.         
            y_start: y-axis coordinate of 1st point in line.         
            y_end: y-axis coordinate of 2nd point in line.
            color: color of line (defaults to white).
            legend: what legend to associate with (defaults to None).

        Returns:
            modified bokeh.plotting.plotting.figure instance with line added.  
        
        Raises:
            None                   
    """
    plot.line([x_start, x_end], [y_start, y_end],
              line_width=2,
              color=color,
              legend=legend)
    return plot


def plot_stream_arrow(plot, line_color,
                      stream_arrow_temp,
                      temp_label_font_size,
                      x_start, x_end, y_start, y_end,
                      stream_name=None):
    """Plot a stream arrow for the heat exchanger network diagram.

        Args:
            plot: bokeh.plotting.plotting.figure instance.
            line_color: color of arrow (defaults to white).
            stream_arrow_temp: Tempreature of the stream to be plotted.
            temp_label_font_size: font-size of the temperature label to be added.
            x_start: x-axis coordinate of arrow base.
            x_end: x-axis coordinate of arrow head.
            y_start: y-axis coordinate of arrow base.
            y_end: y-axis coordinate of arrow head.
            stream_name: Name of the stream to add as a label to arrow (defaults to None).

        Returns:
            modified bokeh.plotting.plotting.figure instance with stream arrow added.
        
        Raises:
            None            
    """
    plot.add_layout(Arrow(line_color=line_color,
                          end=OpenHead(line_color=line_color,
                                       line_width=2, size=10),
                          x_start=x_start, y_start=y_start,
                          x_end=x_end, y_end=y_end))

    temp_end_label = Label(
        x=x_end, y=y_end,
        x_offset=-20, y_offset=8,
        text_font_size=temp_label_font_size,
        text=str(stream_arrow_temp))
    plot.add_layout(temp_end_label)

    if stream_name:
        stream_name_label = Label(
            x=x_end, y=y_end,
            x_offset=-10, y_offset=-12,
            text_font_size=temp_label_font_size,
            text=stream_name)
        plot.add_layout(stream_name_label)
    return plot


def add_exchanger_labels(plot, x, y_start, y_end, label_font_size, exchanger,
                         module_marker_line_color,
                         module_marker_fill_color,
                         mark_modules_with_tooltips):
    """Plot exchanger labels for an exchanger (for Q and A) on a heat exchanger network diagram and add module markers (if needed).

        Args:
            plot: bokeh.plotting.plotting.figure instance. 
            label_font_size: font-size for labels.
            x: x-axis coordinate of exchanger (exchangers are vertical lines so we just need 1 x-value)
            y_start: y-axis coordinate of exchanger start.         
            y_end: y-axis coordinate of exchanger end.            
            exchanger: exchanger dictionary of the form:

                       .. code-block:: python

                           {'hot': 'H2', 'cold': 'C1', 'Q': 1400, 'A': 159, 'annual_cost': 28358, 
                            'stg': 2}

            
            module_marker_line_color : color of border of the module marker.
            module_marker_fill_color : color inside the module marker.
            mark_modules_with_tooltips : whether to add tooltips to plot or not (currently not utilized).

        Returns:
            modified bokeh.plotting.plotting.figure instance with labels added.     

        Raises:
            None                   
    """
    label_q = Label(
        x=x, y=y_start,
        x_offset=-10,
        text_font_size=label_font_size,
        text='{0} kW'.format(exchanger['Q']))
    plot.add_layout(label_q)
    label_a = Label(x=x, y=y_end,
                    y_offset=-20, x_offset=-20,
                    text_font_size=label_font_size,
                    text='{0} m^2'.format(exchanger['A']))
    plot.add_layout(label_a)
    # Mark modules:
    if 'modules' in exchanger:
        plot = add_module_markers_to_heat_exchanger_plot(plot=plot,
                                                             x=x,
                                                             y=y_end,
                                                             modules=exchanger['modules'],
                                                             line_color=module_marker_line_color,
                                                             fill_color=module_marker_fill_color,
                                                             mark_modules_with_tooltips=mark_modules_with_tooltips)
    return plot    