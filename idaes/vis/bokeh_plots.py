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
"""
Bokeh plots.
"""
# stdlib
import warnings

# third-party
from IPython.core.display import HTML
from IPython.core.display import display as idisplay
from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.models import LabelSet

# package
from idaes.vis.plotbase import PlotBase
from idaes.vis.plotutils import *


class BokehPlot(PlotBase):
    def __init__(self, current_plot=None):
        """Set the current plot in this constructor.

        Args:
            current_plot: A bokeh plot object.

        """
        active_dev_warning = '''The visualization library is still in active development and we hope
        to improve on it in future releases. Please use its functionality at your own discretion.'''
        warnings.warn(active_dev_warning, stacklevel=3)
        super().__init__(current_plot)
        self.html = ''
        self.filename = ''

    def annotate(self, x, y, label):
        """Annotate a plot with a given point and a label.

        Args:
            x: Value of independent variable.
            y: Value of dependent variable.
            label: Text label.

        Returns:
            None

        Raises:
            None
        """
        source = ColumnDataSource(data=dict(x=[x], y=[y], labels=[label]))
        self.current_plot.circle(x, y, line_color='red', fill_color='red')
        labels = LabelSet(
            x='x',
            y='y',
            text='labels',
            level='glyph',
            source=source,
            render_mode='canvas',
            text_font_size='9pt',
        )
        self.current_plot.add_layout(labels)

    def resize(self, height=-1, width=-1):
        """Resize a plot's height and width.

        Args:
            height: Height in screen units.
            width: Width in screen units.

        Returns:
            None

        Raises:
            None
        """
        if not self.current_plot:
            return
        if height != -1:
            self.current_plot.plot_height = height
        if width != -1:
            self.current_plot.plot_width = width

    def save(self, destination):
        """Save the current plot object to HTML in filepath provided by destination.

        Args:
            destination: Valid file path to save HTML to.

        Returns:
            filename where HTML is saved.

        Raises:
            None
        """
        if not self.html:
            self.html = file_html(self.current_plot, CDN, None)
        with open(destination, 'w') as output_file:
            output_file.write(self.html)
            self.filename = destination
        return self.filename

    def show(self, in_notebook=True):
        """Display plot in a Jupyter notebook.

        Args:
            self: Plot object.
            in_notebook: Display in Jupyter notebook or generate HTML file.

        Returns:
            None

        Raises:
            None
        """
        if in_notebook:
            self.html = file_html(self.current_plot, CDN, None)
            idisplay(HTML(self.html))
        else:
            show(self.current_plot)


class HeatExchangerNetwork(BokehPlot):
    def __init__(
        self,
        exchangers,
        stream_list,
        mark_temperatures_with_tooltips=False,
        mark_modules_with_tooltips=False,
        stage_width=2,
        y_stream_step=1,
    ):
        """Plot a heat exchanger network diagram.

        Args:
            exchangers: List of exchangers where each exchanger is a dict of the form:

                        .. code-block:: python

                            {'hot': 'H2', 'cold': 'C1', 'Q': 1400, 'A': 159, 'annual_cost': 28358,
                            'stg': 2}


                        where `hot` is the hot stream name, `cold` is the cold stream name, `A` is the area (in m^2), `annual_cost`
                        is the annual cost in $, `Q` is the amount of heat transferred from one stream
                        to another in a given exchanger and `stg` is the stage the exchanger belongs to.
                        The `utility_type`, if present, will specify if we plot the cold stream as water (:class:`idaes.vis.plot_utils.HENStreamType.cold_utility`)
                        or the hot stream as steam (:class:`idaes.vis.plot_utils.HENStreamType.hot_utility`).

                        Additionally, the exchanger could have the key `modules`, like this:

                        .. code-block:: python

                            {'hot': 'H1', 'cold': 'C1', 'Q': 667, 'A': 50, 'annual_cost': 10979, 'stg': 3,
                            'modules': {10: 1, 20: 2}}


                        The value of this key is a dictionary where each key is a module area and each value is how many
                        modules of that area are in the exchanger. It's indicated as a tooltip on the resulting diagram.

                        If a stream is involved in multiple exchanges in teh same stage, the stream will split into multiple
                        sub-streams with each sub-stream carrying one of the exchanges.

            stream_list: List of dicts representing streams where each item is a dict of the form:

                         .. code-block:: python

                            {'name':'H1', 'temps': [443, 435, 355, 333], 'type': HENStreamType.hot}


            mark_temperatures_with_tooltips: if True, we plot the stream temperatures and assign
                                             hover tooltips to them. Otherwise, we label them with text labels.

            mark_modules_with_tooltips: if True, we plot markers for modules (if present) and assign
                                             hover tooltips to them. Otherwise, we don't add module info.

            stage_width: How many units to use for each stage in the diagram (defaults to 2).

            y_stream_step: How many units to use to separate each stream/sub-stream from the next (defaults to 1).
        """

        total_streams_tuple = (-1, len(stream_list) + 1)
        total_cost = sum([exchanger['annual_cost'] for exchanger in exchangers])

        REGULAR_FONT_SIZE = "10pt"
        SMALL_FONT_SIZE = "5pt"
        COLD_STREAM_COLOR = "blue"
        HOT_STREAM_COLOR = "red"
        STREAM_TEMP_MARKER_FILL_COLOR = "white"
        STREAM_TEMP_MARKER_LINE_COLOR = HOT_STREAM_COLOR

        STAGE_WIDTH = stage_width
        STAGE_SEPARATOR = 0.5
        Y_STREAM_STEP = y_stream_step

        p = figure(
            x_range=total_streams_tuple,
            y_range=total_streams_tuple,
            title='Annualized Investment Cost = ${0}'.format(total_cost),
        )
        p.title.align = 'center'

        hot_streams = [
            {k['name']: k['temps']}
            for k in stream_list
            if k['type'] == HENStreamType.hot
        ]
        cold_streams = [
            {k['name']: k['temps']}
            for k in stream_list
            if k['type'] == HENStreamType.cold
        ]
        stage_set = set(sorted([exchanger['stg'] for exchanger in exchangers]))
        color_stage_dict = get_color_dictionary(stage_set)
        stream_y_values, hot_split_streams, cold_split_streams = get_stream_y_values(
            exchangers, hot_streams, cold_streams, y_stream_step=Y_STREAM_STEP
        )

        stage_x_limits = {}
        max_x = 0
        for stage in stage_set:
            if stage not in stage_x_limits:
                stage_x_limits[stage] = {
                    'x_start': max_x,
                    'x_exchanger_start': max_x + STAGE_SEPARATOR,
                    'x_exchanger_end': max_x + STAGE_WIDTH,
                    'x_true_end': max_x + STAGE_WIDTH + STAGE_SEPARATOR,
                }
                max_x = max_x + STAGE_WIDTH

        hot_stream_names = [list(stream.keys())[0] for stream in hot_streams]
        cold_stream_names = [list(stream.keys())[0] for stream in cold_streams]

        for i, stage in enumerate(stage_set):
            stage_exchangers = [
                exchanger
                for exchanger in exchangers
                if exchanger['stg'] == stage and not is_hot_or_cold_utility(exchanger)
            ]

            hot_split_streams_by_stage = [
                hot_split_stream[0]
                for hot_split_stream in hot_split_streams
                if hot_split_stream[2] == stage
            ]
            cold_split_streams_by_stage = [
                cold_split_stream[0]
                for cold_split_stream in cold_split_streams
                if cold_split_stream[2] == stage
            ]

            used_split_y_values_for_stage = []
            is_first_stage = i == 0
            is_last_stage = i == (len(stage_set) - 1)

            # Figure out distance between each exchanger for this stage
            # leaving some leeway for closing out split stream rectangles:
            start_x = stage_x_limits[stage]['x_exchanger_start']
            end_x = stage_x_limits[stage]['x_exchanger_end'] - STAGE_SEPARATOR
            if len(stage_exchangers) > 0:
                dist_between_exchanger = (end_x - start_x) / len(stage_exchangers)
            max_x = start_x

            for exchanger in stage_exchangers:
                exchanger_x_start = max_x + dist_between_exchanger
                hot_stream = exchanger['hot']
                cold_stream = exchanger['cold']
                if hot_stream in hot_split_streams_by_stage:
                    # Give me the 1st unused hot stream split y value for this stage:
                    for y_value in stream_y_values[hot_stream]['split_y_values']:
                        if y_value not in used_split_y_values_for_stage:
                            used_split_y_values_for_stage.append(y_value)
                            exchanger_y_start = y_value
                            break
                else:
                    # Just assign it to the default y value:
                    exchanger_y_start = stream_y_values[hot_stream]['default_y_value']

                if cold_stream in cold_split_streams_by_stage:
                    # Give me the 1st unused cold stream split y value for this stage:
                    for y_value in stream_y_values[cold_stream]['split_y_values']:
                        if y_value not in used_split_y_values_for_stage:
                            used_split_y_values_for_stage.append(y_value)
                            exchanger_y_end = y_value
                            break
                else:
                    # Just assign it to the default y value:
                    exchanger_y_end = stream_y_values[cold_stream]['default_y_value']

                p = plot_line_segment(
                    p,
                    exchanger_x_start,
                    exchanger_x_start,
                    exchanger_y_start,
                    exchanger_y_end,
                    color=color_stage_dict[stage],
                    legend="stage {0}".format(str(stage)),
                )
                p = add_exchanger_labels(
                    p,
                    exchanger_x_start,
                    exchanger_y_start,
                    exchanger_y_end,
                    SMALL_FONT_SIZE,
                    exchanger,
                    STREAM_TEMP_MARKER_LINE_COLOR,
                    STREAM_TEMP_MARKER_FILL_COLOR,
                    mark_modules_with_tooltips,
                )
                max_x += dist_between_exchanger
            for stream in stream_y_values.keys():
                color = (
                    HOT_STREAM_COLOR
                    if stream in hot_stream_names
                    else COLD_STREAM_COLOR
                )
                split_streams_list = (
                    hot_split_streams
                    if stream in hot_stream_names
                    else cold_split_streams
                )
                if (
                    stream not in hot_split_streams_by_stage
                    and stream not in cold_split_streams_by_stage
                ):
                    # Default y stream:
                    p = plot_line_segment(
                        p,
                        stage_x_limits[stage]['x_start'],
                        stage_x_limits[stage]['x_true_end'],
                        stream_y_values[stream]['default_y_value'],
                        stream_y_values[stream]['default_y_value'],
                        color=color,
                    )
                else:
                    # Split y stream:
                    split_count = [
                        stream_record[3]
                        for stream_record in split_streams_list
                        if stream_record[2] == stage and stream_record[0] == stream
                    ][0]
                    split_y_values = stream_y_values[stream]['split_y_values'][
                        0:split_count
                    ]
                    vertical_closer_y_start = max(
                        stream_y_values[stream]['default_y_value'], max(split_y_values)
                    )
                    vertical_closer_y_end = min(
                        stream_y_values[stream]['default_y_value'], min(split_y_values)
                    )
                    p = plot_line_segment(
                        p,
                        stage_x_limits[stage]['x_exchanger_start'],
                        stage_x_limits[stage]['x_exchanger_start'],
                        vertical_closer_y_start,
                        vertical_closer_y_end,
                        color=color,
                    )

                    p = plot_line_segment(
                        p,
                        stage_x_limits[stage]['x_exchanger_end'],
                        stage_x_limits[stage]['x_exchanger_end'],
                        vertical_closer_y_start,
                        vertical_closer_y_end,
                        color=color,
                    )
                    for y in sorted(split_y_values):
                        p = plot_line_segment(
                            p,
                            stage_x_limits[stage]['x_exchanger_start'],
                            stage_x_limits[stage]['x_exchanger_end'],
                            y,
                            y,
                            color=color,
                        )
            # 1st stage:
            # For hot streams:
            # - Plot line segments, stream names and highest tempreature label.
            # For cold streams:
            # - Plot arrows and lowest temperature label.
            # Plot hot utility exchanges.
            if is_first_stage:
                for stream in stream_y_values.keys():
                    temp_list = [
                        stream_dict['temps']
                        for stream_dict in stream_list
                        if stream_dict['name'] == stream
                    ][0]
                    if stream in hot_stream_names:
                        # Plot line segments that start/end a hot/cold stream
                        # and add relevant labels in the 1st stage:
                        p = plot_line_segment(
                            p,
                            stage_x_limits[stage]['x_start'],
                            stage_x_limits[stage]['x_exchanger_start'],
                            stream_y_values[stream]['default_y_value'],
                            stream_y_values[stream]['default_y_value'],
                            color=HOT_STREAM_COLOR,
                        )

                        temp_start_label = Label(
                            x=stage_x_limits[stage]['x_start'],
                            y=stream_y_values[stream]['default_y_value'],
                            x_offset=-10,
                            y_offset=2,
                            text_font_size=REGULAR_FONT_SIZE,
                            text=str(temp_list[0]),
                        )
                        p.add_layout(temp_start_label)

                        stream_name_label = Label(
                            x=stage_x_limits[stage]['x_start'],
                            y=stream_y_values[stream]['default_y_value'],
                            x_offset=-10,
                            y_offset=-12,
                            text_font_size=REGULAR_FONT_SIZE,
                            text=stream,
                        )
                        p.add_layout(stream_name_label)

                    if stream in cold_stream_names:
                        p = plot_stream_arrow(
                            p,
                            COLD_STREAM_COLOR,
                            temp_list[-1],
                            REGULAR_FONT_SIZE,
                            stage_x_limits[stage]['x_exchanger_start'],
                            stage_x_limits[stage]['x_start'],
                            stream_y_values[stream]['default_y_value'],
                            stream_y_values[stream]['default_y_value'],
                        )
                # plot steam exchange:
                steam_exchanger = [
                    exchanger
                    for exchanger in exchangers
                    if (
                        'utility_type' in exchanger
                        and exchanger['utility_type'] == HENStreamType.hot_utility
                    )
                ]
                if len(steam_exchanger) > 0:
                    cold_stream_name = steam_exchanger[0]['cold']
                    y_steam = stream_y_values[cold_stream_name]['default_y_value']
                    p.add_layout(
                        Arrow(
                            line_color=color_stage_dict[stage],
                            end=OpenHead(
                                line_color=color_stage_dict[stage],
                                line_width=2,
                                size=10,
                            ),
                            x_start=0.25,
                            y_start=y_steam - 0.25,
                            x_end=0.25,
                            y_end=y_steam + 0.25,
                        )
                    )
                    p = plot_line_segment(
                        p,
                        0.25,
                        0.25,
                        y_steam - 0.25,
                        y_steam + 0.25,
                        color=color_stage_dict[stage],
                        legend="stage {0}".format(str(stage)),
                    )
                    label_s = Label(
                        x=0.25,
                        y=y_steam - 0.25,
                        x_offset=-10,
                        y_offset=-15,
                        text_font_size=SMALL_FONT_SIZE,
                        text=steam_exchanger[0]['hot'],
                    )
                    p.add_layout(label_s)
                    label_q = Label(
                        x=0.25,
                        y=y_steam + 0.25,
                        x_offset=-5,
                        text_font_size=SMALL_FONT_SIZE,
                        text='{0} kW'.format(steam_exchanger[0]['Q']),
                    )
                    p.add_layout(label_q)

            # Last stage:
            # For cold streams:
            # - Plot line segments, stream names and highest tempreature label.
            # For hot streams:
            # - Plot arrows and lowest temperature label.
            # Plot cold utility exchanges.
            elif is_last_stage:
                for stream in stream_y_values.keys():
                    temp_list = [
                        stream_dict['temps']
                        for stream_dict in stream_list
                        if stream_dict['name'] == stream
                    ][0]
                    if stream in cold_stream_names:
                        # Plot line segments that start/end a hot/cold stream
                        # and add relevant labels in the 1st stage:
                        p = plot_line_segment(
                            p,
                            stage_x_limits[stage]['x_exchanger_end'],
                            stage_x_limits[stage]['x_true_end'],
                            stream_y_values[stream]['default_y_value'],
                            stream_y_values[stream]['default_y_value'],
                            color=COLD_STREAM_COLOR,
                        )

                        temp_start_label = Label(
                            x=stage_x_limits[stage]['x_true_end'],
                            y=stream_y_values[stream]['default_y_value'],
                            x_offset=-10,
                            y_offset=2,
                            text_font_size=REGULAR_FONT_SIZE,
                            text=str(temp_list[0]),
                        )
                        p.add_layout(temp_start_label)

                        stream_name_label = Label(
                            x=stage_x_limits[stage]['x_true_end'],
                            y=stream_y_values[stream]['default_y_value'],
                            x_offset=-10,
                            y_offset=-12,
                            text_font_size=REGULAR_FONT_SIZE,
                            text=stream,
                        )
                        p.add_layout(stream_name_label)

                    if stream in hot_stream_names:
                        p = plot_stream_arrow(
                            p,
                            HOT_STREAM_COLOR,
                            temp_list[-1],
                            REGULAR_FONT_SIZE,
                            stage_x_limits[stage]['x_exchanger_end'],
                            stage_x_limits[stage]['x_true_end'],
                            stream_y_values[stream]['default_y_value'],
                            stream_y_values[stream]['default_y_value'],
                        )
                # plot water exchange:
                water_exchanger = [
                    exchanger
                    for exchanger in exchangers
                    if (
                        'utility_type' in exchanger
                        and exchanger['utility_type'] == HENStreamType.cold_utility
                    )
                ]
                if len(water_exchanger) > 0:
                    hot_stream_name = water_exchanger[0]['hot']
                    y_water = stream_y_values[hot_stream_name]['default_y_value']
                    x_water = (
                        stage_x_limits[stage]['x_exchanger_end']
                        + (
                            stage_x_limits[stage]['x_true_end']
                            - stage_x_limits[stage]['x_exchanger_end']
                        )
                        / 2
                    )
                    p.add_layout(
                        Arrow(
                            line_color=color_stage_dict[stage],
                            end=OpenHead(
                                line_color=color_stage_dict[stage],
                                line_width=2,
                                size=10,
                            ),
                            x_start=x_water,
                            y_start=y_water + 0.25,
                            x_end=x_water,
                            y_end=y_water - 0.25,
                        )
                    )
                    p = plot_line_segment(
                        p,
                        x_water,
                        x_water,
                        y_water + 0.25,
                        y_water - 0.25,
                        color=color_stage_dict[stage],
                        legend="stage {0}".format(str(stage)),
                    )
                    label_s = Label(
                        x=x_water,
                        y=y_water - 0.25,
                        x_offset=-10,
                        y_offset=-15,
                        text_font_size=SMALL_FONT_SIZE,
                        text=water_exchanger[0]['cold'],
                    )
                    p.add_layout(label_s)
                    label_q = Label(
                        x=x_water,
                        y=y_water + 0.25,
                        x_offset=-5,
                        text_font_size=SMALL_FONT_SIZE,
                        text='{0} kW'.format(water_exchanger[0]['Q']),
                    )
                    p.add_layout(label_q)
            # Otherwise, plot line segments after the exchangers to finish off a stage:
            else:
                for stream in stream_y_values.keys():
                    color = (
                        COLD_STREAM_COLOR
                        if stream in cold_stream_names
                        else HOT_STREAM_COLOR
                    )
                    p = plot_line_segment(
                        p,
                        stage_x_limits[stage]['x_exchanger_end'],
                        stage_x_limits[stage]['x_true_end'],
                        stream_y_values[stream]['default_y_value'],
                        stream_y_values[stream]['default_y_value'],
                        color=color,
                    )
            # Plot stage temperatures for each stream
            for stream in stream_y_values.keys():
                temp_list = [
                    stream_dict['temps']
                    for stream_dict in stream_list
                    if stream_dict['name'] == stream
                ][0]
                if len(temp_list) > 1:
                    temp_list = temp_list[1:-1]
                if stage <= len(temp_list):
                    temp_x = (
                        stage_x_limits[stage]['x_exchanger_end']
                        + (
                            stage_x_limits[stage]['x_true_end']
                            - stage_x_limits[stage]['x_exchanger_end']
                        )
                        / 2
                    )
                    temp_to_mark = temp_list[stage - 1]
                    marker_temp = CircleCross(
                        x=temp_x,
                        y=stream_y_values[stream]['default_y_value'],
                        line_color=STREAM_TEMP_MARKER_LINE_COLOR,
                        fill_color=STREAM_TEMP_MARKER_FILL_COLOR,
                        size=10,
                    )
                    if mark_temperatures_with_tooltips:
                        source = ColumnDataSource(
                            data=dict(
                                x=[temp_x],
                                y=[stream_y_values[stream]['default_y_value']],
                                temps=['{0}'.format(str(temp_to_mark))],
                            )
                        )
                        glyph_marker = p.add_glyph(
                            source_or_glyph=source, glyph=marker_temp
                        )
                        hover_temps = HoverTool(
                            renderers=[glyph_marker],
                            tooltips=[('Temperature', '@temps')],
                        )
                        p.add_tools(hover_temps)
                    else:
                        glyph_marker = p.add_glyph(marker_temp)
                        temp_label = Label(
                            x=temp_x,
                            y=stream_y_values[stream]['default_y_value'],
                            text_font_size=SMALL_FONT_SIZE,
                            text='{0}'.format(str(temp_to_mark)),
                        )
                        p.add_layout(temp_label)

        p = turn_off_grid_and_axes_ticks(p)
        super().__init__(current_plot=p)


class ProfilePlot(BokehPlot):
    def __init__(
        self,
        data_frame,
        x='',
        y=None,
        title='',
        xlab='',
        ylab='',
        y_axis_type='auto',
        legend=None,
    ):
        """
        A profile plot includes 2 dependent variables and a single
        independent variable. Based on the Jupyter notebook `here <https://github.com/IDAES/model_contrib/blob/master/examples/mea_simple/mea_example_nb_01.ipynb>`_.

        Args:
            data_frame: a data frame with keys contained in x and y.
            x: Key in data-frame to use as x-axis.
            y: Keys in data-frame to use as y-axis.
            title: Title for a plot.
            xlab: Label for x-axis.
            ylab: Label for y-axis.
            y_axis_type: Specify "log" to pass logarithmic scale.
            legend : List of strings matching y.

        Returns:
            Plot object on success.

        Raises:
            ValueError: Bad input parameters
        """
        # validate inputs
        y = y or []  # coerce to list
        valid, msg = self.validate(data_frame, x, y, legend=legend)
        if not valid:
            raise ValueError(f"Invalid parameters: {msg}")
        # output to static HTML file (commented out for now)
        # output_file("profile_plot.html")
        p = figure(
            tools="pan,box_zoom,reset,save",
            title=title,
            x_axis_label=xlab,
            y_axis_label=ylab,
            y_axis_type=y_axis_type,
        )
        # Plotting both y1 and y2 against x:
        p.line(
            data_frame[x],
            data_frame[y[0]],
            legend=legend[0],
            line_color="green",
        )
        p.line(
            data_frame[x], data_frame[y[1]], legend=legend[1], line_color="blue"
        )
        super().__init__(current_plot=p)


## Not implemented:
# class PropertyModelPlot(BokehPlot):
#     def __init__(self, data_frame, x='', y=[], title='', xlab='',
#                        ylab='', y_axis_type='auto', legend=[]):
#         """Draw pressure/enthalpy plots for different levels of temperature.
#
#         Args:
#             data_frame: a data frame with keys contained in x and y.
#             x: Key in data-frame to plot on x-axis.
#             y: Keys in data-frame to plot on y-axis.
#             title: Title for a plot.
#             xlab: Label for x-axis.
#             ylab: Label for y-axis.
#             y_axis_type: Specify "log" to pass logarithmic scale.
#             legend : List of strings matching y.
#
#         Returns:
#             Plot object on success.
#
#         Raises:
#             MissingVariablesException: Dependent variable or their data
#                                         not passed.
#             BadDataFrameException: No data-frame was generated for
#                                     the model object.
#         """
#         pass
#
#
# class GoodnessOfFitPlot(BokehPlot):
#     def __init__(self, data_frame, x='', y=[], title='', xlab='',
#                         ylab='', y_axis_type='auto', legend=[]):
#         """Draw y against predicted value (y^) and display (calculate?) value of R^2.
#
#         Args:
#             data_frame: a data frame with keys contained in x and y.
#             x: Key in data-frame to use as x-axis.
#             y: Keys in data-frame to plot on y-axis.
#             title: Title for a plot.
#             xlab: Label for x-axis.
#             ylab: Label for y-axis.
#             y_axis_type: Specify "log" to pass logarithmic scale.
#             legend : List of strings matching y.
#
#         Returns:
#             Plot object on success.
#         Raises:
#             MissingVariablesException: Dependent variable or their data
#                                         not passed.
#             BadDataFrameException: No data-frame was generated for
#                                     the model object.
#         """
#         pass
#
#
# class TradeoffPlot(BokehPlot):
#     def __init__(self, data_frame, x='', y=[], title='', xlab='',
#                  ylab='', y_axis_type='auto', legend=[]):
#         """Draw some parameter varying and the result on the objective value.
#
#         Args:
#             data_frame: a data frame with keys contained in x and y.
#             x: Key in data-frame to use as x-axis.
#             y: Keys in data-frame to plot on y-axis.
#             title: Title for a plot.
#             xlab: Label for x-axis.
#             ylab: Label for y-axis.
#             y_axis_type: Specify "log" to pass logarithmic scale.
#             legend : List of strings matching y.
#
#         Returns:
#             Plot object on success.
#
#         Raises:
#             MissingVariablesException: Dependent variable or their data
#                                         not passed.
#             BadDataFrameException: No data-frame was generated for
#                                     the model object.
#         """
#         pass
#
#
# class SensitivityPlot(BokehPlot):
#     def __init__(self, data_frame, x='', y=[], title='', xlab='',
#                     ylab='', y_axis_type='auto', legend=[]):
#         """Need more information.
#         """
#         pass
#
#
# class ResidualPlot(BokehPlot):
#     def __init__(self, data_frame, x='', y=[], title='', xlab='',
#                  ylab='', y_axis_type='auto', legend=[]):
#         """Plot x, some continuous value (e.g: T, P), against Y (% residual value).
#         Is this %-value calculated from variables in the idaes_model_object?
#
#         Args:
#             data_frame: a data frame with keys contained in x and y.
#             x: Key in data-frame to use as x-axis.
#             y: Keys in data-frame to plot on y-axis.
#             title: Title for a plot.
#             xlab: Label for x-axis.
#             ylab: Label for y-axis.
#             y_axis_type: Specify "log" to pass logarithmic scale.
#             legend : List of strings matching y.
#
#         Returns:
#             Plot object on success.
#         Raises:
#             MissingVariablesException: Dependent variable or their data
#                                         not passed.
#             BadDataFrameException: No data-frame was generated for
#                                     the model object.
#         """
#         pass
#
#
# class IsobarPlot(BokehPlot):
#     def __init__(self, data_frame, x='', y=[], title='', xlab='',
#                ylab='', y_axis_type='auto', legend=[]):
#         """Need more information.
#         """
#         pass
#
#
# class StreamTablePlot(BokehPlot):
#     def __init__(self, data_frame, title=''):
#         """Display a table for all names in the idaes_model_object_names
#         indexing rows according to row_start and row_stop.
#
#         Args:
#             data_frame: a data frame with keys contained in x and y.
#             title: Title for a plot.
#
#         Returns:
#             Plot object on success.
#         Raises:
#             MissingVariablesException: Dependent variable or their data
#                                         not passed.
#             BadDataFrameException: No data-frame was generated for
#                                     the model object.
#         """
#         pass
