from tempfile import NamedTemporaryFile
from bokeh.plotting import figure, show
from bokeh.io import output_file
import numpy as np
from trendfilter import bokeh_theme


def plot_model(result, title='', file=None, show_base=False, show_extrap=False,
               show_plot=True, plot_x_min=0, plot_x_max=None,
               plot_y_min=0, plot_y_max=None,
               extrap_min=0, extrap_max=40, extrap_stride=1):

    if file is None:
        file = NamedTemporaryFile().name+'.html'

    output_file(file)

    x = result['x']
    y = result['y']

    if plot_x_max is None:
        plot_x_max = max(x)
        if show_extrap:
            plot_x_max += extrap_max

    if plot_y_max is None:
        plot_y_max = max(y) * 1.05

    plot = figure(title=title, width=900, height=600,
                  y_range=(plot_y_min, plot_y_max),
                  x_range=(plot_x_min, plot_x_max))

    plot.circle(x, y, legend_label='data')
    plot.line(x, y)

    plot.line(x, result['y_fit'], color='red', legend_label='model')

    if show_extrap:
        # over-plot the function, showing the extrapolation too
        f = result['function']
        x_min = max(x) + extrap_min
        x_max = max(x) + extrap_max
        xx = np.arange(x_min, x_max, extrap_stride)
        plot.line(xx, f(xx), color='green', legend_label='model extrapolation')

    if show_base:
        f_base = result['function_base']
        x_max = max(x)
        if extrap_max:
            x_max += extrap_max

        xxx = np.arange(x.min(), x_max)
        plot.line(xxx, f_base(xxx), color='black', legend_label='base model')

    if show_plot:
        show(plot)

    return plot
