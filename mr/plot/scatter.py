import matplotlib.pyplot as plt
import pandas as pd
from . import format, labels, text
from .reference import fonts, plotcolours, modelnames

def plot_points(x, y, ax, xerr=None, yerr=None, colour='lightgrey', marker='.', size=3, width=.5):
    '''
    Combine processes of plotting points and their error bars into one convenience function
    '''
    ax.scatter(x, y, marker=marker, c=colour, s=size)
    if not ((xerr is None) & (yerr is None)):
        ax.errorbar(x, y, xerr=se_x, yerr=se_y, fmt='none', ecolor=colour, elinewidth=width)
    
def draw_plot(x: pd.Series,
                 y: pd.Series,
                 xlabel: str, 
                 ylabel: str,
                 res: dict,
                 ax: matplotlib.axes._subplots.AxesSubplot
                 highlight=None
                 reglines = ['ivw','egger','wm', 'ivw_steig_filtered'],
                 xerr=None,
                 yerr=None):
    '''
    Multipurpose scatter plot
    To be used for both standard (filtered) MR plots and the radial plot
    '''
    # Plot variants (.=matches, square=proxies)
    data['plot_symbol'] = np.where((data['rsid_x']==data['rsid_y']), '.', 'square')
    plot_points(x, y, xerr, yerr, marker=data['plot_symbol'], ax=ax)
    if highlight is not None:
        plot_points(x[highlight], y[highlight], xerr[highlight], yerr[highlight], colour='indianred', ax=ax)
    
    # Overlay specified model lines
    for model in reglines:
        ax.plot(x, res[model][0].fittedvalues, linewidth=1, 
                label=labels.generate_lenend_label(res[model], model), 
                color=plotcolours[model])

    format.style_axes(xlabel, ylabel, ax)
    format.add_horizontal_zeroline(ax)
    format.add_legend(ax)
    format.add_textbox(text.get_variant_counts(data, highlight), text.get_qstats(data, reglines), ax)



    