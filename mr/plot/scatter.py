#import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import Union
from . import format, labels, text
from .reference import fonts, plotcolours, modelnames

def plot_points(ax, x, y, xerr, yerr, colour='lightgrey', marker='.', size=3, width=.5):
    '''
    Combine processes of plotting points and their error bars into one convenience function
    '''
    ax.scatter(x, y, marker=marker, c=colour, s=size)
    if not ((xerr is None) & (yerr is None)):
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='none', ecolor=colour, elinewidth=width)
    
def draw_plot(data: pd.DataFrame,
              xcol: str,
              ycol: str,
              xerr: Union[str, None],
              yerr: Union[str, None],
              highlight: str,
              res: dict,
              xlabel: str,
              ylabel: str,
              title: str,
              ax, 
              reglines = ['ivw','egger','wm', 'ivw_steig_filtered']):
    '''
    Multipurpose scatter plot
    To be used for both standard (filtered) MR plots and the radial plot
    '''
    proxy = np.where(data['rsid_x']!=data['rsid_y'], True, False)
    plot_points(ax, data[~proxy][xcol], data[~proxy][ycol], 
                (data[~proxy][xerr] if xerr is not None else None), 
                (data[~proxy][yerr] if yerr is not None else None))
    plot_points(ax, data[proxy][xcol], data[proxy][ycol], 
                (data[proxy][xerr] if xerr is not None else None), 
                (data[proxy][yerr] if yerr is not None else None), 
                marker='s')
    if highlight is not None:
        plot_points(ax, data[((~proxy) & (data[highlight]))][xcol], data[((~proxy) & (data[highlight]))][ycol],
                    (data[((~proxy) & (data[highlight]))][xerr] if xerr is not None else None), 
                    (data[((~proxy) & (data[highlight]))][yerr] if yerr is not None else None), 
                    colour='indianred')
        plot_points(ax, data[((proxy) & (data[highlight]))][xcol], data[((proxy) & (data[highlight]))][ycol],
                    (data[((proxy) & (data[highlight]))][xerr] if xerr is not None else None), 
                    (data[((proxy) & (data[highlight]))][yerr] if yerr is not None else None), 
                    colour='indianred', marker='s')
    
    # Overlay specified model lines
    for model in [x for x in reglines if not ((x=='wm') | ('filtered' in x))]:
        ax.plot(data[xcol], res[model].fittedvalues, linewidth=1, 
                label=labels.generate_legend_label(res[model], model), 
                color=plotcolours[model])
    for model in [x for x in reglines if 'filtered' in x]:
        ax.plot(data[xcol][~data[highlight]], res[model].fittedvalues, linewidth=1, 
                label=labels.generate_legend_label(res[model], model), 
                color=plotcolours[model])
    if 'wm' in reglines:
        ax.plot(data[xcol], np.abs(res['wm']['beta']*data[ycol].sort_values()), linewidth=1, 
                label=labels.generate_legend_label_wm(res['wm']), 
                color=plotcolours[model])

    format.style_axes(xlabel, ylabel, title, ax)
    format.add_horizontal_zeroline(ax)
    format.add_legend(ax)
    format.add_textbox(text.get_variant_counts(data, highlight), text.get_qstats(data, reglines), ax)
    



    