import matplotlib.pyplot as plt
import pandas as pd

fonts = {'general': {'family': 'sans', 'color': 'black', 'size': 8}}

def split_proxies(data):
    matches = data[data['rsid_x']==data['rsid_y']]
    proxies = data[data['rsid_x']!=data['rsid_y']]
    return matches, proxies

def plot_points(x, y, colour, marker, size=3, se_x=None, se_y=None):
    plt.scatter(x, y, marker=marker, c=colour, s=3)
    if ((se_x is not None) & (se_y is not None)):
        plt.errorbar(x, y, xerr=se_x, yerr=se_y, fmt='none', ecolor=colour, elinewidth=.5)

def draw_scatter(data, xtrait, ytrait, reg=None, fonts=fonts):
    '''
    Draw MR scatter plot with IVW, WM, Egger regression lines
    '''
    matches, proxies = split_proxies(data)
    plot_points(x=matches['beta_x'], y=matches['beta_y'], se_x=matches['se_x'], se_y=matches['se_y'], colour='cyan', marker='.')
    if not proxies.empty:
        plot_points(x=proxies['beta_x'], y=proxies['beta_y'], se_x=proxies['se_x'], se_y=proxies['se_y'], colour='grey', marker='.')

    # Build a dictionary with long legend names for each regression, including beta and pvalue
    regnames = {k:f'{k}: {v.params['beta_x']} (p={v.pvalues['beta_x']})' for k, v in reg.items()} 
        
    # Plot regression lines from dict of regression objects
    plt.plot(data['beta_x'], reg['ivw'].fittedvalues, linewidth=1, color='blue', name=regnames['ivw'])
    plt.plot(data['beta_x'], reg['wm'].fittedvalues, linewidth=1, color='orange', name=regnames['wm'])
    plt.plot(data['beta_x'], reg['egger'].fittedvalues, linewidth=1, color='red', name=regnames['egger'])
    
    # Plot fixed lines, add legend, format plot
    plt.axhline(y=0, linestyle='dotted', linewidth=.5, color='grey')
    plt.tick_params(labelsize=6)
    plt.xlabel(xtrait, fontdict=fonts['axis_title'])
    plt.ylabel(ytrait, fontdict=fonts['axis_title'])
    #plt.text(f'# variants = {}', fontdict=fonts['general'], 'loc='best')
    plt.legend(fondict=fonts['general'], loc='best')