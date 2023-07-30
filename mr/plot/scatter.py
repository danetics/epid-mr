import matplotlib.pyplot as plt
import pandas as pd

def plot_points(x, y, xerr=None, yerr=None, colour='lightgrey', marker='.', size=3, width=.5):
    '''
    Combine processes of plotting points and their error bars into one convenience function
    '''
    plt.scatter(x, y, marker=marker, c=colour, s=size)
    if not ((xerr is None) & (yerr is None)):
        plt.errorbar(x, y, xerr=se_x, yerr=se_y, fmt='none', ecolor=colour, elinewidth=width)
    
def draw_plot(x: pd.Series,
                 y: pd.Series,
                 xlabel: str, 
                 ylabel: str,
                 res: dict,
                 highlight=None
                 reglines = ['ivw','egger','wm'],
                 xerr=None,
                 yerr=None):
    '''
    Multipurpose scatter plot
    To be used for both standard (filtered) MR plots and the radial plot
    '''
    # Plot variants (.=matches, square=proxies)
    data['plot_symbol'] = np.where((data['rsid_x']==data['rsid_y']), '.', 'square')
    plot_points(x, y, xerr, yerr, marker=data['plot_symbol'])
    if highlight is not None:
        plot_points(x[highlight], y[highlight], xerr[highlight], yerr[highlight], colour='indianred')
    
    # Overlay specified model lines
    for model in reglines:
        plt.plot(x, res[model][0].fittedvalues, linewidth=1, label=res[model][1], color=res[model][2])
    

    
    apply_standard_format(xlabel, ylabel)



    