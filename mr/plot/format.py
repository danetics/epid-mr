'''
'''

from .reference import fonts

def apply_format(xlabel, ylabel, counts, qstats):
    '''
    Apply standard font formatting, create legend and display text
    '''
    plt.legend(fontsize=6, labelcolor='black', fancybox=True)
    plt.text(0, 1, '\n'.join([counts, qstats], fontdict=fonts['normal'])
    plt.axhline(y=0, linestyle='dotted', linewidth=.5, color='grey')
    plt.tick_params(labelsize=6)
    plt.xlabel(xlabel, fontdict=plotfont)
    plt.ylabel(ylabel, fontdict=plotfont)
    
    