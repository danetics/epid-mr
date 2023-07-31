'''
'''

from .reference import fonts

def style_axes(xlabel, ylabel, ax):
    ax.tick_params(labelsize=6)
    ax.xlabel(xlabel, fontdict=plotfont)
    ax.ylabel(ylabel, fontdict=plotfont)

def add_horizontal_zeroline(ax):
    ax.axhline(y=0, linestyle='dotted', linewidth=.5, color='grey')
                   
def add_legend(ax):
    ax.legend(fontsize=6, labelcolor='black', fancybox=True)

def add_textbox(counts, qstat, ax):
    ax.text(0, 1, '\n'.join([counts, qstats], fontdict=fonts['normal'])
    
    
    