'''
'''

from .reference import fonts

def style_axes(xlabel, ylabel, title, ax):
    ax.tick_params(labelsize=6)
    ax.set_xlabel(xlabel, fontdict=fonts['normal'])
    ax.set_ylabel(ylabel, fontdict=fonts['normal'])
    ax.set_title(title, fontdict=fonts['bold'], loc='left')

def add_horizontal_zeroline(ax):
    ax.axhline(y=0, linestyle='dotted', linewidth=.5, color='grey')
                   
def add_legend(ax):
    ax.legend(fontsize=5, labelcolor='black', labelspacing=.3, frameon=False)

def add_textbox(counts, qstats, ax):
    ax.text(0.01, 0.99, '\n'.join([counts, qstats]), transform=ax.transAxes, ha='left', va='top', fontdict=fonts['normal'], linespacing=1)