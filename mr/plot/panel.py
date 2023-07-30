from . import scatter

def mr_panel(data_xy, data_yx, res_xy, res_yx):
    '''
    Wrapper to get 2x2 plot panel
    '''
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # Top right (ax1) is the standard MR scatter with Steiger failures flagged
    scatter.draw_plot()