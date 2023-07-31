from . import scatter

def mr_panel(dataxy, datayx, resxy, resyx):
    '''
    Wrapper to get 2x2 plot panel
    '''
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    # Top left (ax1) is the standard MR scatter with Steiger failures flagged
    # Plus main model and IWW (steiger-filtered) reglines
    scatter.draw_plot(x=dataxy['beta_x'], y=dataxy['beta_y'], 
                      xerr=dataxy['se_x'], yerr=dataxy['se_y'], 
                      xlabel=dataxy.attrs['xname'], ylabel=dataxy.attrs['yname'],
                      res=resxy, # default reglines are set for this plot 
                      highlight=dataxy['steiger_fail'],
                      axes=ax1)
    # Top right (ax1) is the radial plot with radial failures flagged and radial-scale reglines
    scatter.draw_plot(x=1/dataxy['ratio_se'], y=dataxy['ratio']/dataxy['ratio_se'], 
                      xlabel='1/ratio_se', ylabel='ratio_z',
                      res=resxy, reglines=['radial','egger_radial'],
                      highlight=(dataxy['radial_fail_coch'] | dataxy['radial_fail_ruck']),
                      axes=ax2)
    # Bottom left is standard scatter with radial failures flagged and radial-filtered reglines
    scatter.draw_plot(x=dataxy['beta_x'], y=dataxy['beta_y'], 
                      xerr=dataxy['se_x'], yerr=dataxy['se_y'], 
                      xlabel=dataxy.attrs['xname'], ylabel=dataxy.attrs['yname'],
                      res=resxy, reglines=['ivw_radial_filtered','egger_radial_filtered'],
                      highlight=(dataxy['radial_fail_coch'] | dataxy['radial_fail_ruck']),
                      axes=ax3)
    # Bottom right is the reverse MR standard scatter.  As for ax1 but inverting x and y
    scatter.draw_plot(x=datayx['beta_x'], y=datayx['beta_y'], 
                      xerr=datayx['se_x'], yerr=datayx['se_y'], 
                      xlabel=datayx.attrs['xname'], ylabel=datayx.attrs['yname'],
                      res=resyx, # default reglines are set for this plot 
                      highlight=datayx['steiger_fail'],
                      axes=ax4)
    plt.show()
    