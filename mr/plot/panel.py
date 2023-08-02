import os
from . import scatter
import matplotlib.pyplot as plt

def mr_panel(dataxy, datayx, resxy, resyx, outdir='.', outname=None):
    '''
    Wrapper to get 2x2 plot panel
    '''
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12,7.5))

    # Top left (ax1) is the standard MR scatter with Steiger failures flagged
    # Plus main model and IWW (steiger-filtered) reglines
    scatter.draw_plot(dataxy, xcol='beta_x', ycol='beta_y', xerr='se_x', yerr='se_y', highlight='steiger_fail', res=resxy,
                      xlabel=dataxy.attrs['xname'], 
                      ylabel=dataxy.attrs['yname'],
                      title='Main Models',
                      ax = ax1)
    # Top right (ax2) is the radial plot with radial failures flagged and radial-scale reglines
    scatter.draw_plot(dataxy, xcol='ratio_inv_se', ycol='ratio_z', xerr=None, yerr=None, highlight='radial_fail', res=resxy,
                      xlabel='Inverse Ratio SE', 
                      ylabel='Ratio Z-score',
                      title='Radial Plot',
                      reglines=['radial','egger_radial'],
                      ax=ax2)
    # Bottom left (ax3) is standard scatter with radial failures flagged and radial-filtered reglines
    scatter.draw_plot(dataxy, xcol='beta_x', ycol='beta_y', xerr='se_x', yerr='se_y', highlight='radial_fail', res=resxy,
                      xlabel=dataxy.attrs['xname'], 
                      ylabel=dataxy.attrs['yname'],
                      title='Heterogeneity-filtered Models',
                      reglines=['ivw_radial_filtered','egger_radial_filtered'],
                      ax=ax3)
    # Bottom right (ax4) is the reverse MR standard scatter.  As for ax1 but inverting x and y
    scatter.draw_plot(datayx, xcol='beta_x', ycol='beta_y', xerr='se_x', yerr='se_y', highlight='steiger_fail', res=resyx,
                      xlabel=datayx.attrs['xname'], 
                      ylabel=datayx.attrs['yname'],
                      title='Reversed Main Models',
                      ax = ax4)
    plt.tight_layout()
    outname = f"{dataxy.attrs['xname'].lower().replace(':','.')}.{dataxy.attrs['yname'].lower().replace(':','.')}.unimr.plots" if outname is None else outname
    plt.savefig(os.path.join(outdir, f'{outname}.png'), bbox_inches='tight')
    