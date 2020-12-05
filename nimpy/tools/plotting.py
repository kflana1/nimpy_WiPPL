import numpy as np
from math import atan2, degrees
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import h5py

'''
this module contains nimpy plotting routines and helpers
It depends on:
    numpy (as np)
    matplotlib.pyplot (as plt)
'''
def get_lims(data, num=50, frac=0.01):
    ''' handy function to guess levels and cmap for data

    get_lims chooses levels based on the min and max
    of data input and then either a divergent colormap
    (if +/- values exist) or a simple rainbow colormap 
    for one sign only values

    Args:
        data (np.ndarray): the data to use
        num (int, default=50): the number of levels to return
        frac (float, default=0.01): the fraction of (max-min)
            to use to consider the data as having a negative
            or positive value

    Returns:
        (cmap, levels): cmap (str): the colormap name
                        levels (np.ndarray): the contour levels
    '''
    mx = np.nanmax(data)
    mn = np.nanmin(data)
    rg = mx - mn
    ### if the minimum is less than some tiny negative number
    ### and the maximum is greater than some tiny number
    if mn < -frac*rg and mx > frac*rg:
        cmap = 'coolwarm'
        ### want the levels to be symmetric about 0
        edge = max([np.abs(mn), mx]) 
        ### ensures 0 is a level if odd number of levels
        if num %2 == 0:
            num += num
        levels = np.linspace(-edge, edge, num, endpoint=True)
    ### if the maximum magnitude is positive, use rainbow
    elif np.abs(mx) > np.abs(mn):
        cmap = 'rainbow'
        levels = np.linspace(mn, mx, num, endpoint=True)
    ### if the maximum magnitude is positive use reversed
    ### rainbow colormap
    else:
        cmap = 'rainbow_r'
        levels = np.linspace(mn, mx, num, endpoint=True)
    return cmap, levels

def contourf(data, grid, fax=None, cmap_levels=None, title=None):
    ''' simple filled contourf plot

    contourf is a simple wrapper for matplotlib.pyplot.contourf
    it includes a colorbar, variable cmap and variable contour
    levels.

    Args:
        data (np.ndarray ndim=2): data you wish to plot
        grid (tuple of np.ndarray's): the X and Y arrays for the data
        fax (tuple, default=None): the figure and axes to use
            if not provided, this function will create it's own
        cmap_levels (tuple (str, np.ndarray)): the cmap and
            contour levels to use. If not provided get_lims is
            called to get a good cmap and levels
        title (str, default is None): If you want to add a title

    Returns:
        f (matplotlib.pyplot.figure): the figure
        ax (matplotlib.pyplot.axes): the axes

    Note: does not actually create plot window. to do so, call
        matplotlib.pyplot.show()
    '''
    if cmap_levels is None:
        cmap, levels = get_lims(data)
    else:
        cmap, levels = cmap_levels
    if fax is None:
        f, ax = plt.subplots()
    else:
        f, ax = fax
    try: 
        cf = ax.contourf(grid[0], grid[1], data, levels=levels, cmap=cmap, extend='both')
    except ValueError:
        cf = ax.contourf(grid[0], grid[1], data)
    f.colorbar(cf, ax=ax, fraction=0.05)
    ax.set_aspect('equal')
    if title is not None:
        ax.set_title(title)
    return f,ax

def vec_contourf(tor, pol, grid, psi=None, title=None, tor_cmap_levels=None, pol_cmap_levels=None, psi_levels=None):
    ''' simple filled contourf plot for vector fields

    vec_contourf plots the toroidal component and poloidal magnitude
    for a vector field. If psi is included, it also plots streamlines
    over the poloidal magnitude plot.

    Args:
        tor (np.ndarray): toroidal component (f2)
        pol (np.ndarray): poloidal magnitude (np.sqrt(f1**2+f3**2))
        grid (tuple of np.ndarray's): the X and Y arrays for the data
        psi (np.ndarray, optional): the streamfunction
        title (str, optional): a title for the plot
        tor_cmap_levels (tuple, optional): the toroidal
            colormap (str) and contour levels (np.ndarray)
        pol_camp_levels (tuple, optional): the poloidal
            colormap (str) and contour levels (np.ndarray)
        psi_levels (np.ndarray): the contour levels for the 
            streamfunction
        if any of these cmaps and levels aren't provided, get_lims()
        will be called to find best cmap and levels

    Note: does not actually create plot window. to do so, call
        matplotlib.pyplot.show()
    '''
    if psi is not None and psi_levels is None:
        _, psi_levels = get_lims(psi, num=12)

    f, axs = plt.subplots(1, 2, sharex=True, sharey=True)
    contourf(tor, grid, (f, axs[0]), tor_cmap_levels, title='toroidal component') 
    contourf(pol, grid, (f, axs[1]), pol_cmap_levels, title='poloidal magnitude')

    if psi is not None:
        if psi_levels is None:
            _, psi_levels = get_lims(psi, num=12)
        axs[1].contour(grid[0], grid[1], psi, levels=psi_levels, colors='k', linewidths=1.5)
    
    if title is not None:
        f.suptitle(title)

def plot_WiPAL_vec(nimobj, fieldname, mult=1., fax=None, levels=None, units=None, psiname=None,cbar=True, cbfrac=0.05):

    comps = nimobj.comps
    R, Z = nimobj.grid()
    fieldobj = getattr(nimobj, fieldname)
    tor = mult * getattr(getattr(fieldobj, comps[1]), 'field')
    p1 = getattr(getattr(fieldobj, comps[0]), 'field')
    p3 = getattr(getattr(fieldobj, comps[2]), 'field')
    pol = mult * np.sqrt(p1**2+p3**2)
    
    if psiname is None:
        if hasattr(fieldobj, 'Psi'):
            psi = getattr(getattr(fieldobj, 'Psi'), 'field')
            _, psi_levels = get_lims(psi, num=12)
        else:
            psi = None
    else:
        psi = getattr(getattr(getattr(nimobj, psiname), 'Psi'), 'field')
        _, psi_levels = get_lims(psi, num=10)

    if nimobj.time == 0:
        psi = None

    cmap, levs = get_lims(np.concatenate((tor, pol)), num=50)
    if levels is None:
        levels = levs

    if fax is None:
        f, ax = plt.subplots()
    else:
        f, ax = fax

    cf_tor = ax.contourf(R, Z, tor, levels=levels, cmap=cmap, extend='both')
    if cbar:
        cb = f.colorbar(cf_tor, ax=ax, fraction=cbfrac, ticks=np.linspace(levels[0], levels[-1], 7, endpoint=True))
        if units is not None:
            cb.set_label(units, fontsize=16)
    cf_pol = ax.contourf(-R, Z, pol, levels=levels, cmap=cmap, extend='both')
    if psi is not None:
        c_psi = ax.contour(-R, Z, psi, levels=psi_levels, colors='k', linewidths=1)
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_ylim([-1.55, 1.55])
    ax.set_xlim([-1.55, 1.55])
    ax.add_artist(plt.Circle((0,0), 1.5, color='gray', lw=5, fill=False, zorder=15))
    ax.plot([0]*2,[-1.5,1.5],'gray',lw=5)
    ax.text(-1.5, 1.5, r'${0}_{{pol}}$'.format(fieldname), fontsize=22)
    ax.text(1.0, 1.5, r'${0}_{{\phi}}$'.format(fieldname), fontsize=22)

def plot_WiPAL_VBJ(nimobj, title=None, saveit=None, tite=True, block=True, units=[r'km/s', r'G', r'kA/m$^{2}$'], mults=[1.e-3, 1.e4, 1.e-3], levels=[None]*3, time=True):
    
    nimobj.add_VectorField('Rho', nimobj.n.field * nimobj.V.R.field, nimobj.n.field * nimobj.V.Phi.field, nimobj.n.field * nimobj.V.Z.field)

    f, axs = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(16, 9))
    
    plot_WiPAL_vec(nimobj, 'V', mult=mults[0], fax=(f,axs[0]), units=units[0], levels=levels[0], psiname='Rho')
    plot_WiPAL_vec(nimobj, 'B', mult=mults[1], fax=(f,axs[1]), units=units[1], levels=levels[1])
    plot_WiPAL_vec(nimobj, 'J', mult=mults[2], fax=(f,axs[2]), units=units[2], levels=levels[2])
    
    if time:
        axs[0].text(-1.5, -2, 't={0:.2f} ms'.format(nimobj.time*1.e3), fontsize=22)
    if title is not None:
        f.suptitle(title, fontsize=24, y=0.88)

    if tite:
        plt.tight_layout()
    
    if saveit is None:
        plt.show(block=block)
    else:
        plt.savefig(saveit, bbox_inches='tight')
        plt.close(f)

def plot_PCX_vec(nimobj, fieldname, mult=1., fax=None, levels=None, units=None, psiname=None,cbar=True, cbfrac=0.03):

    comps = nimobj.comps
    R, Z = nimobj.grid()
    fieldobj = getattr(nimobj, fieldname)
    tor = mult * getattr(getattr(fieldobj, comps[1]), 'field')
    p1 = getattr(getattr(fieldobj, comps[0]), 'field')
    p3 = getattr(getattr(fieldobj, comps[2]), 'field')
    pol = mult * np.sqrt(p1**2+p3**2)
    
    if psiname is None:
        if hasattr(fieldobj, 'Psi'):
            psi = getattr(getattr(fieldobj, 'Psi'), 'field')
            _, psi_levels = get_lims(psi, num=12)
        else:
            psi = None
    else:
        psi = getattr(getattr(getattr(nimobj, psiname), 'Psi'), 'field')
        _, psi_levels = get_lims(psi, num=10)

    if nimobj.time == 0:
        psi = None

    cmap, levs = get_lims(np.concatenate((tor, pol)), num=50)
    if levels is None:
        levels = levs

    if fax is None:
        f, ax = plt.subplots()
    else:
        f, ax = fax

    cf_tor = ax.contourf(R, Z, tor, levels=levels, cmap=cmap, extend='both')
    if cbar:
        cb = f.colorbar(cf_tor, ax=ax, fraction=cbfrac, ticks=np.linspace(levels[0], levels[-1], 7, endpoint=True))
        if units is not None:
            cb.set_label(units, fontsize=16)
    cf_pol = ax.contourf(-R, Z, pol, levels=levels, cmap=cmap, extend='both')
    if psi is not None:
        c_psi = ax.contour(-R, Z, psi, levels=psi_levels, colors='k', linewidths=1)
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_ylim([-0.6, 0.6])
    ax.set_xlim([-0.6, 0.6])
    ax.add_patch(Rectangle((-0.5,-0.5), 1, 1, color='gray', lw=5, fill=False, zorder=15))
    ax.plot([0]*2,[-0.5,0.5],lw=5,color='gray')
    ax.text(-0.25, 0.55, r'${0}_{{pol}}$'.format(fieldname), fontsize=22, ha='center')
    ax.text(0.25, 0.55, r'${0}_{{\phi}}$'.format(fieldname), fontsize=22, ha='center')

def plot_PCX_VBJ(nimobj, title=None, saveit=None, tite=True, block=True, units=[r'm/s', r'G', r'kA/m$^{2}$'], mults=[1., 1.e4, 1.e-3], levels=[None]*3, time=True):
    
    nimobj.add_VectorField('Rho', nimobj.n.field * nimobj.V.R.field, nimobj.n.field * nimobj.V.Phi.field, nimobj.n.field * nimobj.V.Z.field)

    f, axs = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(16, 9))
    
    plot_PCX_vec(nimobj, 'V', mult=mults[0], fax=(f,axs[0]), units=units[0], levels=levels[0], psiname='Rho')
    plot_PCX_vec(nimobj, 'B', mult=mults[1], fax=(f,axs[1]), units=units[1], levels=levels[1])
    plot_PCX_vec(nimobj, 'J', mult=mults[2], fax=(f,axs[2]), units=units[2], levels=levels[2])
    
    if time:
        axs[0].text(-0.5, -0.7, 't={0:.2f} ms'.format(nimobj.time*1.e3), fontsize=22)
    if title is not None:
        f.suptitle(title, fontsize=24, y=0.88)

    if tite:
        plt.tight_layout()
    
    if saveit is None:
        plt.show(block=block)
    else:
        plt.savefig(saveit, bbox_inches='tight')
        plt.close(f)


#Label line with line2D label data
def labelLine(line,x,label=None,align=True,**kwargs):

    ax = line.get_axes()
    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if (x < xdata[0]) or (x > xdata[-1]):
        print('x label location is outside data range!')
        return

    #Find corresponding y co-ordinate and angle of the line
    ip = 1
    for i in range(len(xdata)):
        if x < xdata[i]:
            ip = i
            break

    y = ydata[ip-1] + (ydata[ip]-ydata[ip-1])*(x-xdata[ip-1])/(xdata[ip]-xdata[ip-1])

    if not label:
        label = line.get_label()

    if align:
        #Compute the slope
        dx = xdata[ip] - xdata[ip-1]
        dy = ydata[ip] - ydata[ip-1]
        ang = degrees(atan2(dy,dx))
     
	#Transform to screen co-ordinates
        pt = np.array([x,y]).reshape((1,2))
        trans_angle = ax.transData.transform_angles(np.array((ang,)),pt)[0]

    else:
        trans_angle = 0

    #Set a bunch of keyword arguments
    if 'color' not in kwargs:
        kwargs['color'] = line.get_color()

    if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
        kwargs['ha'] = 'center'

    if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
        kwargs['va'] = 'center'

    if 'backgroundcolor' not in kwargs:
        kwargs['backgroundcolor'] = ax.get_axis_bgcolor()

    if 'clip_on' not in kwargs:
        kwargs['clip_on'] = True

    if 'zorder' not in kwargs:
        kwargs['zorder'] = 2.5

    ax.text(x,y,label,rotation=trans_angle,**kwargs)

def labelLines(lines,align=True,xvals=None,**kwargs):

    ax = lines[0].get_axes()
    labLines = []
    labels = []

    #Take only the lines which have labels other than the default ones
    for line in lines:
        label = line.get_label()
        if "_line" not in label:
            labLines.append(line)
            labels.append(label)

    if xvals is None:
        xmin,xmax = ax.get_xlim()
        xvals = np.linspace(xmin,xmax,len(labLines)+2)[1:-1]

    for line,x,label in zip(labLines,xvals,labels):
        labelLine(line,x,label,align,**kwargs)

