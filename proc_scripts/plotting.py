from pyspecdata import *
from sympy import symbols
def expand_limits(thisrange,s,axis='t2'):
    """" Used to expand limits of a range (typically used for slicing) by 3X

    Parameters
    ----------
    thisrange: tuple of 2 floats
        slice range you are interested in extending
    s: nddata
        the data you are planning on slicing
    axis: str
        name of the axis along which the slice range will be applied

    Returns
    -------
    retval: tuple of 2 floats
        the expanded range
    """
    thisrange = list(thisrange)
    full_range = s.getaxis(axis)[r_[0,-1]]
    retval = array([thisrange[j] if thisrange[j] is not
            None else full_range[j] for j in [0,-1]])
    print(repr(retval))
    m = mean(retval)
    s = retval-m
    retval = 3*s+m
    sgn = [-1,1] # greater than or less than
    return tuple(full_range[j] if
            retval[j]*sgn[j] > full_range[j]*sgn[j]
            else
            retval[j] for j in range(2))

def draw_limits(thisrange,s):
    """
    determines the range of the t2 axis and pairs it with the range given
    returns the x axis limits for plotting

    Parameters
    ----------
    thisrange: tuple of 2 floats
        slice range you are interested in extending
    s: nddata
        the data you are planning on slicing
    Returns
    -------
    the limits for the t2 axis for plotting purposes
    """
    full_range = s.getaxis('t2')[r_[0,-1]]
    dw = diff(s.getaxis('t2')[:2]).item()
    print("I find the full range to be",full_range)
    sgn = [-1,1] # add or subtract
    pairs = [[thisrange[j],full_range[j]+0.5*dw*sgn[j]] for j in range(2)]
    pairs[0] = pairs[0][::-1] # flip first
    my_xlim = gca().get_xlim() # following messes w/ xlim for some reason
    for j in pairs:
        if None not in j:
            print("drawing a vspan at",j)
            axvspan(j[0],j[1],color='w',alpha=0.5,linewidth=0)
    gca().set_xlim(my_xlim)

class fl_mod(figlist_var):
    """
    Used to create an image for comparison where two images or plots are 
    side by side in the same window. Takes characteristics from figlist_var
    Parameters
    ==========
    plotname:   string with name of plot 
    s:          nddata being analyzed
    thisrange:  range along x axis to be analyzed

    Returns
    =======
    plot image side by side with the cropped log 

    """
    def real_imag(self,plotname,s):
        thisfig,(ax1,ax2) = subplots(1,2)
        self.next(plotname, fig=thisfig)
        sca(ax1)
        self.image(s.real)
        title('real')
        my_clim = gci().get_clim()
        sca(ax2)
        self.image(s)
        gci().set_clim(my_clim) #to match real
        title('imaginary')
        return
    def side_by_side(self,plotname,s,thisrange):
        """a bit of a hack to get the two subplots into
        the figure list -- also a good test for objective
        figure list -- for each slice out 3x thisrange, and then
        show the lines for thisrange"""
        thisfig,(ax1,ax2) = subplots(1,2)
        self.next(plotname, fig=thisfig)
        sca(ax1)
        forplot = s['t2':expand_limits(thisrange,s)]
        self.image(forplot.C.setaxis('power','#').set_units('power','scan #'))
        draw_limits(thisrange,forplot)
        sca(ax2)
        self.image(forplot.C.cropped_log().C.setaxis(
'power','#').set_units('power','scan #'))
        draw_limits(thisrange,forplot)
        title('cropped log')
        return
    def plot_curve(fl, f, name, guess=None):
        """Plot the data with fit curve and fit equation.

        Parameters
        ----------
        f: fitdata
            data (a fitting instance), on which `f.fit()` has already been run
        guess: None or nddata
            The result of `s.settoguess();s.eval(100)` where 100 can be
            any integer.
            Used to display the guess on the plot as well.
        name: str
            the name of the plot
        """
        fl.next(name)
        fl.plot(f, 'o', label=f.name())
        fl.plot(f.eval(100), label='%s fit'%f.name())
        if guess is not None:
            fl.plot(guess, '-', label='initial guess')
        text(0.75, 0.25, f.latex(), transform=gca().transAxes, size='large',
                horizontalalignment='center',color='k')
