from pyspecdata import *
from sympy import symbols
import matplotlib.pyplot as plt
import numpy as np
import logging

class fl_mod(figlist_var):
    """
    Extends figlist_var with various new convenience functions.
    """
    def next(self, *arg, **kwargs):
    #    kwargs.update({"figsize": (9, 5.56), "legend": True})
        super().next(*arg, **kwargs)

    def real_imag(self,plotname,s):
        thisfig,(ax1,ax2) = plt.subplots(1,2)
        self.next(plotname, fig=thisfig)
        plt.sca(ax1)
        self.image(s.real)
        title('real')
        my_clim = gci().get_clim()
        plt.sca(ax2)
        self.image(s)
        gci().set_clim(my_clim) #to match real
        title('imaginary')
        return

    def side_by_side(self,plotname,s,thisrange):
        """a bit of a hack to get the two subplots into
        the figure list -- also a good test for objective
        figure list -- for each slice out 3x thisrange, and then
        show the lines for thisrange"""
        thisfig,(ax1,ax2) = plt.subplots(1,2)
        plt.sca(ax1)
        forplot = s['t2':expand_limits(thisrange,s)]
        if 'power' in s.dimlabels:
            self.image(forplot.C.setaxis('power','#').set_units('power','scan #'))
        else:
            self.image(forplot)
        draw_limits(thisrange,forplot)
        plt.sca(ax2)
        if 'power' in s.dimlabels:
            self.image(forplot.C.cropped_log().C.setaxis('power','#').set_units('power','scan #'))
        else:
            self.image(forplot.C.setaxis('vd','#').set_units('vd','scan #'))
        draw_limits(thisrange,forplot)
        title('cropped log')
        #self.next(plotname, fig=figure())
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
        text(0.75, 0.25, f.latex(), transform=plt.gca().transAxes, size='large',
                horizontalalignment='center',color='k')

    def complex_plot(fl, d, label="", show_phase=False, show_real=True,alpha=0.5,linestyle=None,linewidth=3,color='k'):
        colors = []
        for j in range(ndshape(d)["ch"]):
            chlabel = d.getaxis("ch")[j]
            if j==0:
                l = fl.plot(d["ch", j],
                        linestyle=linestyle,
                        linewidth=linewidth,
                        alpha=alpha,
                        label="reflected pulse abs " + label,
                        color=color
                        )
            else:
                l = fl.plot(
                        abs(d["ch",j]),
                        linestyle=linestyle,
                        linewidth=linewidth,
                        alpha=alpha,
                        label="reflected pulse abs" + label,
                        color=color
                        )
            colors.append(l[0].get_color())
            if show_real:
                if j==0:
                    fl.plot(
                        d["ch", j].real,
                        linewidth=1,
                        color=colors[-1],
                        alpha=alpha,
                        label="reflected real " + label,
                        )
                else:
                    fl.plot(abs(d["ch",j].real),
                            linewidth=1,
                            color=colors[-1],
                            alpha=alpha,
                            label="reflected real" + label,
                            )

                if j==0:
                    fl.plot(
                        d["ch", j].imag,
                        "--",
                        linewidth=1,
                        color=colors[-1],
                        alpha=alpha,
                        label="fwd pulse imag" + label,
                        )
                else:
                    fl.plot(
                        d["ch",j].imag,
                        "--",
                        linewidth=1,
                        colors=colors[-1],
                        alpha=alpha,
                        label="reflected imag" + label,
                        )

            fl.grid()
            if show_phase:
                fl.twinx(orig=False)
                if j==0:
                    fl.plot(
                            d["ch", j].angle/2/pi,
                        ".",
                        linewidth=1,
                        color=colors[-1],
                        alpha=0.3,
                        label="reflected angle " + label,
                        )
                else:
                    fl.plot(d["ch",j].angle/2/pi,
                            ".",
                            linewidth=1,
                            color=colors[-1],
                            alpha=0.3,
                            label="reflected" + label,
                            )
                plt.ylabel("phase / cyc", size=10)
                ax2=plt.gca()
                gridandtick(ax2, use_grid=False)
                ax2.grid(False)
                fl.twinx(orig=True)
        return colors

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
    retval = np.array([thisrange[j] if thisrange[j] is not
            None else full_range[j] for j in [0,-1]])
    m = np.mean(retval)
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
    dw = np.diff(s.getaxis('t2')[:2]).item()
    logging.info(strm("I find the full range to be",full_range))
    sgn = [-1,1] # add or subtract
    pairs = [[thisrange[j],full_range[j]+0.5*dw*sgn[j]] for j in range(2)]
    pairs[0] = pairs[0][::-1] # flip first
    my_xlim = plt.gca().get_xlim() # following messes w/ xlim for some reason
    for j in pairs:
        if None not in j:
            logging.info(strm("drawing a vspan at",j))
            plt.axvspan(j[0],j[1],color='w',alpha=0.5,linewidth=0)
    plt.gca().set_xlim(my_xlim)
