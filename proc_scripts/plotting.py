from pyspecdata import *
from sympy import symbols
def expand_limits(thisrange,s):
    thisrange = list(thisrange)
    full_range = s.getaxis('t2')[r_[0,-1]]
    retval = array([thisrange[j] if thisrange[j] is not
            None else full_range[j] for j in [0,-1]])
    print(repr(retval))
    m = mean(retval)
    s = m-retval
    retval = 3*s+m
    sgn = [-1,1] # greater than or less than
    return tuple(full_range[j] if
            retval[j]*sgn[j] > full_range[j]*sgn[j]
            else
            retval[j] for j in range(2))

def draw_limits(thisrange,s):
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
        self.image(forplot)
        draw_limits(thisrange,forplot)
        sca(ax2)
        self.image(forplot.C.cropped_log())
        draw_limits(thisrange,forplot)
        title('cropped log')
        return

    def recovery(s,rec_curve,f_range):
        fl.next('recovery curve')
        fl.plot(rec_curve,'o')
        f =fitdata(rec_curve)
        M0,Mi,R1,vd = sympy.symbols("M_0 M_inf R_1 indirect",real=True)
        f.functional_form = Mi + (M0-Mi)*sympy.exp(-vd*R1)
        logger.info(strm("Functional form", f.functional_form))
        fl.next('t1 test')
        fl.plot(f, 'o',label=f.name())
        f.fit()
        fl.plot(f.eval(100),label=
                '%s fit'%f.name())
        text(0.75, 0.25, f.latex(), transform=gca().transAxes, size='large',
                horizontalalignment='center',color='k')
        print("output:",f.output())
        print("latex:",f.latex())
        return

    def decay(s, f_range):
        data = s['t2':f_range].sum('t2')
        fl.next('Echo decay')
        fl.plot(data,'o')
        print("starting T2 curve")
        f = fitdata(data.real)
        M0,R2,tE = sympy.symbols("M_0 R_2 tE", real=True)
        f.functional_form = M0*sympy.exp(-tE*R2)
        fl.next('T2 test')
        fl.plot(f,'o',label=f.name())
        f.fit()
        fl.plot(f.eval(100),label='%s fit'%f.name())
        text(0.75,0.25, f.latex(),transform=gca().transAxes, size='large',
                horizontalalignment='center', color= 'k')
        print("output",f.output())
        print("latex",f.latex())
        T2 = 1./f.output('R_2')
        return
    
    def plot_curve(s, rec_curve, f_range, curve):
        if curve == 'recovery':
            recovery(s, rec_curve, f_range)
            return
        if curve == 'decay':
            decay(s, f_range)
            return
