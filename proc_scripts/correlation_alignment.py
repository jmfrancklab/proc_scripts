from pylab import*
from pyspecdata import *
from scipy.optimize import leastsq, minimize, basinhopping
import numpy as np
from matplotlib.ticker import FuncFormatter
@FuncFormatter
def to_percent(y, position):
    s = '%.2g'%(100 * y)
    if rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'
def correl_align(s, align_phases=False,indirect_dim='indirect',fig_title='correlation alignment',maxiter=50,tol=1e-4,fl=None):
    N = ndshape(s)[indirect_dim]
    sig_energy = (abs(s)**2).data.sum().item() / N
    if fl is not None:
        fl.next('before correlation\nsig. energy=%g'%sig_energy + fig_title)
        fl.image(s)
    energy_diff = 1.
    i = 0
    energy_vals = []
    this_E = (abs(s.C.sum(indirect_dim))**2).data.sum().item() / N**2
    energy_vals.append(this_E / sig_energy)
    last_E = None
    for my_iter in range(maxiter):
        i += 1
        print("*** *** ***")
        print("CORRELATION ALIGNMENT ITERATION NO. ",i)
        print("*** *** ***")
        if align_phases:
            ph0 = s.C.sum('t2')
            ph0 /= abs(ph0)
            s /= ph0
        s.ift('t2')
        correl = s * 0
        for q in range(1,ndshape(s)[indirect_dim]):
            temp = s.C
            temp = nddata(roll(temp.data.conj(),q,axis=0),[indirect_dim,'t2']) 
            correl +=s * temp
        #correl.setaxis('t2',s.getaxis('t2'))
        #correl.setaxis(indirect_dim,s.getaxis(indirect_dim))
        correl.ft_clear_startpoints('t2')
        correl.ft('t2', shift=True, pad=2**14)
        if fl is not None:
            fl.next('Correlation function')
            fl.image(correl)
        f_shift = correl.run(real).argmax('t2')
        #f_shift = nddata(f_shift.data,[indirect_dim])
        s *= np.exp(-1j*2*pi*f_shift*s.fromaxis('t2'))
        s.ft('t2')
        if align_phases:
            s *= ph0
        this_E = (abs(s.C.sum(indirect_dim))**2).data.sum().item() / N**2
        energy_vals.append(this_E / sig_energy)
        print('averaged signal energy (per vd):', this_E)
        if last_E is not None:
            energy_diff = (this_E - last_E)/sig_energy
            print(energy_diff)
            if abs(energy_diff) < tol:
                break
        last_E = this_E
    if fl is not None: 
        fl.next('correlation convergence')
        fl.plot(array(energy_vals),'x')
        gca().yaxis.set_major_formatter(to_percent)
    #sig_energy = (abs(s)**2).data.mean().item()
    print("*** *** ***")
    print("after",sig_energy)
    print("*** *** ***")
    if fl is not None:
        fl.next('After correlation\nph0 restored sig energy=%g'%sig_energy + fig_title)
        fl.image(s)
    return s, energy_vals    

