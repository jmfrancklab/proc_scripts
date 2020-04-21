from pyspecdata import *
from utility import dBm2power
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from hermitian_function_test import hermitian_function_test, zeroth_order_ph
from sympy import symbols
rcParams["savefig.transparent"] = True
def expand_limits(thisrange,d):
    thisrange = list(thisrange)
    full_range = d.getaxis('t2')[r_[0,-1]]
    retval = array([thisrange[j] if thisrange[j] is not
            None else full_range[j] for j in [0,-1]])
    print(repr(retval))
    m = mean(retval)
    d = m-retval
    retval = 3*d+m
    sgn = [-1,1] # greater than or less than
    return tuple(full_range[j] if
            retval[j]*sgn[j] > full_range[j]*sgn[j]
            else
            retval[j] for j in range(2))
def draw_limits(thisrange,d):
    full_range = d.getaxis('t2')[r_[0,-1]]
    dw = diff(d.getaxis('t2')[:2]).item()
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
class fl_mod (figlist_var):
    def side_by_side(self,plotname,s,thisrange):
        """a bit of a hack to get the two subplots into
        the figure list -- also a good test for objective
        figure list -- for each slice out 3x thisrange, and then
        show the lines for thisrange"""
        thisfig,(ax1,ax2) = subplots(1,2)
        fl.next(plotname, fig=thisfig)
        sca(ax1)
        forplot = s['t2':expand_limits(thisrange,s)]
        fl.image(forplot)
        print("drawing limits",thisrange)
        draw_limits(thisrange,forplot)
        sca(ax2)
        fl.image(forplot.C.cropped_log())
        draw_limits(thisrange,forplot)
        title('cropped log')
        return
fl = fl_mod()
t2 = symbols('t2')

# notice how I dramatically simplified the following
# going forward, this is how we should deal with filenames
# (that was an earlier commit... then I use the list of
# tuples to specify the ranges for frequency and time)

# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be
for searchstr,freq_range,time_range in [
        ["200306_DNP_lg_probe_w34.*",(-300,300),(None,0.05)]
        ]:
    files = search_filename(searchstr, 'test_equip')
    assert len(files)==1, "I found %d files matching the pattern %s"%(len(files),searchstr)
    dirname, filename = os.path.split(files[0])
    nodename = 'signal'
    s = nddata_hdf5(filename+'/signal',
            directory=dirname)
    prog_power = s.getaxis('power').copy()
    print("programmed powers",prog_power)
    s.setaxis('power',r_[
        0,dBm2power(array(s.get_prop('meter_powers'))+20)]
        ).set_units('power','W')
    print("meter powers",s.get_prop('meter_powers'))
    print("actual powers",s.getaxis('power'))
    print("ratio of actual to programmed power",
                s.getaxis('power')/prog_power)
    nPoints = s.get_prop('acq_params')['nPoints']
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    s.set_units('t','s')
    s.chunk('t',['ph2','ph1','t2'],[2,4,-1])
    s.labels({'ph2':r_[0.,2.]/4,
        'ph1':r_[0.,1.,2.,3.]/4})
    s.reorder(['ph2','ph1'])
    s.ft('t2',shift=True)
    s.ft(['ph1','ph2'])
    fl.next('all data: frequency domain')
    fl.image(s)

    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
            s,freq_range)
    s = s['t2':freq_range]
    s.ift('t2')
    residual,best_shift = hermitian_function_test(s[
        'ph2',-2]['ph1',1])
    fl.next('hermitian test')
    fl.plot(residual)
    s.setaxis('t2',lambda x: x-best_shift)
    s.register_axis({'t2':0}, nearest=False)
    # {{{ implement zeroth-order correction
    # note that it's not going to have only one
    # dimension, b/c we will have at least a power
    # dimension
    ph0 = s['t2':0]['ph2',-2]['ph1',1]

    s /= zeroth_order_ph(ph0, fl=fl)
    if s['t2':0]['ph2',-2]['ph1',1]['power',0].real < 0:
        s *= -1
    # }}}
    fl.side_by_side('time domain (after filtering and phasing)\n$\\rightarrow$ use to adjust time range',
            s,time_range)
    s = s['t2':time_range]
    # {{{ all of this is to check and see if we think
    # we can add the two sides of the echo to increase
    # our SNR -- right now, it doesn't look like it
    fl.next('echo mirror test')
    echo_start = s.getaxis('t2')[0]
    dw = diff(s.getaxis('t2')[r_[0,1]]).item()
    centered_echo = s['t2':(echo_start,-echo_start+dw)]
    plotdata = abs(centered_echo/centered_echo['t2',::-1])
    plotdata[lambda x: x>2] = 2
    fl.image(plotdata)
    # }}}
    fl.next('apodize and zero fill')
    R = 5.0/(time_range[-1]) # assume full decay by end time
    s *= exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=1024)
    fl.image(s)
    # {{{ select FID
    s.ift('t2')
    s = s['ph2',-2]['ph1',1]['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    # }}}
    fl.next('compare highest power to no power')
    idx_maxpower = argmax(s.getaxis('power'))
    fl.plot(s['power',0])
    fl.plot(s['power',idx_maxpower])
    fl.next('full enhancement curve')
    fl.plot(s)
    fl.next('enhancement')
    enhancement = s['t2':(-50,50)].C.sum('t2').real
    enhancement /= enhancement['power',0]
    fl.plot(enhancement['power',:idx_maxpower+1],'ko', human_units=False)
    fl.plot(enhancement['power',idx_maxpower+1:],'ro', human_units=False)

    ylabel('Enhancement')
fl.show();quit()
