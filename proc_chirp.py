from pyspecdata import *
from scipy.optimize import minimize,leastsq
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
from scipy.signal import tukey
from pylab import *
import matplotlib.pyplot as plt
do_slice = False # slice frequencies and downsample -- in my hands, seems to decrease the quality of the fit 
standard_cost = False # use the abs real to determine t=0 for the blip -- this actually doesn't seem to work, so just use the max
show_transfer_func = False # show the transfer function -- will be especially useful for processing non-square shapes
logger = init_logging('info')
# 2to3 JF 1/31

with figlist_var() as fl:    
    # {{{ load data, set units, show raw data
    for searchstr,exp_type,nodename,postproc,corrected_volt,y_lim in [
            ('201113_triwave_control_2','ODNP_NMR_comp','capture1','chirp',True,1.25)
            ]:
        a = find_file(searchstr, exp_type=exp_type, expno=nodename,
                postproc=postproc, lookup=postproc_dict)
        fl.next('Raw signal')
        fl.plot(a['ch',0], alpha=0.5, label='control')    
        fl.plot(a['ch',1], alpha=0.5, label='reflection')
        #}}}
        #{{{FT and plot control/reflection
        a.ft('t',shift=True)
        fl.next('freq domain for RM probe with control')
        fl.plot(a['ch',0],alpha=0.5,label='control')
        fl.plot(a['ch',1], alpha=0.5, label='reflection')
        #}}}
        #{{{throw out low frequencies
        a_control = a['ch',0]['t':(0,None)]
        a_refl = a['ch',1]['t':(0,None)]
        #}}}
        b = a.C
        #{{{divide reflection by control
        a = a_refl/a_control
        fl.next('when divided by control for %s'%searchstr)
        fl.plot(abs(a),label='shorting cap')
        #}}}
        ylim(0,y_lim)
        b.reorder('ch',first=False)
        b.ift('t')
        b *= 2
        ranges = abs(b)['ch',0].contiguous(lambda x: x > 0.09*x.data.max())
        ranges = ranges[0,:].tolist()
        logger.info(strm('Slicing chirp for',searchstr,'from',ranges[0]*1e6,'to',
            ranges[1]*1e6,'us'))
        b = b['t':tuple(ranges)]
        b.setaxis('t',lambda x: x-b.getaxis('t')[0])
        b.setaxis('t',lambda x: 25e6-x*25e6/4096e-8)
        fl.next('$S_{11}$ : analytic amplitude')
        ratio = b['ch',1]/b['ch',0]
        ratio.name('Reflection')
        plot_params = dict(alpha=0.8,
                markersize=2)
        fl.plot(2*abs(ratio),'.',**plot_params)
        axhline(0.425,color='gray',alpha=0.5)
        fl.next('$S_{11}$ : phase')
        fl.plot((ratio).angle/pi,'.',c='blue',**plot_params)

