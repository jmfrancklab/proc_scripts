from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
rcParams["savefig.transparent"] = True

fl = fl_mod()
t2 = symbols('t2')

# slice out the FID from the echoes,
# also frequency filtering, in order to generate the
# list of integrals for ODNP
# to use: as a rule of thumb, make the white boxes
# about 2x as far as it looks like they should be

# leave this as a loop, so you can load multiple files
for searchstr,exp_type,nodename,postproc,freq_range,time_range in [
        ["200306_DNP_lg_probe_w34.*", 'test_equip', 'signal',
            'spincore_ODNP_v1', (-300,300), (None,0.05)]
        ]:
    s = find_file(searchstr, exp_type=exp_type, expno=nodename,
            postproc=postproc,
            lookup=postproc_dict)
    
    #{{{Fourier Transforms coherence channels
    s.ft(['ph1','ph2']) 
    #}}}

    #{{{visualize the frequency limits
    fl.next('all data: frequency domain')
    fl.image(s)
    fl.side_by_side('show frequency limits\n$\\rightarrow$ use to adjust freq range',
            s,freq_range)
    #}}}

    #{{{slice out the frequency range along t2 axis
    s = s['t2':freq_range]
    #}}}
    
    #{{{inverse fourier transform into time domain
    s.ift('t2')
    #}}}
    
    #{{{visualize time domain after filtering and phasing 
    logger = init_logging("info")
    logger.info(strm("THIS IS THE SHAPE"))
    logger.info(strm(ndshape(s)))
    s = slice_FID_from_echo(s,0,1)['t2':(None,0.05)]
    fl.side_by_side('time domain (after filtering and phasing)\n$\\rightarrow$ use to adjust time range', s, time_range)
    #}}}
    
    #{{{slices out time range along t2 axis
    s =s['t2':time_range]
    #}}}
    
    #{{{visualize centered
    echo_start = s.getaxis('t2')[0]
    dw = diff(s.getaxis('t2')[r_[0,1]]).item()
    centered_echo = s['t2':(echo_start,-echo_start+dw)]
    plotdata = abs(centered_echo/centered_echo['t2',::-1])
    fl.next('plot data')
    plotdata[lambda x: x>2] = 2
    fl.image(plotdata)
    #}}}

    #{{{slice FID from echo
    logger.info(strm("THIS IS THE SHAPE"))
    logger.info(strm(ndshape(s)))
    s = slice_FID_from_echo(s)['t2':(None,0.05)]
    #}}}
    
    # {{{ align the peaks
    orig = s.C
    fl.next('before alignment')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    old_order = list(s.dimlabels)
    avg = s['ph1',1]['ph2',-2].C.mean_all_but('t2')
    # {{{ try to use the correlation align
    s.ift(['ph1','ph2'])
    s.smoosh(['ph2','ph1','power'],'transient')
    print(ndshape(s))
    s = s['transient',:20]
    print(ndshape(s))
    s.ft('t2')
    avg.ft('t2')
    fl.next('s when smooshed')
    fl.plot(s,human_units=False)
    fl.next('avg')
    fl.plot(avg)
    s.ift('t2')
    avg.ift('t2')
    fl.basename = 'overlay, rainbow order, black=last'
    s = correlation_align(s,avg,color='r',fl=fl)
    s = correlation_align(s,avg,color='y',fl=fl)
    s = correlation_align(s,avg,color='g',fl=fl)
    s = correlation_align(s,avg,color='b',fl=fl)
    s = correlation_align(s,avg,color='m',fl=fl)
    s = correlation_align(s,avg,fl=fl)
    fl.next('s after correlation')
    s.ft('t2')
    fl.plot(s,human_units=False)
    fl.show();quit()
    s.ft(['ph1','ph2'])
    # }}}
    fl.next('after alignment')
    s.reorder(old_order)
    #s.reorder('ph2','ph1','power','indirect','t2')
    s.ft('t2')
    fl.image(s,human_units=False)
    s.ift('t2')
    #}}}

    #{{{redefine time range along t2
    s.set_units('indirect',None)
    fl.side_by_side('time domain (after filtering and phasing)\n$\\rightarrow$ use to adjust time range',
            orig,time_range)
    s = orig['t2':time_range]
    #}}}
    
    #{{{mirror test to test centered data
    fl.next('echo mirror test')
    echo_start = s.getaxis('t2')[0]
    dw = diff(s.getaxis('t2')[r_[0,1]]).item()
    centered_echo = s['t2':(echo_start,-echo_start+dw)]
    plotdata = abs(centered_echo/centered_echo['t2',::-1])
    fl.next('plot data')
    plotdata[lambda x: x>2] = 2
    fl.image(plotdata)
    #}}}
    
    #{{{apodizing and zero fill
    fl.next('apodize and zero fill')
    R = 5.0/(time_range[-1]) # assume full decay by end time
    s *= exp(-s.fromaxis('t2')*R)
    s.ft('t2',pad=1024)
    fl.image(s)
    #}}}
    
    #{{{select coherence channel in time domain
    s.ift('t2')
    s = s['ph2',-2]['ph1',1]['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    #}}}

    #{{{plotting enhancement curve at lowest and highest powers 
    fl.next('compare highest power to no power')
    idx_maxpower = argmax(s.getaxis('power'))
    fl.plot(s['power',0])
    fl.plot(s['power',idx_maxpower])
    #}}}
    
    #{{{plotting full enhancement curve
    fl.next('full enhancement curve')
    fl.plot(s)
    #}}}
    
    #{{{plotting enhancement vs power
    fl.next('enhancement')
    enhancement = s['t2':(-50,50)].sum('t2').real
    enhancement /= enhancement['power',0]
    fl.plot(enhancement['power',:idx_maxpower+1],'ko', human_units=False)
    fl.plot(enhancement['power',idx_maxpower+1:],'ro', human_units=False)
    ylabel('Enhancement')
    #}}}
fl.show();quit()
