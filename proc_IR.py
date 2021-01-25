from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from proc_scripts.fitting import recovery
from sympy import symbols
from matplotlib import *
fl = fl_mod()
t2 = symbols('t2')
logger = init_logging("info")
# {{{ input parameters
filter_bandwidth = 15e3
coh_sel = {'ph1':0,
        'ph2':-1}
coh_err = {'ph1':1,# coherence channels to use for error
        'ph2':r_[0,2,3]}
# }}}
clock_correction = nddata(np.linspace(-3,3,2500),'clock_correction')

for searchstr,exp_type,nodename, postproc in [
        ('210120_OHTEMPO10mM_sol_probe_IR','test_equip','signal','spincore_IR_v1')
        #('w3_201111','test_equip',2,'ab_ir2h')
        #('CTAB_w15_41mM_201124','test_equip',2,'ab_ir2h'),
        #('freeSL_201007','test_equip',5,'ag_IR2H',None)
        #('w8_200731', 'test_equip', 2, 'ag_IR2H',None),
        #('free4AT_201014','test_equip',3,'ag_IR2H',None)
        #('free4AT100mM_201104', 'test_equip',2,'ab_ir2h',None),
        #('ag_oct182019_w0_8','test_equip',3,'ab_ir2h',None)
        #('freeD2O_201104','test_equip',2,'ab_ir2h',None),
        ]:
    fl.basename = searchstr
    s = find_file(searchstr, exp_type=exp_type,
            expno=nodename,
            postproc=postproc, lookup=postproc_dict,
            dimname='indirect')
            #{{{filter data
    #s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    #}}}

    #fl.show();quit()
    s.ift('t2')
    fl.next('time domain cropped log')
    fl.image(s.C.setaxis('indirect','#').set_units('indirect','scan #').cropped_log())
    #{{{rough centering
    rough_center = abs(s).convolve('t2',10).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.next(' rough centered')
    fl.image(s.C.setaxis('indirect','#'))
    #}}}
    s.ft('t2')
    #fl.next('before clock correction')
    #fl.image(s.C.setaxis('indirect','#'))
    #s_clock = s['ph1',0]['ph2',1].sum('t2')
    #s.ift(['ph1','ph2'])
    #min_index = abs(s_clock).argmin('indirect',raw_index=True).item()
    #s_clock *= np.exp(-1j*clock_correction*s.fromaxis('indirect'))
    #s_clock['indirect',:min_index+1] *= -1
    #s_clock.sum('indirect').run(abs)
    #fl.next('clock correction')
    #fl.plot(s_clock,'.',alpha=0.7)
    #clock_correction = s_clock.argmax('clock_correction').item()
    #pyplot.axvline(x=clock_correction,alpha=0.5,color='r')
    #s *= np.exp(-1j*clock_correction*s.fromaxis('indirect'))
    #fl.next('after auto-clock correction')
    #s.ft(['ph1','ph2'])
    #fl.image(s.C.setaxis('indirect','#'))
    s.ift('t2')
    #{{{hermitian function test and apply best shift
    best_shift = hermitian_function_test(s[
        'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']],fl=fl)
    logger.info(strm("best shift is",best_shift))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    fl.next('time domain after hermitian test')
    fl.image(s.C.setaxis('indirect','#').set_units('indirect','scan #'))
    fl.next('frequency domain after hermitian test')
    s.ft('t2')
    fl.image(s.C.setaxis('indirect','#').set_units('indirect','scan #'))
    s.ift('t2')
    #}}}
    #{{{zeroth order phase correction
    ph0 = s['t2':0]['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
    logger.info(strm(ndshape(ph0)))
    #fl.show();quit()
    if len(ph0.dimlabels) > 0:
        assert len(ph0.dimlabels) == 1, repr(ndshape(ph0.dimlabels))+" has too many dimensions"
        ph0 = zeroth_order_ph(ph0, fl=fl)
        logger.info(strm('phasing dimension as one'))
    else:
        logger.info(strm("there is only one dimension left -- standard 1D zeroth order phasing"))
        ph0 = ph0/abs(ph0)
    s /= ph0
    fl.next('frequency domain -- after hermitian function test and phasing')
    s.ft('t2')
    fl.image(s.C.convolve('t2',10).C.setaxis(
'indirect','#').set_units('indirect','scan #'))
    #fl.show();quit()
    #}}}
    #{{{select t2 axis range and 
    s.convolve('t2',5) 
    s.ift('t2')
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    fl.next('where to cut')
    s = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']] 
    fl.plot(s,)
    #fl.show();quit()
    #{{{If visualizing via ILT
    fl.next('Plotting phased spectra')
    for j in range(ndshape(s)['indirect']):
        fl.plot(s['indirect',j]['t2':(-5e3,20e3)],
            alpha=0.5,
            label='vd=%g'%s.getaxis('indirect')[j])
    #quit()
    #exponential curve
    rec_curve = s['t2':(-5e3,20e3)].C.sum('t2')
    fl.next('recovery curve')
    fl.plot(rec_curve,'o')
    fl.show();quit()
    #{{{recovery curve and fitting
    # below, f_range needs to be defined
    M0,Mi,R1,vd = symbols("M_0 M_inf R_1 indirect",real=True)
    f,T1 = recovery(s, (-60,60),guess=None)
    fl.plot_curve(f,'inversion recovery curve')
    fl.show();quit()
    #attempting ILT plot with NNLS_Tikhonov_190104

    T1 = nddata(logspace(-3,3,150),'T1')
    l = sqrt(logspace(-8.0,0.5,35)) #play around with the first two numbers to get good l curve,number in middle is how high the points start(at 5 it starts in hundreds.)
    plot_Lcurve =False
    if plot_Lcurve:
        def vec_lcurve(l):
            return s.C.nnls('indirect',T1,lambda x,y: 1.0-2*exp(-x/y), l=l)

        x=vec_lcurve(l) 

        x_norm = x.get_prop('nnls_residual').data
        r_norm = x.C.run(linalg.norm,'T1').data

        with figlist_var() as fl:
            fl.next('L-Curve')
            figure(figsize=(15,10))
            fl.plot(log10(r_norm[:,0]),log10(x_norm[:,0]),'.')
            annotate_plot = True
            show_lambda = True
            if annotate_plot:
                if show_lambda:
                    for j,this_l in enumerate(l):
                        annotate('%0.3f'%this_l, (log10(r_norm[j,0]),log10(x_norm[j,0])),
                                ha='left',va='bottom',rotation=45)
                else:
                    for j,this_l in enumerate(l):
                        annotate('%d'%j, (log10(r_norm[j,0]),log10(x_norm[j,0])),
                                ha='left',va='bottom',rotation=45)
        d_2d = s*nddata(r_[1,1,1],r'\Omega')
    #fl.show();quit()
    offset = s.get_prop('proc')['OFFSET']
    this_l = 0.042#pick number in l curve right before it curves up
    o1 = 297.01 #o1 for free D2O
    sfo1 = s.get_prop('acq')['BF1']
    s.setaxis('t2',lambda x:
            x-o1)
    s.setaxis('t2',lambda x:
            x/(-sfo1)).set_units('t2','ppm')
    soln = s.real.C.nnls('indirect',T1, lambda x,y: 1.0-2.*exp(-x/y),l=this_l)
    soln.reorder('t2',first=False)
    soln.rename('T1','log(T1)')
    soln.setaxis('log(T1)',log10(T1.data))
    fl.next('w=6 in hexane')
    fl.image(soln.C.set_units('t2','ppm'))
    fl.show();quit()
    print("SAVING FILE")
    np.savez(searchstr+'_'+str(nodename)+'_ILT_inv',
            data=soln.data,
            logT1=soln.getaxis('log(T1)'),
            t2=soln.getaxis('t2'))               
    print("FILE SAVED")
    quit()


    #}}}

