from pyspecdata import *
from scipy.optimize import leastsq,minimize
from proc_scripts import *
from proc_scripts.load_data import postproc_dict
from proc_scripts.fitting import recovery
from sympy import symbols
fl = fl_mod()
t2 = symbols('t2')
logger = init_logging("info")
# {{{ input parameters
filter_bandwidth = 5e3
coh_sel = {'ph1':0,
        'ph2':-1}
coh_err = {'ph1':1,# coherence channels to use for error
        'ph2':r_[0,2,3]}
# }}}

for searchstr,exp_type,nodename, postproc, clock_correction,this_l in [
        #('freeSL_201007','test_equip',5,'ag_IR2H',None)
        #('w8_200731', 'test_equip', 2, 'ag_IR2H',None),
        #('free4AT_201014','test_equip',3,'ag_IR2H',None)
        #('free4AT100mM_201104', 'test_equip',2,'ab_ir2h',None),
        #('ag_oct182019_w0_8','test_equip',3,'ab_ir2h',None)
        ('w20_201111','test_equip',2,'ab_ir2h',None,None),
        ]:
    fl.basename = searchstr
    if clock_correction is None:
        s = find_file(searchstr, exp_type=exp_type,
            expno=nodename,
            postproc=postproc, lookup=postproc_dict,
            dimname='indirect',fl=fl)
    else:
        s = find_file(searchstr, exp_type=exp_type,
            expno=nodename,
            postproc=postproc, lookup=postproc_dict,
            dimname='indirect',
            clock_correction=clock_correction)
    #fl.show();quit()
        #{{{filter data
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    #}}}
    #{{{hermitian function test and apply best shift
    fl.next('frequency domain before')
    fl.image(s)
    s.ift('t2')
    best_shift = hermitian_function_test(s[
        'ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']],fl=fl)
    logger.info(strm("best shift is",best_shift))
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    fl.next('time domain after hermitian test')
    fl.image(s)
    fl.next('frequency domain after')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    #}}}
    #{{{zeroth order phase correction
    ph0 = s['t2':0]['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]
    logger.info(strm(ndshape(ph0)))
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
    fl.image(s.C.convolve('t2',10))
    #fl.show();quit()
    #}}}
    #{{{select t2 axis range and 
    s.ift('t2')
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    s.ft('t2')
    fl.next('where to cut')
    s = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]*-1 
    fl.plot(s)
    #fl.show();quit()
    #{{{If visualizing via ILT
    fl.next('Plotting phased spectra')
    for j in range(ndshape(s)['indirect']):
        fl.plot(s['indirect',j]['t2':(-20,20)],
            alpha=0.5,
            label='vd=%g'%s.getaxis('indirect')[j])

    #exponential curve
    rec_curve = s['t2':(-20,20)].C.sum('t2')
    fl.next('recovery curve')
    fl.plot(rec_curve,'o')
    #fl.show();quit()
    #attempting ILT plot with NNLS_Tikhonov_190104

    T1 = nddata(logspace(-3,3,150),'T1')
    l = sqrt(logspace(-6.0,0.05,35)) #play around with the first two numbers to get good l curve,number in middle is how high the points start(at 5 it starts in hundreds.)
    plot_Lcurve = False
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
    sfo1 = 289.69
    arbitrary_reference = s.get_prop('acq')['BF1'] # will eventually be 
    print("SFO1 is",sfo1)
    s.setaxis('t2',lambda x:x + sfo1 - arbitrary_reference)
    this_l = 0.033 #pick number in l curve right before it curves up
    soln = s.real.C.nnls('indirect',T1, lambda x,y: 1.0-2.*exp(-x/y),l=this_l)
    soln.reorder('t2',first=False)
    soln.rename('T1','log(T1)')
    soln.setaxis('log(T1)',log10(T1.data))
    fl.next('water loading 20')
    fl.image(soln['t2':(100,300)])


    #fl.show();quit()
    print("SAVING FILE")
    np.savez(searchstr+'_'+str(nodename)+'_ILT_inv',
            data=soln.data,
            logT1=soln.getaxis('log(T1)'),
            t2=soln.getaxis('t2'))               
    print("FILE SAVED")
    quit()


    #}}}
    #{{{recovery curve and fitting
    s_sliced = s['ph2',coh_sel['ph2']]['ph1',coh_sel['ph1']]*-1 # bc inverted at high powers
    # below, f_range needs to be defined
    s_sliced.ift('t2')
    fl.next('recovery curve')
    fl.plot(s_sliced)
    fl.show();quit()
    M0,Mi,R1,vd = symbols("M_0 M_inf R_1 indirect",real=True)
    f,T1 = recovery(s_sliced, (-50,50),guess=None)
    fl.plot_curve(f,'inversion recovery curve')
        #{{{attempting ILT plot with NNLS_Tikhonov_190104
    T1 = nddata(logspace(-3,3,150),'T1')
    l = sqrt(logspace(-8.0,0.1,35)) #play around with the first two numbers to get good l curve,number in middle is how high the points start(at 5 it starts in hundreds.)
    if this_l is None:
        def vec_lcurve(l):
            return s.nnls('indirect',T1,lambda x,y: 1.0-2*exp(-x/y), l=l)

        x=vec_lcurve(l) 

        x_norm = x.get_prop('nnls_residual').data
        r_norm = x.run(linalg.norm,'T1').data

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
    #}}}
    #{{{setting axis to incorporate SFO1 based off of 
    # acqu file of data
    sfo1 = 251.76
    arbitrary_reference = s.get_prop('acq')['BF1'] 
    logger.info(strm("SFO1 is",sfo1))
    s.setaxis('t2',lambda x:x + sfo1 - arbitrary_reference)
    #}}}
    #{{{creating plot off of solution to L curve
    if this_l is not None:
        soln = s.real.nnls('indirect',T1, lambda x,y: 1.0-2.*exp(-x/y),l=this_l)
        soln.reorder('t2',first=False)
        soln.rename('T1','log(T1)')
        soln.setaxis('log(T1)',log10(T1.data))
        fl.next('solution')
        fl.image(soln['t2':(100,300)])
fl.show()

