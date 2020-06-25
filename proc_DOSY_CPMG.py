from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import *
from sympy import symbols
import os
init_logging('debug')

t2 = symbols('t2')
expno = 3
filename = 'ab_jun172020_w0_5'
with figlist_var(file_name=filename+'.pdf') as fl:
    d = find_file(filename,
                  expno=expno,
                 exp_type='NMR_Data_AAB')
    # {{{ all of this would be your "preprocessing" and would be tied to the name of your pulse sequence
    l22 = int(d.get_prop('acq')['L'][22]) # b/c the l are integers by definition
    l25 = int(d.get_prop('acq')['L'][25])
    d12 = d.get_prop('acq')['D'][12]
    d11 = d.get_prop('acq')['D'][11]
    p1 = d.get_prop('acq')['P'][1]
    ppg = d.get_prop('pulprog')
    # {{{ these are explanatory -- maybe comment them out?
    m = re.search('(.*dwdel1=.*)',ppg,flags=re.IGNORECASE)
    print(m.groups()) # show the line that sets dwdel1
    # then look for de and depa
    print([(j,d.get_prop('acq')[j]) for j in d.get_prop('acq').keys() if 'de' in j.lower()])
    # I actually can't find depa
    # }}}
    m = re.search('\ndefine list<grad_scalar> gl1 = {(.*)}',ppg)
    grad_list = array([float(j.group()) for j in re.finditer('([0-9.]+)',m.groups()[0])])
    m = re.search('([0-9.]+) G/mm', d.get_prop('gradient_calib'))
    grad_list *= float(m.groups()[0])*0.1
    dwdel1 = 3.5e-6 # where does this come from? DE is actually larger than this?
    # {{{ find anavpt without hard-setting
    m = re.search('"anavpt=([0-9]+)"',ppg)
    if m is None:
        raise ValueError("I can't find anavpt in the pulse sequence")
    anavpt = int(m.groups()[0])
    # }}}
    dwdel2 = (anavpt*0.05e-6)/2
    TD = d.get_prop('acq')['TD2']
    quadrature_points = TD/2
    num_points_per_echo = quadrature_points/l25
    acq_time = dwdel2*num_points_per_echo*2
    # {{{ so, in principle, later, we can/should do what I did above (w/ eval),
    # but it's getting crazy now, so I stop for now
    tau_extra = 20e-6
    tau_pad = tau_extra-6e-6
    tau_pad_start = tau_extra-dwdel1-6e-6
    tau_pad_end = tau_extra-6e-6
    tE = dwdel1 + 5e-6 + tau_pad_start + 1e-6 + num_points_per_echo*(dwdel2*2) + tau_pad_end
    # }}}

    d.chunk('indirect',['indirect','phcyc'],[l22,-1])
    d.chunk('phcyc',['ph8','ph4','m','n'],[2,2,2,2])
    d.setaxis('ph8',r_[0.,2.]/4)
    d.setaxis('ph4',r_[0.,2.]/4)
    d.setaxis('m',r_[0,2.]/4)
    d.setaxis('n',r_[0,2.]/4)
    d.ft(['ph8','ph4','m','n'])
    d.reorder(['m','n','ph4','ph8','indirect','t2'])
    d.setaxis('indirect',grad_list)
    fl.next('abs raw data')
    fl.image(abs(d))
    d.chunk('t2',['echo','t2'],[l25,-1])
    d.reorder(['m','n','ph4','ph8','indirect','echo','t2'])
    d.ft('t2', shift=True).ift('t2') # this is overkill -- need a pyspecdata function that does this w/out the fft
    # }}}
    s = d['ph8',0]['ph4',1]['m',1]['n',0]

    # coarse phasing before hermitian_function_test doesn't seem to be required
    echo_center = hermitian_function_test(s['indirect',0],fl=fl)
    print("echo center is",echo_center)
    s.setaxis('t2', lambda x: x-echo_center)
    fl.next('request 5 (after subtracting best shift)')
    s.register_axis({'t2':0})
    #s /= zeroth_order_ph(s['t2',0], fl=fl)
    ph0 = s['t2',0].C.mean('t2')
    ph0 /= abs(ph0)
    s /= ph0
    fl.image(s)
    time_bound = min(abs(s.getaxis('t2')[r_[0,-1]]))
    s = s['t2':(-time_bound,time_bound)]
    assert isclose(s.getaxis('t2')[0],-s.getaxis('t2')[-1]), strm("if this doesn't work, talk to John b/c it means the pyspecdata slice notation isn't inclusive and I think we want it to be (?)",s.getaxis('t2')[r_[0,-1]])
    s.ft('t2')
    fl.next('abs: request 3')
    fl.image(abs(s))
    fl.next('request 3')
    fl.image(s)
    perform_fitting = False
    if perform_fitting:
        f_axis = s.fromaxis('t2')
        def costfun(p):
            zeroorder_rad,firstorder = p
            phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
            phshift *= exp(-1j*2*pi*zeroorder_rad)
            corr_test = phshift * s
            return (abs(corr_test.data.imag)**2)[:].sum()
        iteration = 0
        def print_fun(x, f, accepted):
            global iteration
            iteration += 1
            logger.info(strm(iteration, x, f, int(accepted)))
            return
        sol = basinhopping(costfun, r_[0.,0.],
                minimizer_kwargs={"method":'L-BFGS-B'},
                callback=print_fun,
                stepsize=100.,
                niter=100,
                T=1000.
                )
        zeroorder_rad, firstorder = sol.x
        phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
        phshift *= exp(-1j*2*pi*zeroorder_rad)
        s *= phshift
        print("RELATIVE PHASE SHIFT WAS %0.1f\\us and %0.1f$^\circ$", firstorder, angle(zeroorder_rad)/pi*180)
        fl.next('after phased - real ft')
        fl.image(s.real)
        fl.next('after phased - imag ft')
        fl.image(s.imag)
        if s['echo',0].data[:].sum().real < 0:
            s *= -1
        fl.next('plot')
        fl.plot(s)
        print("Ok")
        save_data = False
        if save_data:
            np.savez('proc_DOSY_CPMG_'+filename+'_expno'+str(expno),
                    data = final_spec.data,
                    indirect = final_spec.getaxis('indirect'),
                    echo = final_spec.getaxis('echo'),
                    t2 = final_spec.getaxis('t2'),
                    )
