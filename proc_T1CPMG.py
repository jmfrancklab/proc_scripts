from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import *
from proc_scripts import postproc_dict
from sympy import symbols
fl = fl_mod()
mpl.rcParams['figure.figsize'] = [8.0, 6.0]
rcParams["savefig.transparent"] = True
# {{{ input parameters
clock_correction = 0
filter_bandwidth = 5e3
t2 = symbols('t2')
# }}}
for searchstr, exp_type, nodename, postproc in [
        ('w8_200731','test_equip',5,'ag_T1CPMG_2h')
        #('200303','T1CPMG_AER')
        ]:
    s = find_file(searchstr,exp_type=exp_type,
            expno=nodename, postproc=postproc,
            lookup=postproc_dict)
    fl.show();quit()
    #{{{ pulling acq params
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    nScans = s.get_prop('acq_params')['nScans']
    p90_s = s.get_prop('acq_params')['p90_us']*1e-6
    deadtime_s = s.get_prop('acq_params')['deadtime_us']*1e-6
    deblank_s = s.get_prop('acq_params')['deblank_us']*1e-6
    marker_s = s.get_prop('acq_params')['marker_us']*1e-6
    tau1_s = s.get_prop('acq_params')['tau1_us']*1e-6
    pad_start_s = s.get_prop('acq_params')['pad_start_us']*1e-6
    pad_end_s = s.get_prop('acq_params')['pad_end_us']*1e-6
    #}}}
    s.set_units('t','s')
    fl.next('raw data - no clock correction')
    fl.image(s)
    fl.next('raw data - clock correction')
    s *= exp(-1j*s.fromaxis('vd')*clock_correction)
    fl.image(s)
    #{{{ set up of CPMG axis -- chunking
    orig_t = s.getaxis('t')
    acq_time_s = orig_t[nPoints]
    twice_tau = deblank_s + 2*p90_s + deadtime_s + pad_start_s + acq_time_s + pad_end_s + marker_s
    t2_axis = linspace(0,acq_time_s,nPoints)
    tE_axis = r_[1:nEchoes+1]*twice_tau
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('tE',tE_axis)
    s.setaxis('t2',t2_axis).set_units('t2','s')
    vd_list = s.getaxis('vd')
    s.setaxis('vd',vd_list).set_units('vd','s')
    #}}}
    fl.next('chunked data - coherence channels')
    s.ft(['ph1'])
    fl.image(s)
    fl.next('filtered + rough centered data')
    # centering here means about 0 - this is mostly just an axis adjustment
    s.ft('t2',shift=True)
    s = s['t2':(-filter_bandwidth/2,filter_bandwidth/2)]
    s.ift('t2')
    rough_center = abs(s).convolve('t2',0.01).mean_all_but('t2').argmax('t2').item()
    s.setaxis(t2-rough_center)
    fl.image(s)
    s.ft('t2')
    # here we would align the frequencies which is hereto unresolved
    s.ift('t2')
    residual,best_shift = hermitian_function_test(s[
        'ph1',1])
    fl.next('hermitian test')
    fl.plot(residual)
    print("best shift is",best_shift)
    # {{{ slice out the FID appropriately and phase correct
    # it
    s.setaxis('t2', lambda x: x-best_shift).register_axis({'t2':0})
    fl.next('time domain after hermitian test')
    fl.image(s)
    ph0 = s['t2':0]['tE',0]['ph1',1]
    print(ndshape(ph0))
    ph0 = zeroth_order_ph(ph0, fl=fl)
    print('phasing dimension as one')
    s /= ph0
    fl.next('frequency domain -- after hermitian function test and phasing')
    s.ft('t2')
    fl.image(s)
    s.ift('t2')
    fl.next('check phasing -- real')
    fl.plot(s['ph1',1]['tE',0])
    gridandtick(gca())
    fl.next('check phasing -- imag')
    fl.plot(s['ph1',1]['tE',0].imag)
    gridandtick(gca())
    # not sure that I really like how real vs imag looks here
    # looks like signal at t2 = 0 gets phased properly, but anything else not
    s = s['t2':(0,None)]
    s['t2',0] *= 0.5
    fl.next('phased and FID sliced')
    fl.image(s)
    fl.next('phased and FID sliced -- frequency domain')
    s.ft('t2')
    fl.image(s)
    fl.next('signal vs. vd')
    s_sliced = s['ph1',1]
    s_forerror = s['ph1',0]
    # variance along t2 gives error for the mean, then average it across all
    # other dimension, then sqrt for stdev
    # could not get error to work, given the echo dimension
    fl.image(s_sliced.sum('t2'))
    d = s_sliced.real
    print("CONSTRUCTING KERNELS...")
    Nx = 100
    Ny = 100
    Nx_ax = nddata(logspace(-5,3,Nx),'T1')
    Ny_ax = nddata(logspace(-5,3,Ny),'T2')
    data = d.C
    data.rename('vd','tau1').setaxis('tau1',vd_list)
    data.rename('tE','tau2').setaxis('tau2',tE_axis)
    x = data.C.nnls(('tau1','tau2'),
           (Nx_ax,Ny_ax),
           (lambda x1,x2: 1.-2*exp(-x1/x2),
            lambda y1,y2: exp(-y1/y2)),
                     l='BRD')

    x.setaxis('T1',log10(Nx_ax.data)).set_units('T1',None)
    x.setaxis('T2',log10(Ny_ax.data)).set_units('T2',None)
    figure()
    title(r'$T_{1} - T_{2} distribution$')
    image(x)
    xlabel(r'$log(T_2/$s$)$')
    ylabel(r'$log(T_1/$s$)$')
    np.savez('proc_'+date+'_'+id_string+'_1',
            data = x.data,
            logT1 = x.getaxis('T1'),
            logT2 = x.getaxis('T2'))
fl.show();quit()
