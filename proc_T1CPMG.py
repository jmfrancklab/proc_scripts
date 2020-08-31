from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import fl_mod,find_echo_center, center_echo
from proc_scripts import postproc_dict
from sympy import symbols
from proc_scripts.fitting import decay
logger = init_logging("debug")
logger.info("this is a test")
logger.info("this is a test debug")
fl = fl_mod()
mpl.rcParams['figure.figsize'] = [8.0, 6.0]
rcParams["savefig.transparent"] = True
# {{{ input parameters
clock_correction = 0
filter_bandwidth = 5e3
t2 = symbols('t2')
# }}}
for searchstr, exp_type, nodename,this_l in [
        ('w8_200731','NMR_Data_AG',5,0.088)
        #('200303','T1CPMG_AER')
        ]:
    s = find_file(searchstr,exp_type=exp_type,
            expno=nodename,lookup=postproc_dict, fl=fl)
    centers = []
    for j in range(ndshape(s)['nScans']):
        s_slice = s['nScans',j]
        this_center = find_echo_center(s_slice,fl=fl)
        centers.append(this_center)
    logger.info(centers)
    avg_center = sum(centers)/len(centers)
    s = center_echo(s, avg_center,fl=fl)
    print(ndshape(s))
    fl.next('abs vs imag',legend=True)
    fl.plot(abs(s['tE',20]),'-',label='abs')
    fl.plot(s.imag['tE',20],'--',label='imag')
    fl.show();quit()
    fl.next('s centered')
    fl.image(s)
    fl.next('s centered in freq domain')
    s.ft('t2')
    fl.image(s)
    fl.show();quit()

    #}}}
    #{{{slice out signal and sum along t2
    s = s['t2':(-30,30)]
    s.sum('t2')
    fl.next('summed along t2')
    fl.image(s)
    fl.next('sliced data in time domain')
    s.ift('t2')
    fl.image(s)

    fl.show();quit()
    #}}}
    #{{{save to hdf5 file
    #s.name('w8_200731')
    #s.hdf5_write('w8_200731.h5')
    #{{{attempting ILT plot with NNLS_Tikhonov_190104
    vd_list = s.getaxis('indirect')
    tE_axis = s.getaxis('echoes')
    Nx = 10
    Ny = 10
    Nx_ax = nddata(logspace(-5,3,Nx),'T1')
    Ny_ax = nddata(logspace(-5,3,Ny),'T2')
    data = s.C
    data.rename('indirect','tau1').setaxis('tau1',vd_list)
    data.rename('echoes','tau2').setaxis('tau2',tE_axis)
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
    np.savez('proc_'+searchstr+'_1',
            data = x.data,
            logT1 = x.getaxis('T1'),
            logT2 = x.getaxis('T2'))

    fl.show();quit()
