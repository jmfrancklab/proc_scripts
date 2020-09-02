from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import fl_mod,find_echo_center, center_echo
from proc_scripts import postproc_dict
from sympy import symbols
from proc_scripts.fitting import decay
logger = init_logging("debug")
fl = fl_mod()
mpl.rcParams['figure.figsize'] = [8.0, 6.0]
rcParams["savefig.transparent"] = True
# {{{ input parameters
clock_correction = 0
filter_bandwidth = 5e3
t2 = symbols('t2')
# }}}
for searchstr, exp_type, nodename,flat_echo in [
        ('w8_200731','NMR_Data_AG',5,True)
        #('w8_1AT2RM_200731','NMR_Data_AG',4,True)
        #('w8_1AT4RM_200731','NMR_Data_AG',4,True)
        #('200303','T1CPMG_AER')
        ]:
    s = find_file(searchstr,exp_type=exp_type,
            expno=nodename,lookup=postproc_dict, fl=fl)
    fl.next('selected coherence')
    s = s['ph2',-1]['ph1',0]
    fl.image(s)
    #fl.show();quit()
    #this section is hard coded for flat echoes. I print the shape of s
    #to get the length of t2 and ensure it is an odd number. I then take 
    #the middle index and set this to 0. We will find a way to not have
    #this hard coded but for now this is what we have. 9/1/20
    if flat_echo:
        s['t2',16]=0
        center=find_echo_center(s,fl=fl)
        s = center_echo(s,center,fl=fl)
    else:    
        centers = []
        for j in range(ndshape(s)['indirect']):
            s_slice = s['indirect',j]
            this_center = find_echo_center(s_slice,fl=fl)
            centers.append(this_center)
        logger.info(centers)
        avg_center = sum(centers)/len(centers)
        s = center_echo(s, avg_center, fl=fl)
    #{{{Used to test if echo is flat or not
    #s = s['tE',20]['indirect',1]
    #fl.next('abs vs imag',legend=True)
    #fl.plot(abs(s),'-',label='abs')
    #fl.plot(s.imag,'--',label='imag')
    #fl.show();quit()
    #}}}
    fl.next('s centered')
    fl.image(s)
    fl.next('s centered in freq domain')
    s.ft('t2')
    fl.image(s)
    #fl.show();quit()

    #}}}
    #{{{slice out signal and sum along t2
    s = s['t2':(-154,154)]
    s.sum('t2')
    fl.next('summed along t2')
    fl.image(s)
    fl.show();quit()
    #}}}
    #{{{save to hdf5 file
    #s.name('w8_200731')
    #s.hdf5_write('w8_200731.h5')
    #{{{attempting ILT plot with NNLS_Tikhonov_190104
    tE_axis = s.getaxis('tE')
    Nx = 50
    Ny = 50
    Nx_ax = nddata(logspace(-5,3,Nx),'T1')
    Ny_ax = nddata(logspace(-5,3,Ny),'T2')
    data = s.C
    data.rename('indirect','tau1').setaxis('tau1',vd_list)
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
    np.savez('proc_'+searchstr+'_1',
            data = x.data,
            logT1 = x.getaxis('T1'),
            logT2 = x.getaxis('T2'))

    fl.show();quit()
