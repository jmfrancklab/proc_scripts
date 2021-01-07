from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import fl_mod, find_echo_center, center_echo, postproc_dict
from sympy import symbols
import sympy as sp
from proc_scripts.fitting import decay
logger = init_logging("debug")
fl = fl_mod()
mpl.rcParams['figure.figsize'] = [8.0, 6.0]
rcParams["savefig.transparent"] = True
# {{{ input parameters
filter_bandwidth = 5e3
t2 = symbols('t2')
test_for_flat_echo = False# test for flat echo and exit
write_h5 = True
read_h5 = False 
# }}}
for searchstr, exp_type, nodename, flat_echo, clock_correction, h5_name, h5_dir in [
        #('w8_2RM1AT_201008','test_equip',4,False,0,'T1CPMG_201008_w8_2RM1AT.h5','process_data_AG')
        #('w8_201008','test_equip',3,False,0,'T1CPMG_201008_w8.h5','process_data_AG')
        #('free4AT_201008','test_equip',6,False,0,'T1CPMG_201008_FreeAT.h5','process_data_AG'),
        ('free4AT100mM_201104','test_equip',3,True,0,'T1CPMG_201104_Free4AT100mM.h5','process_data_AG'),
        #('free4AT_201014','test_equip',7,True,0,'T1CPMG_201014_FreeAT_1.h5','process_data_AG')
        #('w8_200731','test_equip',5,False,0,'T1CPMG_200731.h5','process_data_AG')
        #('w8_1AT2RM_200731','test_Equip',4,True,0,'T1CPMG_0920.h5','AG_processed_data')
        #('w8_1AT4RM_200731','NMR_Data_AG',4,True)
        #('200303','T1CPMG_AER','signal',False,1.785)
        ]:
    if write_h5:
        #before running go into the preprocessing in load_data as some parameters are hardcoded. Double check these
        s = find_file(searchstr,exp_type=exp_type,
                expno=nodename,lookup=postproc_dict, fl=fl)
        #fl.show();quit()
        fl.next('selected coherence')
        s = s['ph2',-1]['ph1',0]
        s *= exp(-1j*s.fromaxis('indirect')*clock_correction)
        fl.image(s)
        # this section is hard coded for flat echoes. I print the shape of s
        # to get the length of t2 and ensure it is an odd number. I then take 
        # the middle index and set this to 0. We will find a way to not have
        # this hard coded but for now this is what we have. 9/1/20
        s.ift('t2')
        if flat_echo:
           s['t2',16]=0
           avg_center=find_echo_center(s,fl=fl)
        else:    
            centers = []
            for j in range(ndshape(s)['indirect']):
                s_slice = s['indirect',j]
                this_center = find_echo_center(s_slice,fl=fl)
                centers.append(this_center)
            logger.info(centers)
            avg_center = sum(centers)/len(centers)
        s = center_echo(s, avg_center, fl=fl)
        if test_for_flat_echo:
            #{{{Used to test if echo is flat or not
            s = s['tE',20]['indirect',1]
            fl.next('abs vs imag',legend=True)
            fl.plot(abs(s),'-',label='abs')
            fl.plot(s.imag,'--',label='imag')
            fl.show();quit()
            #}}}
        fl.next('s centered')
        fl.image(s)
        fl.next('s centered in frequency domain')
        s.ft('t2')
        fl.image(s)
        #fl.show();quit()
        #}}}
        #{{{slice out signal and sum along t2
        s = s['t2':(-300,300)]
        s.sum('t2')
        fl.next('summed along t2')
        fl.image(s)
        s *= -1
        #}}}
        #{{{CPMG decay curve
        CPMG = s['indirect',-1]
        fl.next('CPMG curve')
        fl.plot(CPMG,'o')
        fit_CPMG = fitdata(CPMG)
        M0,R2,vd = symbols("M_0 R_2 tE",real=True)
        fit_CPMG.functional_form = (M0)*sp.exp(-vd*R2)
        logger.info(strm("Functional Form", fit_CPMG.functional_form))
        logger.info(strm("Functional Form", fit_CPMG.functional_form))
        fit_CPMG.fit()
        print("output:",fit_CPMG.output())
        print("latex:",fit_CPMG.latex())
        T2 = 1./fit_CPMG.output('R_2')
        fl.next('CPMG fit of T1CPMG for free D2O')
        fl.plot(CPMG,'o',label='data')
        fl.plot(fit_CPMG.eval(100),label='fit')
        print("T2 IS:",T2)
        fl.show();quit()
        #{{{IR recovery curve
        IR = s['tE',0]
        fl.next('IR')
        fl.plot(IR,'o')
        f = fitdata(IR)
        M0,Mi,R1,vd = symbols("M_0 M_inf R_1 indirect",real=True)
        f.functional_form = Mi + (M0-Mi)*sp.exp(-vd*R1)
        logger.info(strm("Functional Form", f.functional_form))
        logger.info(strm("Functional Form", f.functional_form))
        f.fit()
        print("output:",f.output())
        print("latex:",f.latex())
        T1 = 1./f.output('R_1')
        fl.next('fit for IR portion of T1CPMG')
        fl.plot(IR,'o',label='fake data')
        fl.plot(f.eval(100),label='fit',human_units=False)
        T1 = 1./f.output('R_1')
        #}}}
        fl.show();quit()
        #}}}
        #{{{save to hdf5 file
        s.name(searchstr) # use searchstr as the node name withing the HDF5 file
        s.hdf5_write(h5_name, directory=getDATADIR(h5_dir))
        print(ndshape(s))
        #}}}
    if read_h5:
        s = nddata_hdf5(h5_name+'/'+searchstr,getDATADIR(exp_type=h5_dir)) 
        #{{{attempting ILT plot with NNLS_Tikhonov_190104
        tE_axis = s.getaxis('tE')
        vd_list = s.getaxis('indirect')
        Nx = 50
        Ny = 50
        Nx_ax = nddata(logspace(-5,3,Nx),'T1')
        Ny_ax = nddata(logspace(-5,3,Ny),'T2')
        s_ILT=s.C
        s_ILT.rename('indirect','tau1').setaxis('tau1',vd_list)
        s_ILT.rename('tE','tau2').setaxis('tau2',tE_axis)
        s_ILT = s_ILT.C.nnls(('tau1','tau2'),
               (Nx_ax,Ny_ax),
               (lambda x1,x2: 1.-2*exp(-x1/x2),
                lambda y1,y2: exp(-y1/y2)),
                         l='BRD')
        s_ILT.setaxis('T1',log10(Nx_ax.data)).set_units('T1',None)
        s_ILT.setaxis('T2',log10(Ny_ax.data)).set_units('T2',None)
        #figure()
        #title(r'$T_{1} - T_{2} distribution$ for free AT in soln')
        fl.next('free 20 mM 4-AT in solution')
        fl.image(s_ILT)
        xlabel(r'$log(T_2/$s$)$')
        ylabel(r'$log(T_1/$s$)$')
        fl.show();quit()
        s_ILT.name(searchstr+'_ILT') # use searchstr as the node name withing the HDF5 file
        #s_ILT.hdf5_write(h5_name, directory=getDATADIR(h5_dir))
    fl.show()
