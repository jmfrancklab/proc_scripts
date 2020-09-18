from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping
from proc_scripts import fl_mod, find_echo_center, center_echo, postproc_dict
from sympy import symbols
from proc_scripts.fitting import decay
logger = init_logging("debug")
fl = fl_mod()
mpl.rcParams['figure.figsize'] = [8.0, 6.0]
rcParams["savefig.transparent"] = True
# {{{ input parameters
filter_bandwidth = 5e3
t2 = symbols('t2')
test_for_flat_echo = False # test for flat echo and exit
write_h5 = True
read_h5 = True
# }}}
for searchstr, exp_type, nodename, flat_echo, clock_correction, h5_name, h5_dir in [
        #('w8_200731','NMR_Data_AG',5,True)
        ('w8_1AT2RM_200731','NMR_Data_AG',4,True,0,'T1CPMG_0920.h5','AG_processed_data')
        #('w8_1AT4RM_200731','NMR_Data_AG',4,True)
        #('200303','T1CPMG_AER','signal',False,1.785)
        ]:
    if write_h5:
        s = find_file(searchstr,exp_type=exp_type,
                expno=nodename,lookup=postproc_dict, fl=fl)
        fl.next('selected coherence')
        s = s['ph2',-1]['ph1',0]
        s *= exp(-1j*s.fromaxis('vd')*clock_correction)
        fl.image(s)
        # this section is hard coded for flat echoes. I print the shape of s
        # to get the length of t2 and ensure it is an odd number. I then take 
        # the middle index and set this to 0. We will find a way to not have
        # this hard coded but for now this is what we have. 9/1/20
        if flat_echo:
            # if the echo is flat, why were you using find_echo_center?
            # here I just manually tell it to use the middle point
            avg_center = s.getaxis('t2')[0].item() + diff(s.getaxis('t2')[r_[0,-1]]).item()
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
        fl.next('s centered in freq domain')
        s.ft('t2')
        fl.image(s)
        #}}}
        #{{{slice out signal and sum along t2
        s = s['t2':(-154,154)]
        s.sum('t2')
        fl.next('summed along t2')
        fl.image(s)
        #}}}
        #{{{save to hdf5 file
        s.name(searchstr) # use searchstr as the node name withing the HDF5 file
        s.hdf5_write(h5_name, directory=getDATADIR(h5_dir))
        #}}}
    if read_h5:
        s = nddata_hdf5(h5_name+'/'+searchstr, directory=getDATADIR(h5_dir))
        #{{{attempting ILT plot with NNLS_Tikhonov_190104
        tE_axis = s.getaxis('tE')
        vd_list = s.getaxis('indirect')
        Nx = 50
        Ny = 50
        x_name = r'$log(T_2/$s$)$'
        y_name = r'$log(T_1/$s$)$'
        Nx_ax = nddata(logspace(-5,3,Nx),x_name)
        Ny_ax = nddata(logspace(-5,3,Ny),y_name)
        s.rename('indirect','tau1').setaxis('tau1',vd_list)
        s.rename('tE','tau2').setaxis('tau2',tE_axis)
        s_ILT = s.C.nnls(('tau1','tau2'),
               (Nx_ax,Ny_ax),
               (lambda x1,x2: 1.-2*exp(-x1/x2),
                lambda y1,y2: exp(-y1/y2)),
                         l='BRD')

        s_ILT.set_units(x_name, None).set_units(y_name, None)
        fl.next('distributions')
        title(r'$T_{1} - T_{2} distribution$ for water loading 8 1AT/2RM')
        fl.image(s_ILT)
        s_ILT.name(searchstr+'_ILT') # use searchstr as the node name withing the HDF5 file
        s_ILT.hdf5_write(h5_name, directory=getDATADIR(h5_dir))
    fl.show()
