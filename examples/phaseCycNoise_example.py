from pylab import *
from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
from pyspecProcScripts.DCCT_func import DCCT
from pyspecProcScripts import *
save_fig = False
fl = figlist_var()
#this_figsize = (4.5,12)
this_figsize = (9,5.56)
for date,id_string,node_name in [
        ('210607','TEMPOL_100mM_cap_probe_DNP','FIR_27dBm'),
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    nodename = node_name
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'ODNP_NMR_comp/ODNP'))
    s = s['nScans',0]
    print(ndshape(s))
    s.reorder(['ph2','ph1','vd','t2'])
    s.setaxis('vd','#')
    s.ft(['ph1','ph2'])
    fl.next('raw data %s'%node_name)
    s['vd',0] *= -1
    s.setaxis('t2', lambda x: x - 1.06e-3)
    s.ft('t2',shift=True)
    s.ift(['ph1','ph2'])
    s *= exp(s.fromaxis('ph2')*2*pi*1j)
    fl.image(s)
    s.ift('t2')
    s /= exp(s.fromaxis('ph2')*2*pi*1j)
    s.ft(['ph1','ph2'])
    #temp = abs(s)**2
    #temp['ph2',1]['ph1',0] *= 0
    #temp.mean_all_but('t2')
    #s = (abs(s['ph2',1]['ph1',0])**2).mean_all_but('t2')
    #fl.next('time')
    #fl.plot(s/100)
    #fl.plot(temp)
    s.ift(['ph1','ph2'])
    #DCCT(s,fl.next('',figsize=this_figsize))
    signal_pathway = {"ph1":0, "ph2":1}
    s.ft('t2')
    opt_shift,sigma,mask_func = correl_align(s,
            indirect_dim='vd',
            signal_pathway=signal_pathway,
            sigma=3000 / 2.355,
            max_shift=300,
            fl=fl)
fl.show();quit()
