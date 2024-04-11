from pyspecdata import *
from pyspecProcScripts.simple_functions import select_pathway
from pyspecProcScripts import *
from pyspecProcScripts import fid_from_echo, lookup_table
init_logging(level='debug')

signal_pathway = {'ph1':1,'ph2':-2}
with figlist_var() as fl:
    # in the following line, there is a problem because the postproc should
    # always be automatically selected! 
    # 
    # in this case, it uses the function proc_spincore_IR
    s = find_file('240227_E37_6-MSL_A1_Rasbatch240220_ODNP_1',exp_type='ODNP_NMR_comp/ODNP',expno='FIR_noPower',lookup=lookup_table)
    # the FID slice makes the assumption that tau (from the ppg) is subtracted this should be done in the postproc function!
    s.squeeze()
    s *= s.shape['nScans'] # should be done in postproc function!
    s.mean('nScans')
    s.reorder(['ph1','ph2','vd','t2']) # should be done in postproc function! (I think it is?)
    fl.next('raw')
    fl.image(s)
    print(s.shape) # why?
    s.ift('t2')
    s.set_units('t2','s') # this is done in the postproc
    s.ft('t2')
    fl.basename = nodename
    s = s['vd',0]
    # autoslice, phase and take FID slice
    s = fid_from_echo(s,signal_pathway,fl=fl)#,peak_lower_thresh = 0.01)
    s = select_pathway(s,signal_pathway)
    fl.next('phased and FID sliced')
    fl.plot(s)
