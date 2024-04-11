from pyspecdata import *
from pyspecProcScripts.simple_functions import select_pathway
from pyspecProcScripts import *
from pyspecProcScripts import fid_from_echo, lookup_table
init_logging(level='debug')

signal_pathway = {'ph1':1,'ph2':-2}
with figlist_var() as fl:
    # in the following line, there is a problem because the postproc should
    # always be automatically selected! 
    s = find_file('240227_E37_6-MSL_A1_Rasbatch240220_ODNP_1',exp_type='ODNP_NMR_comp/ODNP',expno='FIR_noPower',lookup=lookup_table)
    s.squeeze()
    s.ft('t2',shift = True)
    s.ft('ph1',unitary = True)
    s *= s.shape['nScans']
    s.mean('nScans')
    s.reorder(['ph1','ph2','vd','t2'])
    fl.next('raw')
    fl.image(s)
    print(s.shape)
    s.ift('t2')
    s.set_units('t2','s')
    s.ft('t2')
    fl.basename = nodename
    s = s['vd',0]
    # autoslice, phase and take FID slice
    s = fid_from_echo(s,signal_pathway,fl=fl)#,peak_lower_thresh = 0.01)
    s = select_pathway(s,signal_pathway)
    fl.next('phased and FID sliced')
    fl.plot(s)
