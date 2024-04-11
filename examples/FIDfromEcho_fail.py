# In[1]:

get_ipython().magic(u'load_ext pyspecdata.ipy')
from IPython.display import display, Markdown, Latex
import IPython.display as id
from pyspecdata import *
from pyspecProcScripts.simple_functions import select_pathway
from pyspecProcScripts import *
from pyspecProcScripts import fid_from_echo

# In[2]:

signal_pathway = {'ph1':1,'ph2':-2}
with figlist_var() as fl:
    for searchstr,exptype,nodename,postproc,freq_slice in [
        ['240227_E37_6-MSL_A1_Rasbatch240220_ODNP_1','ODNP_NMR_comp/ODNP','FIR_noPower',
            'spincore_nutation_v2',(-800,800)]
        ]:
        s = find_file(searchstr,exp_type=exptype,expno=nodename,postproc=None)
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
        s = fid_from_echo(s,signal_pathway,fl = fl)
        s = select_pathway(s,signal_pathway)
        fl.next('phased and FID sliced')
        fl.plot(s)

# In[3]:


