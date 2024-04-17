# In[1]:

from pyspecdata import init_logging, find_file,figlist_var
from pyspecProcScripts.simple_functions import select_pathway
from pyspecProcScripts import fid_from_echo, lookup_table

init_logging(level="debug")

with figlist_var() as fl:
    s = find_file(
        "240227_E37_6-MSL_A1_Rasbatch240220_ODNP_1",
        exp_type="ODNP_NMR_comp/ODNP",
        expno="FIR_noPower",
        lookup=lookup_table,
    )
    s.mean("nScans")
    fl.next("raw")
    fl.image(s)
    fl.next("raw time domain")
    s.ift('t2')
    fl.image(s['t2':(None,80e-3)])
    s.ft('t2')
    # autoslice, phase and take FID slice
    s = fid_from_echo(s, signal_pathway = s.get_prop('coherence_pathway'), fl=fl) 
    s = select_pathway(s, s.get_prop('coherence_pathway'))
    fl.next("phased and FID sliced")
    fl.image(s)

# In[2]:


