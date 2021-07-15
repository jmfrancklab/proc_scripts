from pyspecdata import *
from pyspecProcScripts import *
from pyspecProcScripts import lookup_table
from pyspecProcScripts.third_level.process_enhancement import process_enhancement
from sympy import symbols, Symbol, latex,limit,init_printing
#{{{input parameters
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  # clear
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.9),  # 90% transparent white
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # clear
})
logger = init_logging("info")
t2 = symbols('t2')
signal_pathway = {'ph1':1}
# shouldn't be part of the pull request
freq_range= (-4e3,4e3)
t_range = (0,0.083)
fl = fl_mod()
#}}}
for filename,nodename,file_location,postproc in [
        ('210614_TEMPOL_500uM_cap_probe_DNP','enhancement',
            'ODNP_NMR_comp/test_equipment','spincore_ODNP_v1'),
        ]:
    s = find_file(filename,exp_type=file_location,expno=nodename,
            postproc=postproc,lookup=lookup_table,fl=fl)
    myslice = s['t2':freq_range]
    mysign = select_pathway(myslice, signal_pathway).real.sum('t2').run(np.sign)
    enhancement = process_enhancement(s,freq_range=freq_range, 
            t_range=t_range,sgn=mysign,avg_dim='nScans',fl=fl)
    fl.show()

