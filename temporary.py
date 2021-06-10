generate_image_files = False
if generate_image_files:
    import os; os.environ['pyspecdata_figures']='latex'
from pyspecdata import *
from proc_scripts import *
from proc_scripts import postproc_dict

init_logging(level="debug")
d = find_file(
    "210528_TEMPOL_3uM_cap_probe",
    exp_type="ODNP_NMR_comp",
    expno="enhancement",
    postproc="spincore_ODNP_v1",
    lookup=postproc_dict,
)
fl = figlist_var()
fl.next("freq")
fl.image(d.C.setaxis("power", "#").set_units("power", "scan #"))
d.ift("t2")
fl.next("time")
fl.image(d.C.setaxis("power", "#").set_units("power", "scan #"))
d.register_axis({'t2':0})
#d.ft('t2')
#d = d['t2':(-3e3,3e3)]
#d.ift('t2')
best_shift, window_size = hermitian_function_test(d["ph1", 1]["ph2", 0])
logger.info(strm("hermitian best shift",best_shift))
best_shift = 0.33e-3
logger.info(strm("manual best shift",best_shift))
d.setaxis("t2", lambda x: x - best_shift)
d.register_axis({"t2": 0})
d.ft('t2')
d = d['t2':(-3e3,3e3)]
fl.next("freq -- after phasing")
fl.image(d['t2':(-500,500)].C.setaxis("power", "#").set_units("power", "scan #"))
d.ift(['ph1','ph2'])
phasing = d['t2',0].C
phasing.data *= 0
phasing.ft(['ph1','ph2'])
phasing['ph1',1]['ph2',0] = 1
phasing.ift(['ph1','ph2'])
fl.next("ph,freq -- apply receiver phase")
d /= phasing
fl.image(d.C.setaxis("power", "#").set_units("power", "scan #"))
fl.next("ph,freq -- transients in phase")
ph0 = d.C.sum('t2')
ph0 /= abs(ph0)
d /= ph0
fl.image(d['t2':(-500,500)].C.setaxis("power", "#").set_units("power", "scan #"))
fl.next(r'after correlation alignment, $\varphi$ domain')
opt_shift,sigma = correl_align(d, indirect_dim='power',
        ph1_selection=0,ph2_selection=0,sigma=1e3)
d.ift('t2')
d *= np.exp(-1j*2*pi*opt_shift*d.fromaxis('t2'))
d.ft('t2')
fl.basename = None
fl.image(d['t2':(-500,500)].C.setaxis('power','#').set_units('power','scan #'))
fl.next(r'after masked correlation alignment, $\varphi$ domain')
opt_shift,sigma = correl_align(d, indirect_dim='power',
        ph1_selection=0,ph2_selection=0,sigma=30e3/2.35)
d.ift('t2')
d *= np.exp(-1j*2*pi*opt_shift*d.fromaxis('t2'))
d.ft('t2')
fl.image(d['t2':(-500,500)].C.setaxis('power','#').set_units('power','scan #'))
d.ft(['ph1','ph2'])
fl.next(r'after correlation alignment, $\Delta p$ domain')
fl.image(d['t2':(-500,500)].C.setaxis('power','#').set_units('power','scan #'))
fl.show("temp_correl_210601.png")
