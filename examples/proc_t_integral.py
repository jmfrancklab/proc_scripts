import numpy as np
from numpy import r_
import pyspecdata as psd
import pyspecProcScripts as psdpr
import matplotlib.pyplot as plt
from cycler import cycler
from pylab import rcParams

thisdict = {}
default_colors = rcParams["axes.prop_cycle"].by_key()["color"]
thisls = [":", (0, (3, 1, 1, 1, 1, 1))]
prop_cycle = (cycler(ls=thisls)*3 + cycler(color = default_colors[0:6]))()
f_range = (-300,300)
excluded_pathway = [(0, 0), (0,3)]
with psd.figlist_var() as fl:
    data = psd.find_file(
            "241010_ssProbe_toroid_pmp90_16p7_echo",
            exp_type = "ODNP_NMR_comp/Echoes",
            expno="echo_1",
            lookup = psdpr.lookup_table)
    # {{{ Apply phase corrections
    data.ift("t2")
    data["t2"] -= data.getaxis("t2")[0]
    best_shift = psdpr.hermitian_function_test(
            psdpr.select_pathway(data.C.mean("nScans"),data.get_prop("coherence_pathway")))
    data.setaxis("t2", lambda x: x - best_shift).register_axis({"t2":0})
    data /= psdpr.zeroth_order_ph(psdpr.select_pathway(data["t2":0],data.get_prop("coherence_pathway")))
    # }}}
    # {{{ FID slice
    data = data["t2":(0,None)]
    data *= 2
    data["t2":0] *= 0.5
    # }}}
    fl.next("FID sliced data")
    data.ft("t2")
    fl.image(data)
    # {{{ Normalize so integral = 1
    thisdict['f_stop_orig'] = data.getaxis("t2")[-1]
    data.ift("t2")
    thisdict['t_stop_orig'] = data.getaxis("t2")[-1]
    thisdict["og_data"] = data.C
    # {{{ zero fill and make real/symmetric data
    N = len(data["nScans",0]["ph1",1].data)
    data.ft("t2", pad = 3 * N)
    data.run(np.real)
    data.ft_new_startpoint("t2","t")
    data.set_ft_prop("t2",None)
    data.ift("t2",shift=True)
    # }}}
    # {{{ apodize with a Lorenzian
    apo_fn =np.exp(-abs(data.fromaxis("t2")) / 0.01) 
    data *= apo_fn
    # }}}
    data.ft("t2")
    padded_data = data.C
    padded_data.ift("t2")
    nonzero_data = data.C
    nonzero_data.ift("t2")
    nonzero_data = nonzero_data["t2":(-thisdict["t_stop_orig"],thisdict["t_stop_orig"])]
    nonzero_data.ft("t2")
    nonzero_data.ift("t2")
    data.ift("t2")
    t_int = psdpr.t_integrate(og_data = thisdict["og_data"], 
            nonzero_data= nonzero_data, padded_data = padded_data, 
            apo_fn = apo_fn, frq_slice = f_range)
    # thisdict["og_data"].C.getaxis("t2")[-1] = nonzero_data.C.getaxis("t2"[-1]
    # apo_data = padded_data
    # {{{ make analytical sinc fn and multiply
