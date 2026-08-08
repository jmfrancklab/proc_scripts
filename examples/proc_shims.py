"""
This script takes in the datasets acquired with run_shimming in spincore_apps.
After properly phasing, the signal energy for each voltage and current is
calculated to print the final optimal voltage setting for the tested channel.
"""
from pylab import r_, exp, subplot
from pyspecdata import figlist_var, find_file

with figlist_var() as fl:
    for filename, nodename, color_tuple, label_str in [
        ("221213_ShimOn3", "shims_addr5_Z_yOpt", ("teal", "pink"), 
            ("opt y (run 4)", "opt y + opt z1 (run 4)")),
    ]:
        s = find_file(filename + ".h5", exp_type="ODNP_NMR_comp/echoes", expno=nodename)
        s = s.C.mean("nScans")
        tau = s.get_prop("acq_params")["tau_us"] - 10.0
        if "t2" not in s.dimlabels:
            s.set_units("t", "s")
            s.chunk("t", ["ph1", "ph2", "t2"], [4, 2, -1])
            s.setaxis("ph1", r_[0, 1, 2, 3] / 4)
            s.setaxis("ph2", r_[0, 2] / 4)
        s.reorder("t2", first=False)
        # {{{ DC offset correction
        max_t2 = s.getaxis("t2")[-1]
        rx_offset = s["t2" : (0.75 * max_t2, None)].mean("t2")
        s -= rx_offset
        s.ft("t2", shift=True)
        s.ift("t2")
        s.ft(["ph2", "ph1"])
        # }}}
        # center echo
        s.setaxis("t2", lambda x: x - tau * 1e-6).register_axis({"t2": 0})
        fl.next("y shims - NiSO4 %s" % nodename)
        ax = subplot(1, 2, 1)
        s.name("")
        voltage_setting = s.get_prop("acq_params")["voltage_setting"]
        for V_index, V_val in enumerate(voltage_setting):
            fl.plot(
                abs(s)["ph1", 1]["ph2", 0]["t2":(None, 0.08)]["I", V_index],
                alpha=0.5,
                label="voltage: %0.4f V" % V_val,
                ax=ax,
            )
        # {{{ gaussian apo to improve noise before
        #     plotting
        sigma = 0.09
        s = s["t2":(0, None)]
        s["t2", 0] *= 0.5
        g_filter = exp(-s.fromaxis("t2") ** 2 / 2 / sigma**2)
        fl.plot(g_filter["t2":(None, 0.08)] * abs(s).max(), ax=ax)
        s *= g_filter
        # }}}
        # {{{ zeroth order phase correction
        ph0 = s["t2", 0]
        ph0 /= abs(ph0)
        s /= ph0
        # }}}
        # {{{ calculate signal energy
        signal_E = (abs(s["t2":(None, 0.09)]["ph1", 1]["ph2", 0]) ** 2).sum("t2")
        signal_E_max = signal_E.argmax().item()
        signal_E_max_idx = list(s.getaxis("I")).index(signal_E_max)
        print("Signal energy maximum (current):", signal_E_max)
        print(
            "Signal energy maximum (voltage):",
            s.get_prop("acq_params")["voltage_setting"][signal_E_max_idx],
        )
        print("Index for signal energy maximum:", signal_E_max_idx)
        # }}}
        s.ft("t2")
        # {{{ plot data
        ax = subplot(1, 2, 2)
        for V_index, V_val in enumerate(voltage_setting):
            l = fl.plot(
                s["ph1", 1]["ph2", 0]["I", V_index],
                alpha=0.5,
                label="Voltage: %0.4f V" % V_val,
                ax=ax,
            )
            fl.plot(
                s["ph1", 1]["ph2", 0]["I", V_index].imag,
                alpha=0.1,
                label=None,
                ax=ax,
                color=l[-1].get_color(),
            )
        fl.next("NiSO4 - freq")
        fl.plot(
            s["t2":(-300, 300)]["ph1", 1]["ph2", 0]["I", 0],
            alpha=0.4,
            c=color_tuple[0],
            label=label_str[0],
        )
        fl.plot(
            s["t2":(-300, 300)]["ph1", 1]["ph2", 0]["I", signal_E_max_idx],
            alpha=0.4,
            c=color_tuple[1],
            label=label_str[1],
        )
        s.ift("t2")
        fl.next("NiSO4 - time domain")
        fl.plot(
            abs(s)["t2":(None, 0.1)]["ph1", 1]["ph2", 0]["I", 0],
            alpha=0.4,
            c=color_tuple[0],
            label=label_str[0],
        )
        fl.plot(
            abs(s)["t2":(None, 0.1)]["ph1", 1]["ph2", 0]["I", signal_E_max_idx],
            alpha=0.4,
            c=color_tuple[1],
            label=label_str[1],
        )
        # }}}
