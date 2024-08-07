from pyspecdata import *
from scipy.optimize import leastsq, minimize, basinhopping
from pyspecProcScripts import (
    fl_mod,
    center_echo,
    hermitian_function_test,
    lookup_table,
    DCCT,
)
from sympy import symbols
import sympy as sp
from pyspecProcScripts.fitting import decay
import matplotlib.pyplot as plt

logger = init_logging("debug")
fl = figlist_var()
plt.rcParams["figure.figsize"] = [8.0, 6.0]
rcParams["savefig.transparent"] = True
# {{{ input parameters
filter_bandwidth = 5e3
t2 = symbols("t2")
test_for_flat_echo = False  # test for flat echo and exit
# Below we give the option of JUST writing the processed
# data to a file (write_h5) and then the option of reading
# and plotting/imaging the processed result (read_h5) as
# time saver.
clock_correction = False
write_h5 = True  # writes the hdf5 file
read_h5 = False  # reads the completed hdf5 file
# }}}
for (
    searchstr,
    exp_type,
    nodename,
    flat_echo,
    clock_correction,
    freq_slice,
    t_slice,
    h5_name,
    h5_dir,
) in [
    # ('w8_2RM1AT_201008','test_equip',4,False,0,'T1CPMG_201008_w8_2RM1AT.h5','process_data_AG')
    # ('w8_201008','test_equip',3,False,False,(None,None),'T1CPMG_201008_w8.h5','process_data_AG')
    # ('free4AT_201008','test_equip',6,False,0,'T1CPMG_201008_FreeAT.h5','process_data_AG'),
    # ('free4AT100mM_201104','test_equip',3,True,None,(-2e3,2e3),
    #    'T1CPMG_201104_Free4AT100mM.h5','process_data_AG'),
    # ('free4AT_201014','test_equip',7,True,None,(None,None),'T1CPMG_201014_FreeAT_1.h5','process_data_AG')
    # ('w8_200731','test_equip',5,False,0,(None,None),'T1CPMG_200731.h5','process_data_AG')
    # ('w8_1AT2RM_200731','test_Equip',4,True,0,'T1CPMG_0920.h5','AG_processed_data')
    # ('w8_1AT4RM_200731','NMR_Data_AG',4,True)
    (
        "w3_201111",
        "test_equip",
        3,
        False,
        False,
        (-76, 76),
        (None, 0.1),
        "w3_T1CPMG",
        "process_data_AG",
    )
]:
    # If file is not processed yet, write_h5 should be True above in the input
    # parameters
    # If file is already written one can declare write_h5 as False to save time
    if write_h5:
        # before running go into the preprocessing in load_data as some parameters are hardcoded. Double check these
        s = find_file(
            searchstr,
            exp_type=exp_type,
            expno=nodename,
            lookup=lookup_table,
            fl=fl,
        )
        s.ift("t2")
        # {{{ clock correction
        if clock_correction:
            clock_corr = nddata(np.linspace(-3, 3, 2500), "clock_corr")
            s.ft("t2")
            fl.next("before clock correction")
            fl.image(as_scan_nbr(s))
            s_clock = s["ph1", 0]["ph2", -1].sum("t2")
            s.ift(["ph1", "ph2"])
            min_index = abs(s_clock).argmin("vd", raw_index=True).item()
            s_clock *= np.exp(-1j * clock_corr * s.fromaxis("vd"))
            s_clock["vd", : min_index + 1] *= -1
            s_clock.sum("vd").run(abs)
            fl.next("clock correction")
            fl.plot(s_clock, ".", alpha=0.7)
            clock_corr = s_clock.argmax("clock_corr").item()
            pyplot.axvline(x=clock_corr, alpha=0.5, color="r")
            s *= np.exp(-1j * clock_corr * s.fromaxis("vd"))
            s.ft(["ph1", "ph2"])
            fl.next("after auto-clock correction")
            fl.image(s.C.setaxis("vd", "#"))
            s.ift("t2")
            # }}}
        if test_for_flat_echo:
            # {{{Used to test if echo is flat or not
            s = s["tE", 20]["indirect", 1]
            fl.next("abs vs imag", legend=True)
            fl.plot(abs(s), "-", label="abs")
            fl.plot(s.imag, "--", label="imag")
            fl.show()
            quit()  # this flag is used to ask for a partial test, and so this is the *only* valid way to use quit
        if flat_echo:
            s["t2", 16] = 0
            avg_center = hermitian_function_test(s)
        else:
            centers = []
            for j in range(ndshape(s)["indirect"]):
                s_slice = s["indirect", j]
                this_center = hermitian_function_test(s_slice)
                centers.append(this_center)
            logger.info(centers)
            avg_center = sum(centers) / len(centers)
        s = center_echo(s, avg_center, fl=fl)
        # }}}
        fl.next("s centered")
        fl.image(s.C.setaxis("indirect", "#").set_units("indirect", "scan #"))
        # {{{slice out signal and sum along t2
        s = s["ph1", 0]["ph2", -1]
        s.sum("t2")
        fl.next("summed along t2")
        fl.image(s.C.setaxis("indirect", "#").set_units("indirect", "scan #"))
        s *= -1

        # }}}
        # {{{CPMG decay curve
        CPMG = s["indirect", -1]
        fl.next("CPMG curve")
        fl.plot(CPMG, "o")
        fit_CPMG = fitdata(CPMG)
        M0, R2, vd = symbols("M_0 R_2 tE", real=True)
        fit_CPMG.functional_form = (M0) * sp.exp(-vd * R2)
        logger.info(strm("Functional Form", fit_CPMG.functional_form))
        logger.info(strm("Functional Form", fit_CPMG.functional_form))
        fit_CPMG.fit()
        logger.info(strm(" CPMG output:", fit_CPMG.output()))
        logger.info(strm("CPMG latex:", fit_CPMG.latex()))
        T2 = 1.0 / fit_CPMG.output("R_2")
        fl.next("CPMG fit of T1CPMG for free D2O")
        fl.plot(CPMG, "o", label="data")
        fl.plot(fit_CPMG.eval(100), label="fit")
        logger.info(strm("T2 IS:", T2))
        # {{{IR recovery curve
        IR = s["tE", 0]
        fl.next("IR")
        fl.plot(IR, "o")
        f = fitdata(IR)
        M0, Mi, R1, vd = symbols("M_0 M_inf R_1 indirect", real=True)
        f.functional_form = Mi + (M0 - Mi) * sp.exp(-vd * R1)
        logger.info(strm("Functional Form", f.functional_form))
        logger.info(strm("Functional Form", f.functional_form))
        f.fit()
        logger.info(strm("IR output:", f.output()))
        logger.info(strm("IR latex:", f.latex()))
        T1 = 1.0 / f.output("R_1")
        fl.next("fit for IR portion of T1CPMG")
        fl.plot(IR, "o", label="fake data")
        fl.plot(f.eval(100), label="fit", human_units=False)
        T1 = 1.0 / f.output("R_1")
        fl.show()
        quit()
        # }}}
        # }}}
        # {{{save to hdf5 file
        s.name(
            searchstr
        )  # use searchstr as the node name withing the HDF5 file
        s.hdf5_write(h5_name, directory=getDATADIR(h5_dir))
        logger.info("saving as hdf5 file with shape:", strm(ndshape(s)))
        # }}}
    # if read_h5 is True there must already be a written h5 file for the file that
    # is to be plotted/imaged that already exists.
    if read_h5:
        s = nddata_hdf5(h5_name + "/" + searchstr, getDATADIR(exp_type=h5_dir))
        # {{{attempting ILT plot with NNLS_Tikhonov_190104
        tE_axis = s.getaxis("tE")
        vd_list = s.getaxis("indirect")
        Nx = 50
        Ny = 50
        Nx_ax = nddata(np.logspace(-5, 3, Nx), "T1")
        Ny_ax = nddata(np.logspace(-5, 3, Ny), "T2")
        s_ILT = s.C
        s_ILT.rename("indirect", "tau1").setaxis("tau1", vd_list)
        s_ILT.rename("tE", "tau2").setaxis("tau2", tE_axis)
        s_ILT = s_ILT.C.nnls(
            ("tau1", "tau2"),
            (Nx_ax, Ny_ax),
            (
                lambda x1, x2: 1.0 - 2 * np.exp(-x1 / x2),
                lambda y1, y2: np.exp(-y1 / y2),
            ),
            l="BRD",
        )
        s_ILT.setaxis("T1", np.log10(Nx_ax.data)).set_units("T1", None)
        s_ILT.setaxis("T2", np.log10(Ny_ax.data)).set_units("T2", None)
        # figure()
        # title(r'$T_{1} - T_{2} distribution$ for free AT in soln')
        fl.next("%s" % searchstr)
        fl.image(s_ILT)
        plt.xlabel(r"$log(T_2/$s$)$")
        plt.ylabel(r"$log(T_1/$s$)$")
fl.show()
