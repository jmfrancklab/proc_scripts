# {{{ Imports & Initialization
from pyspecdata import *
import os
from pyspecProcScripts import *
import time
import matplotlib.pyplot as plt

plt.rcParams.update(
    {
        "figure.facecolor": (1.0, 1.0, 1.0, 0.0),
        "axes.facecolor": (1.0, 1.0, 1.0, 0.25),
        "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),
    }
)
fl = fl_mod()
# }}}
# Input parameters
date = time.strftime("%y%m%d")
filename = "210707_Q183R1a_pR_DDM_ODNP"
exp_type = "odnp"
postproc = "spincore_IR_v1"
coherence_pathway = {"ph1": 0, "ph2": 1}
ph2_idx = 1; ph1_idx = 0
for filename, nodename, f_range, t_range, rep, clock_correction, IR, ILT in [
    (
        "210811_S179R1a_pR_DDM_KH2PO4_ODNP",
        "FIR_0dBm",
        (-200, 175),
        (None, 83e-3),
        4,
        True,
        False,
        False,
    ),
]:
    outname = "%s_%s" % (date, "_".join(filename.split("_")[1:-1]))
    titl = "%s %s" % (date, " ".join(filename.split("_")[1:-1]))
    fl.basename = titl
    s = find_file(filename, exp_type=exp_type, expno=nodename,
            postproc=postproc, lookup=lookup_table)

    if 'indirect' in s.dimlabels:
        s.rename('indirect','vd')

    fl.next('\nraw data')
    fl.image(s)
    myslice = s['t2':f_range]

    mysgn = [select_pathway(myslice,coherence_pathway).real['t2',np.argmax(i)].run(np.sign) for i in abs(select_pathway(myslice,coherence_pathway).real.data)][0]
    if len(s.getaxis('vd')) == 12:
        mysgn = nddata(np.r_[-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',s.getaxis('vd'))
    if len(s.getaxis('vd')) == 16:
        mysgn = nddata(np.r_[-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],'vd').labels('vd',s.getaxis('vd'))
    fl.next('\nshowing coherence channel - zoom')
    fl.image(s.C['ph1',ph1_idx]['ph2',ph2_idx]['t2':(-750,750)]*mysgn)
    fl.show()

    # {{{ Process_IR -- adapted from process_IR.py (uses as_scan_nbr and select_pathway from process_IR.py from latest commit on master from 8/25/21)
    # Parameters
    sign = mysgn; signal_pathway = coherence_pathway; excluded_pathways = [(0,0)]
    flip=False; fl = fl; clock_correction = False; 

    # Code
    s *= sign # fixing the sign of each fid along the 'vd' dimension
    s['ph2',0]['ph1',0]['t2':0] = 0 # kill the axial noise
    s.ift('t2') # ift into the time dimension
    s.reorder(['ph1','ph2','vd','t2']) # reordering the data so that the time points
                                       # come last after the variable delay points
    #   Apply DC Offset
    s.ift(['ph1','ph2']) # ift the phase dimensions -- what does this do???
    t_start = t_range[-1] * (3/4) # taking the last time point and multiplying by 3/4                                  # as a starting time point to cut off by for corrections
    rx_offset_corr = s['t2':(t_start,None)].data.mean() # defining the offset correction
                                                   # as the mean along the time axis,
                                                   # given the starting cutoff
    s -= rx_offset_corr # correcting any vertical/signal offset 
    s.ft('t2') # ft back into frequency domain
    s.ft(['ph1','ph2']) #ft back from the ift dimenstion of phase

    #   Zero Crossing
    if flip:
        zero_crossing = abs(select_pathway(s,signal_pathway)).sum('t2').argmin('vd',raw_index=True).item()
        s['vd',:zero_crossing] *= -1

    # Crop Data Along Frequency Domain
    s = s['t2':f_range]
    s.ift('t2') # ift into time domain

    # Clock Correction
    if clock_correction:
        clock_corr = nddata(np.linspace(-3,3,2500),'clock_corr')
        s.ft('t2') # ft into freq domain

        if fl is not None:
            fl.next('\nbefore clock correction')
            fl.image(as_scan_nbr(s))

        s_clock = s['ph1',1]['ph2',0].sum('t2') # sum of all f-domain signal in signal coherence channel
        s.ift(['ph1','ph2'])
        min_index = abs(s_clock).argmin('vd',raw_index=True).item()
        s_clock *= np.exp(-1j * clock_corr * s.fromaxis('vd')) # multiply by exponential containing clock correction
        s_clock['vd', :min_index+1] *= -1
        s_clock.sum('vd').run(abs)

        if fl is not None:
            fl.next('\nclock correction')
            fl.image(s_clock,'.',alpha=0.7)
        clock_corr = s_clock.argmax('clock_corr').item()
        plt.axvline(x=clock_corr, alpha=0.5, color='r')
        s *= np.exp(-1j * clock_corr * s.fromaxis('vd'))
        s.ft(['ph1','ph2'])

        if fl is not None:
            fl.next('\nafter clock correction')
            fl.image(s.C.setaxis('vd','#'))
        s.ift('t2') # ift into time domain

    # Phasing
    best_shift, max_shift = hermitian_function_test(select_pathway(s.C.mean('vd'),signal_pathway))
    quit() 







    #  Error Path
    print(set(excluded_pathways))
    print([(signal_pathway['ph1'],signal_pathway['ph2'])],)
    error_path = (set(((j,k) for j in range(ndshape(s)['ph1']) for k in range(ndshape(s)['ph2'])))
            - set(excluded_pathways)
            - set([(signal_pathway['ph1'],signal_pathway['ph2'])]))
    error_path = [{'ph1':j,'ph2':k} for j,k in error_path]

    # Integrating with Error from Excluded Pathways
    s_int, frq_slice, mystd = integral_w_errors(s, signal_pathway, error_path,
            fl=fl, return_frq_slice=True)
    
    # }}}
