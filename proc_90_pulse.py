from pyspecdata import *
from pyspecdata.load_files.bruker_nmr import bruker_data
from proc_scripts import zeroth_order_ph,calc_baseline,ph1_real_Abs
from proc_scripts.load_data import postproc_dict
from scipy.optimize import minimize,curve_fit,least_squares
rcParams['lines.linewidth'] = 0.5
matplotlib.rcParams["figure.figsize"] = [8.0,5.0]

fl=figlist_var()
for searchstr, exp_type, which_exp, postproc, manual_phcyc, w0 in [

