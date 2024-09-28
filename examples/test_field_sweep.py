import pyspecdata as psd
import pyspecProcScripts as prscr

signal_pathway = {"ph1": 1}
with psd.figlist_var() as fl:
    thisfile, exptype, nodename, post_proc, label_str = (
        "240924_13p5mM_TEMPOL_field.h5",
        "ODNP_NMR_comp/field_dependent",
        "field_1",
        "field_sweep_v4",
        "240924 13.5 mM TEMPOL field sweep",
    )
    s = psd.find_file(
        thisfile,
        exp_type=exptype,
        expno=nodename,
        postproc=post_proc,
        lookup=prscr.lookup_table,
    )
    use_freq = False
    if use_freq:
        s["indirect"] = s["indirect"]["carrierFreq"]*1e6
        s.set_units("indirect","Hz")
    else:
        # I wanted to use the carrier (use_frep=True), but it seems like the
        # first point isn't stored properly
        s["indirect"] = s["indirect"]["Field"]
    prscr.rough_table_of_integrals(s, fl=fl)
