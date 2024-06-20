"""
Show data with postproc
====================
`py proc_raw_data.py NODENAME FILENAME EXP_TYPE`

Fourier transforms (and any needed data corrections for older data) are performed according to the `postproc_type` attribute of the data node.
This script plots the result, as well as signal that's averaged along the `nScans` dimension.
"""
from pyspecProcScripts.load_data import lookup_table
from pyspecProcScripts import select_pathway
from pyspecdata import *
import sys
import os
assert len(sys.argv) == 4
d = find_file(sys.argv[2], exp_type=sys.argv[3], expno=sys.argv[1],
                      lookup=lookup_table)
with figlist_var() as fl:
    d.squeeze()
    fl.next("raw data")
    def image_or_plot(d):
        if len(d.dimlabels) > 2:
            fl.image(d)
        else:
            fl.plot(d)
    image_or_plot(d)
    if 'nScans' in d.dimlabels:
        d.mean('nScans')
        fl.next("signal averaged along nScans")
        image_or_plot(d)
    if d.get_prop('coherence_pathway') is not None:
        d = select_pathway(d, d.get_prop('coherence_pathway'))
        fl.next("with coherence pathway selected")
        image_or_plot(d)
