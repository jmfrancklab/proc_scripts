"""(*e.g.*) process an ODNP dataset, process a :sup:`2`\ H
CPMG relaxometry dataset, *etc.*

These are the *top level* routines, because they use
second-level routines.  For the sake of organization,
there are **only three levels** of functions,
and all the top level functions live in ``toplevel.py``.
(The top level functions will replace the current
independent processing scripts).

The idea is that these **only** accept file names and
simple numbers and true-false values are arguments.
They should consist **mostly or entirely** of
second-level functions.
The idea is that you would call these either from a GUI
or from your notebook directory.

Scripts to be converted to top level:
    proc_CPMG.py
    proc_echo_mw.py
    proc_Hahn_echoph.py
    proc_IR.py
    proc_IR_noecho.py
    proc_nutation.py
    proc_square_refl.py
"""
