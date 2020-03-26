The processing will be structured into multiple levels.

Lowest level
============

At the lowest level will be pySpecData, which provides an object-oriented
structure (nddata) for storing all information about code (*e.g.* data, error, axes, units, and names of dimensions)
in a single "container" (technically called an "object"), as well as routines
("methods") for manipulating that data (doing things like Fourier transforming,
taking sums, etc. etc.).

First level in this repo
========================

Currently, we have the following modules that provide functions intended to be reused.  We will refer to these as "first level" routines:

-   :doc:`Phasing module <./phasing>` (zero and first order phasing)
-   :doc:`Alignment module <./align>` (aligning peaks that are drifting in frequency)

Second level in this repo
=========================

We are currently in the process of identifying the most up to date code that does each of the following tasks, and
each will be separated into a function inside an appropriate module, and added to the previous list.
At least parts of this should be relevant
to most things.  We will refer to these as "second level" routines, since they use the "first level" routines:

-  load raw data, and put into appropriately shaped nddata (whether
      Bruker or Spincore) → most up to date in the processing code
      associated with ODNP aquisition, which currently lives in
      ``run_Hahn_echo_mw.py`` in ``spincore_apps/SpinCore_pp``

-  Convert short spin echoes into appropriately sliced and phased FIDs
      → most up to date in ``proc_echo_mw.py``

-  Integrate and phase the integrals → most up to date in ``proc_IR.py``

-  Fit the integrals with a relevant curve → Sam's version of
      ``proc_echo_mw.py``
-  Any potentially reusable code in the Q test code → ``proc_square_refl.py``

-   (planned, see below) :doc:`High level <./high_level>`

Top level
=========

Finally, we can construct higher-level functions, all in a single module, that
utilize these to *e.g.* process an ODNP run.
These are the *top level* routines, because they use second-level routines.
For the sake of organization, there are **only three levels** of functions.

The idea is that these **only** accept file names and simple numbers and true-false values are arguments.
They should consist **mostly or entirely** of second-level functions.
The idea is that you would call these either from a GUI or from your notebook directory.
