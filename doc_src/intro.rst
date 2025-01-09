Intro
.....

This set of scripts utilizes
`pySpecData <http://jmfrancklab.github.io/pyspecdata>`_
to perform various signal processing tasks that we use for Overhauser Dynamic Nuclear Polarization.
The examples show how to utilize several functions that include:

-   signal alignment
-   integration

    -   including propagation of associated error

-   phase correction of spin echoes (as in Beaton, Guinness, Franck JCP 2022)
-   several examples of utilizing the pySpecData fitting functionality

    -   nutation curves (SE and FID)
    -   fast inversion recovery
    -   :math:`E(p)` curves for :math:`k_{\sigma}` (as of Oct '24, we are in the process of revamping/simplifying this)

Most of the signal is currently native pySpecData-style HDF5 data, although we
have plans to include some high-field Bruker processing here as well.
