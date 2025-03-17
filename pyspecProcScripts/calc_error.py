import pyspecdata as psd
import numpy as np
from .simple_functions import select_pathway

def calc_error(s,int_bounds=(-500,500), signal_pathway = {"ph1":1}, direct= "t2",indirect="nScans",use_mask = False, error_path = {}):
    f = s.getaxis(direct)
    df = f[1]-f[0]
    collected_variance = psd.ndshape(
        [psd.ndshape(s)[indirect], len(error_path)], [indirect, "pathways"]
    ).alloc().setaxis("pathways",error_path)
    if use_mask == False:
        assert error_path is not None, "If you don't want to use the mask you need to tell me what pathways to use for calculating the noise!"
        for j in range(len(error_path)):
            # calculate N₂ Δf² σ², which is the variance of the integral
            # (by error propagation) where N₂ is the number of points in
            # the indirect dimension
            s_forerror = select_pathway(s, error_path[j])
            # previous line wipes everything out and starts over -- why not use
            # collected_variance above, as I had originally set up --> part of
            # issue #44
            if j == 0:
                N2 = psd.ndshape(s_forerror)[direct]
            # mean divides by N₁ (indirect), integrate multiplies by Δf, and the
            # mean sums all elements (there are N₁N₂ elements)
            s_forerror -= s_forerror.C.mean_all_but([indirect, direct]).mean(
                direct
            )
            s_forerror.run(lambda x: abs(x) ** 2 / 2).mean_all_but(
                [direct, indirect]
            ).mean(direct)
            s_forerror *= df**2  # Δf
            s_forerror *= N2
            collected_variance["pathways",j] = s_forerror
        collected_variance = psd.sqrt(collected_variance.sum("pathways")/len(error_path))
        return collected_variance
