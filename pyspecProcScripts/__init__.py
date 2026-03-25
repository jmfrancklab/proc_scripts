from importlib import import_module

# Keep package import light so submodule-based console entry points start
# quickly. Public names are loaded the first time they are requested.
_EXPORTS = {
    "zeroth_order_ph": (".phasing", "zeroth_order_ph"),
    "hermitian_function_test": (
        ".phasing",
        "hermitian_function_test",
    ),
    "ph1_real_Abs": (".phasing", "ph1_real_Abs"),
    "determine_sign": (".phasing", "determine_sign"),
    "fid_from_echo": (".phasing", "fid_from_echo"),
    "find_peakrange": (".phasing", "find_peakrange"),
    "rough_table_of_integrals": (
        ".third_level.rough_table_of_integrals",
        "rough_table_of_integrals",
    ),
    "align_esr": (".third_level.align_esr", "align_esr"),
    "lookup_table": (".load_data", "lookup_table"),
    "expand_limits": (".plotting", "expand_limits"),
    "draw_limits": (".plotting", "draw_limits"),
    "fl_mod": (".plotting", "fl_mod"),
    "slice_FID_from_echo": (".slice_FID_from_echo", "slice_FID_from_echo"),
    "calc_baseline": (".baseline", "calc_baseline"),
    "dBm2power": (".Utility", "dBm2power"),
    "Vpp2power": (".Utility", "Vpp2power"),
    "power2dBm": (".Utility", "power2dBm"),
    "center_echo": (".CPMG_phasing", "center_echo"),
    "QESR_apply_scalefactor": (
        ".first_level.QESR_rescale",
        "QESR_apply_scalefactor",
    ),
    "QESR_scalefactor": (
        ".first_level.QESR_rescale",
        "QESR_scalefactor",
    ),
    "show_pickle_table": (
        ".first_level.show_pickle_table",
        "show_pickle_table",
    ),
    "recovery": (".fitting", "recovery"),
    "decay": (".fitting", "decay"),
    "fwhm_calculator": (".fwhm_calculate", "fwhm_calculator"),
    "integrate_limits": (".integrate_limits", "integrate_limits"),
    "correl_align": (".correlation_alignment", "correl_align"),
    "select_pathway": (".simple_functions", "select_pathway"),
    "logobj": (".simple_functions", "logobj"),
    "find_apparent_anal_freq": (
        ".simple_functions",
        "find_apparent_anal_freq",
    ),
    "DCCT": (".DCCT_func", "DCCT"),
    "frequency_domain_integral": (
        ".integral_w_error",
        "frequency_domain_integral",
    ),
    "active_propagation": (
        ".integral_w_error",
        "active_propagation",
    ),
    "apod_matched_filter": (
        ".apod_matched_filter",
        "apod_matched_filter",
    ),
    "QESR": (".third_level.QESR", "QESR"),
    "L2G": (".envelope", "L2G"),
    "fit_envelope": (".envelope", "fit_envelope"),
    "convert_to_power": (".convert_to_power", "convert_to_power"),
    "clock_correct": (".clock_correct", "clock_correct"),
    "calc_masked_variance": (".calc_error", "calc_masked_variance"),
}

__all__ = list(_EXPORTS)


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attr_name)
    globals()[name] = value
    return value
