from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


def _load_impl():
    module_path = Path(__file__).with_name("pyspecProcScripts") / "hack_acq_params.py"
    spec = spec_from_file_location("_pyspecProcScripts_hackacq_impl", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load hack-acq module from {module_path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


main = _load_impl().main


if __name__ == "__main__":
    raise SystemExit(main())
