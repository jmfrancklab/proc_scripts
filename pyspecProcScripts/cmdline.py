# PYTHON_ARGCOMPLETE_OK

import argparse
import configparser
import os
import re
import shutil
import sys
from pathlib import Path

from .command_registry import _COMMAND_SPECS, register_command

try:
    import argcomplete
    from argcomplete.completers import SuppressCompleter
except ImportError:
    argcomplete = None
    SuppressCompleter = None

_EXP_TYPE_FILENAME_NODE_HELP = {
    "exp_type": "Configured pyspecdata experiment type.",
    "filename": "HDF5 file name inside the experiment directory.",
    "node": "HDF5 node path to edit.",
}
_ROOT_COMMAND = "pyspecProcScripts"
_AUTOCOMPLETE_CANDIDATES = (
    "python",
    "python3",
    "py",
    "pip",
    "hfss_py",
    "pydifft",
    _ROOT_COMMAND,
)

_RAW_DESCRIPTION = """
Show data with postproc
=======================
`pyspecProcScripts raw EXP_TYPE FILENAME NODENAME`

Fourier transforms (and any needed data corrections for older data) are
performed according to the `postproc_type` attribute of the data node.  This
script plots the result, as well as signal that's averaged along the `nScans`
dimension.

Tested with:

``pyspecProcScripts raw ODNP_NMR_comp/Echoes 240620_200uM_TEMPOL_pm_echo.h5 echo_6``

``pyspecProcScripts raw ODNP_NMR_comp/Echoes 240620_200uM_TEMPOL_pm_generic_echo.h5 \
echo_8``

``pyspecProcScripts raw ODNP_NMR_comp/Echoes 240620_200uM_TEMPOL_pm_generic_CPMG.h5 \
CPMG_9``

``pyspecProcScripts raw ODNP_NMR_comp/field_dependent 240920_27mM_TEMPOL_debug_field \
field_3``

``pyspecProcScripts raw ODNP_NMR_comp/ODNP K42.*A1_kRasbatch240814 ODNP``

``pyspecProcScripts raw ODNP_NMR_comp/ODNP K42.*A1_kRasbatch240814 FIR_34dBm``
"""

_CONFIG = configparser.ConfigParser()
for _config_name in (".pyspecdata", "_pyspecdata"):
    _config_path = Path.home() / _config_name
    if _config_path.exists():
        _CONFIG.read(_config_path, encoding="utf-8")
        break
if _CONFIG.has_section("ExpTypes"):
    _EXP_TYPE_NAMES = sorted(
        [key for key, _ in _CONFIG.items("ExpTypes")],
        key=str.casefold,
    )
    _EXP_TYPE_DIRS = {
        key.casefold(): Path(value).expanduser()
        for key, value in _CONFIG.items("ExpTypes")
    }
else:
    _EXP_TYPE_NAMES = []
    _EXP_TYPE_DIRS = {}
if _CONFIG.has_section("General") and _CONFIG.has_option(
    "General", "data_directory"
):
    _DATA_DIRECTORY = Path(
        _CONFIG.get("General", "data_directory")
    ).expanduser()
else:
    _DATA_DIRECTORY = None

# Cache directory scans and HDF5 top-level nodes so repeated completion stays
# responsive while tabbing through the same path.
_DIRECTORY_CACHE = {}
_NODE_CACHE = {}


class HelpfulArgumentParser(argparse.ArgumentParser):
    """Print full help on parse errors so users can see subcommands/options."""

    def error(self, message):
        self.exit(2, f"{self.prog}: error: {message}\n\n{self.format_help()}")


def _commands_on_path(candidates):
    commands = []
    for name in candidates:
        if name == _ROOT_COMMAND or shutil.which(name):
            commands.append(name)
    return commands


def _build_cli_epilog(include_subcommands=False):
    lines = []
    if include_subcommands:
        lines.append("Available subcommands:")
        for name, spec in _COMMAND_SPECS.items():
            lines.append(f"  {name:<8} {spec['help']}")
        lines.append("")
    interpreter_commands = _commands_on_path(("python", "python3", "py"))
    if not interpreter_commands:
        interpreter_commands = ["python"]
    completion_commands = _commands_on_path(_AUTOCOMPLETE_CANDIDATES)
    lines.extend(
        [
            "AUTOCOMPLETE: install argcomplete and add the following to your "
            "shell rc file (.bashrc, .zshrc, etc.)\nThe syntax for bashrc is:",
            "",
            'eval "$(register-python-argcomplete --shell bash '
            + " ".join(interpreter_commands)
            + ')"',
            "complete -o default -o nospace -F _python_argcomplete "
            + " ".join(completion_commands),
            "",
            "(tailor the list of commands to the different commands you have "
            "available on your computer, and you might need to use AI to refactor for zshrc, etc.)",
        ]
    )
    return "\n".join(lines)


def resolve_exp_type_dir(exp_type):
    if exp_type.casefold() in _EXP_TYPE_DIRS:
        return _EXP_TYPE_DIRS[exp_type.casefold()]
    if _DATA_DIRECTORY is not None:
        return _DATA_DIRECTORY / exp_type
    return None


def resolve_hdf_filename(exp_type, filename, allow_search=True):
    exp_dir = resolve_exp_type_dir(exp_type)
    if exp_dir is None:
        return None
    direct_candidate = Path(filename).expanduser()
    if direct_candidate.is_absolute() and direct_candidate.exists():
        return direct_candidate
    direct_candidate = exp_dir
    for part in Path(filename).parts:
        if part in ("", "."):
            continue
        candidate = direct_candidate / part
        if candidate.exists():
            direct_candidate = candidate
            continue
        try:
            matches = [
                name
                for name in os.listdir(direct_candidate)
                if name.casefold() == part.casefold()
            ]
        except OSError:
            matches = []
        if not matches:
            direct_candidate = None
            break
        direct_candidate = direct_candidate / matches[0]
    if direct_candidate is not None and direct_candidate.exists():
        return direct_candidate
    basename = os.path.basename(filename)
    if basename != filename:
        direct_candidate = exp_dir / basename
        if direct_candidate.exists():
            return direct_candidate
        try:
            matches = [
                name
                for name in os.listdir(exp_dir)
                if name.casefold() == basename.casefold()
            ]
        except OSError:
            matches = []
        if matches:
            return exp_dir / matches[0]
    if not allow_search:
        return None
    try:
        from pyspecdata import search_filename

        return Path(
            search_filename(
                re.escape(basename),
                unique=True,
                exp_type=exp_type,
                print_result=False,
            )
        )
    except Exception:
        return None


def hackacq_completer(prefix, action, parsed_args=None, **kwargs):
    prefix_cf = prefix.casefold()
    if action.dest == "exp_type":
        matches = []
        for value in _EXP_TYPE_NAMES:
            if value.casefold().startswith(prefix_cf):
                matches.append(value)
        return matches
    if parsed_args is None or not getattr(parsed_args, "exp_type", None):
        return []
    if action.dest == "filename":
        base_dir = resolve_exp_type_dir(parsed_args.exp_type)
        if base_dir is None or not base_dir.is_dir():
            return []
        if prefix.endswith("/"):
            prefix_dirname = prefix.rstrip("/")
        else:
            prefix_dirname = os.path.dirname(prefix)
        resolved_parts = []
        if prefix_dirname:
            for part in prefix_dirname.split("/"):
                if not part:
                    continue
                try:
                    directory_mtime = base_dir.stat().st_mtime_ns
                except OSError:
                    return []
                cache_key = str(base_dir)
                if (
                    cache_key not in _DIRECTORY_CACHE
                    or _DIRECTORY_CACHE[cache_key][0] != directory_mtime
                ):
                    with os.scandir(base_dir) as scan_it:
                        _DIRECTORY_CACHE[cache_key] = (
                            directory_mtime,
                            tuple(
                                sorted(
                                    (
                                        (entry.name, entry.is_dir())
                                        for entry in scan_it
                                    ),
                                    key=lambda item: item[0].casefold(),
                                )
                            ),
                        )
                matches = [
                    name
                    for name, is_dir in _DIRECTORY_CACHE[cache_key][1]
                    if is_dir and name.casefold() == part.casefold()
                ]
                if not matches:
                    return []
                resolved_parts.append(matches[0])
                base_dir = base_dir / matches[0]
        try:
            directory_mtime = base_dir.stat().st_mtime_ns
        except OSError:
            return []
        cache_key = str(base_dir)
        if (
            cache_key not in _DIRECTORY_CACHE
            or _DIRECTORY_CACHE[cache_key][0] != directory_mtime
        ):
            with os.scandir(base_dir) as scan_it:
                _DIRECTORY_CACHE[cache_key] = (
                    directory_mtime,
                    tuple(
                        sorted(
                            (
                                (entry.name, entry.is_dir())
                                for entry in scan_it
                            ),
                            key=lambda item: item[0].casefold(),
                        )
                    ),
                )
        candidates = []
        for name, is_dir in _DIRECTORY_CACHE[cache_key][1]:
            relative_name = "/".join(resolved_parts + [name])
            if is_dir:
                candidates.append(relative_name + "/")
            elif name.casefold().endswith((".h5", ".hdf5")):
                candidates.append(relative_name)
        matches = []
        for value in candidates:
            if value.casefold().startswith(prefix_cf):
                matches.append(value)
        return matches
    if action.dest == "node" and getattr(parsed_args, "filename", None):
        try:
            hdf_filename = resolve_hdf_filename(
                parsed_args.exp_type,
                parsed_args.filename,
                allow_search=False,
            )
            if hdf_filename is None or not hdf_filename.exists():
                return []
            file_mtime = hdf_filename.stat().st_mtime_ns
            cache_key = str(hdf_filename)
            if (
                cache_key not in _NODE_CACHE
                or _NODE_CACHE[cache_key][0] != file_mtime
            ):
                import h5py

                with h5py.File(hdf_filename, "r") as hdf_file:
                    _NODE_CACHE[cache_key] = (
                        file_mtime,
                        tuple(sorted(hdf_file.keys(), key=str.casefold)),
                    )
            matches = []
            for value in _NODE_CACHE[cache_key][1]:
                if value.casefold().startswith(prefix_cf):
                    matches.append(value)
            return matches
        except Exception:
            return []
    return []


@register_command(
    "Edit labeling fields inside a saved HDF5 acquisition.",
    help=_EXP_TYPE_FILENAME_NODE_HELP,
)
def hackacq(exp_type, filename, node):
    from .hack_acq_params import run_hackacq

    return run_hackacq(exp_type, filename, node)


@register_command(
    "Show raw data with postproc handling applied.",
    description=_RAW_DESCRIPTION,
    help=_EXP_TYPE_FILENAME_NODE_HELP,
)
def raw(exp_type, filename, node):
    from .raw import run_raw

    return run_raw(exp_type, filename, node)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = HelpfulArgumentParser(
        prog=_ROOT_COMMAND,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_build_cli_epilog(),
    )
    subparsers = parser.add_subparsers(
        dest="command",
        title="available subcommands",
    )
    subparsers.required = True
    for name, spec in _COMMAND_SPECS.items():
        subparser = subparsers.add_parser(
            name,
            help=spec["help"],
            description=spec["description"],
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=_build_cli_epilog(include_subcommands=True),
        )
        for argument in spec["arguments"]:
            action = subparser.add_argument(
                *argument["flags"], **dict(argument["kwargs"])
            )
            if action.dest in {
                "exp_type",
                "filename",
                "node",
            }:
                action.completer = hackacq_completer
        subparser.set_defaults(_handler=spec["handler"])
    if argcomplete is not None:
        # Let argcomplete accept the canonical completion even when the user
        # typed a different case, e.g. `B27` -> `b27/echoes`.
        argcomplete.autocomplete(
            parser,
            always_complete_options=False,
            validator=lambda completion, prefix: completion.casefold().startswith(
                prefix.casefold()
            ),
            default_completer=SuppressCompleter(),
        )
    if not argv:
        parser.print_help()
        return 1
    namespace = parser.parse_args(argv)
    handler = namespace._handler
    handler_kwargs = dict(vars(namespace))
    handler_kwargs.pop("_handler", None)
    handler_kwargs.pop("command", None)
    return handler(**handler_kwargs)


if __name__ == "__main__":
    raise SystemExit(main())
