import inspect


class CommandRegistrationError(Exception):
    """Raised when a subcommand is registered more than once."""


# Keep the registry tiny and explicit so the CLI setup stays easy to follow.
_COMMAND_SPECS = {}


def register_command(
    help_text,
    description=None,
    help=None,
):
    """Register a subcommand handler for the main CLI."""

    def decorator(func):
        name = func.__name__.replace("_", "-")
        if name in _COMMAND_SPECS:
            raise CommandRegistrationError(
                f"Command '{name}' already registered"
            )
        _COMMAND_SPECS[name] = {
            "handler": func,
            "help": help_text.strip(),
            "description": (
                description if description is not None else help_text
            ).strip(),
            "arguments": [],
        }
        signature = inspect.signature(func)
        argument_help = help if help is not None else {}
        for parameter in signature.parameters.values():
            if parameter.kind in (
                inspect.Parameter.VAR_POSITIONAL,
                inspect.Parameter.VAR_KEYWORD,
            ):
                continue
            flags = []
            kwargs = {}
            if parameter.default is inspect._empty:
                flags.append(parameter.name)
            else:
                dash_prefix = "-" if len(parameter.name) == 1 else "--"
                flags.append(dash_prefix + parameter.name.replace("_", "-"))
                kwargs["default"] = parameter.default
                if isinstance(parameter.default, bool):
                    kwargs["action"] = (
                        "store_false" if parameter.default else "store_true"
                    )
                elif parameter.default is not None:
                    kwargs["type"] = type(parameter.default)
            if parameter.annotation not in (inspect._empty, None):
                if "type" not in kwargs and parameter.annotation is not bool:
                    kwargs["type"] = parameter.annotation
            if parameter.name in argument_help:
                kwargs["help"] = argument_help[parameter.name].strip()
            _COMMAND_SPECS[name]["arguments"].append(
                {"flags": flags, "kwargs": kwargs}
            )
        return func

    return decorator
