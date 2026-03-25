# Repository Instructions

- When operating locally in this repository, always use the base environment at `/home/jmfranck/base`.
- Do not create or use a separate virtual environment for local work in this repo.
- Use `/home/jmfranck/base/bin/python` and `/home/jmfranck/base/bin/pip` for installs, editable installs, and verification commands.
- Run editable installs from the repository root with `/home/jmfranck/base/bin/python -m pip install -e . --no-build-isolation`.
- If that editable install fails because a build tool is missing from the base environment, install the missing package into `/home/jmfranck/base` and rerun the same command.
