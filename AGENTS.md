# Repository Instructions

- When running locally: Before running tests or making code changes, always use the shared base environment at `/Users/atahan/.venv` by running `source /Users/atahan/.venv/bin/activate`. Do not reinstall the project locally unless explicitly asked; the local checkout is already installed in editable mode.
- When running on cloud: Before running tests or making code changes, always install the project in editable mode with Meson by running `pip install -e . --no-build-isolation` from the repository root.
- When running locally: Do not create or switch to any other virtual environment for this repository.
- Follow the linting and validation rules enforced by the workflows and scripts in `.github` when making changes.
- On cloud: Ensure all runtime and testing dependencies are installed so that `pytest` can execute without missing modules. If the editable install fails because build tools are absent (e.g., Meson, ninja, numpy), install the missing packages and rerun the command.
- On cloud: After installing dependencies, run the relevant `pytest` targets to validate changes.
- Locally: In the base environment, run the relevant `pytest` targets to validate changes.

## Project coding rules

- Prefer `pyspecdata` idioms and objects over raw array-level replacements.
- Use the most up-to-date scripts in the repository as implementation references.
- Keep code as legible as possible.
- Use self-explanatory variable names.
- Avoid defining intermediate variables unless they are reused or materially improve clarity.
- Do not introduce helper functions when they are unnecessary.
- Avoid helper functions when the logic is used only once.
- Avoid helper functions when the operation is simple.
- Inline simple one-use logic when that keeps the code clearer.
- Helper functions are acceptable when reuse, complexity, testing, or clarity justifies them.
- Avoid `try`/`except` blocks unless genuinely necessary.
- Prefer elegant, optimized, direct implementations where possible.
- Keep diffs minimal to make PR review and conflict resolution easier.
- Comments are encouraged when they explain implementation logic, signal-processing reasoning, or repo conventions.
- Update comments when behavior changes.
- Update docstrings when behavior, API, or physical interpretation changes.
- Do not remove useful tests to make code simpler.
- Do not silently weaken tests when they encode desired behavior.
- Preserve intended test coverage.

## Branch and workflow rules

- User-provided branch/workflow instructions take priority over later style preferences.
- Every implementation step, branch merge, and commit requires separate user approval.
- Do not mutate a control/reference branch just to create a baseline.
- If a control/reference branch is broken, repair it only in a separately approved commit.
- Repairs to control/reference branches should include tests, comments, and docstrings when appropriate.
- Use the destination branch for final integration work.
- Do not merge or cherry-pick unrelated changes wholesale from reference branches.
- Consult reference branches, commits, and scripts for behavior without copying unrelated changes.

## Testing rules

- Include testing steps in every implementation plan.
- Behavioral changes should include corresponding tests.
- API changes should include corresponding tests.
- Prefer focused tests for the changed behavior first.
- Run broader validation after focused tests when risk justifies it.
- Do not remove useful tests to make code simpler.
- Do not silently weaken tests when they encode desired behavior.
- If a test requires a small patch point, keep the production code minimal and mark the helper with `SINGLE_USE_EXCEPTION`.
- For unit parsing/formatting changes, verify both display behavior and algebraic equivalence.

## Dependency and CI rules

- CI installs the package with `pip install .`.
- CI then installs `pytest`.
- CI runs `pytest`.
- Runtime dependencies imported by installed package modules must be declared in `pyproject.toml` under `[project].dependencies`.
- Do not rely on locally installed packages unless they are declared dependencies or explicitly installed by CI.

## Flake8 and formatting rules

- Flake8 runs only on changed Python files.
- `__init__.py` files are excluded from Flake8 changed-file checks.
- `demos/` files are excluded from Flake8 changed-file checks.
- Flake8 uses `--max-line-length=79`.
- Flake8 uses `--inline-quotes='"'`.
- Flake8 ignores `E203,E225,E226,E231,E261,E265,E293,E302,E305,E306,E401,E731,E741,Q002,W291,W293,W503`.
- Prefer keeping changed lines at or below 79 characters.
- `black` line length is 79.
- `ruff` line length is 79.
- `ruff` extends `E501`, so line length is explicitly checked when Ruff is used.

## Single-use checker rules

- The single-use checker runs only on changed Python files.
- `__init__.py` files are excluded from the single-use checker.
- The single-use checker flags module-level variables that are defined once and used once.
- The single-use checker flags top-level functions that are defined once and directly called once.
- The preferred fix for a flagged single-use variable is to inline the expression at its only use.
- The preferred fix for a flagged single-use top-level function is to move the logic to the call site.
- Do not introduce helper functions for simple one-use logic unless there is a clear reason.
- If a single-use helper or variable is intentionally kept, place `# SINGLE_USE_EXCEPTION -- <concise reason>` immediately above its definition.
- Valid `SINGLE_USE_EXCEPTION` reasons include public API, exported API, downstream compatibility, test patch point, dispatch boundary, and control-flow clarity.
- Tuple/list unpacking variables are allowed by the single-use checker.
- `with ... as name` aliases are allowed by the single-use checker.
- `for` loop targets are allowed by the single-use checker.
- `main()` called only under `if __name__ == "__main__":` is allowed by the single-use checker.
- Single-use module variables are allowed inside a top-of-file vim fold block for editable parameters.
- Editable-parameter fold blocks must start with `# {{{`.
- Editable-parameter fold blocks must end with `# }}}`.
- Editable-parameter fold block comments must include the phrase `changeable parameters`.

## TODO policy

- Changed files must not contain the exact marker `TODO ☐`.
- Use a different wording or resolve the TODO before pushing.
