# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This project provides Python recipes for the CRIRES+ data reduction pipeline, built on ESO's **pyesorex** framework and **PyCPL** library.

**Project Structure:**
```
pycr2res/
  pycr2res/           # Python package (installed via uv)
    __init__.py
    wavecorr.py       # Wavelength correction utilities
  pyrecipes/          # Recipe files for pyesorex discovery
    cr2res_util_hello.py
    cr2res_util_wavecorr.py
```

**Key Design:** Recipes in `pyrecipes/` are simple CPL interfaces. Business logic goes in the `pycr2res/` package.

## Nomenclature

- To _extract_ means to get the spectrum flux from pixel-space to a 1D array by collapsing one dimension (which does not need to align with pixel grid!). Don't use _extract_ or _extraction_ in other contexts to avoid confusiosn. Find other wordings or synonyms instead.
- A _spectrum_ has three arrays, flux, error and wavelength (wl). _Shifting_ a spectrum does not affect flux or error, only wl, i.e. the spectral bins stay the same, just get assigned new wl. Shift is not just just on offset by usually a low-order polynomial.
- _Resampling_ is when a spectrum's flux gets re-distributed from one wl-scale to another, i.e. the binning changes.
- CRIRES has three _detectors_ that each see a (single-digit) number of _spectral orders_. The spectra are saved separately for each _detector-order_ (aka _segment_) since they are non-contiguous.

## Development Setup

### Prerequisites
**CRITICAL:** This project requires the ESO **CPL library** (Common Pipeline Library) to be installed at the system level:
```bash
# Ubuntu/Debian
apt-get install libcpl-dev

# The build will fail without this!
```

### Python Environment
This project uses **uv** for fast, modern Python package management. Python **3.12** is required (3.13 has matplotlib compatibility issues with pyesorex).

```bash
# Install dependencies (use uv, not pip)
uv sync

# Install with development dependencies
uv sync --all-extras

# Install pre-commit hooks (IMPORTANT: run this once after cloning)
uv run --with pre-commit pre-commit install
```

The `.python-version` file pins Python 3.12 to avoid segfaults.

## Common Commands

### Using uv
**IMPORTANT: Always use `uv run` to execute Python commands.** This ensures the correct environment and dependencies.

```bash
# Run tests
uv run py.test
```

### Running Recipes with pyesorex
**CRITICAL:** Point `PYESOREX_PLUGIN_DIR` to the `pyrecipes/` subdirectory, NOT the repo root:

```bash
# Correct - avoids .venv scanning crash
PYESOREX_PLUGIN_DIR=/home/user/pycr2res/pyrecipes uv run pyesorex --recipes

# Run a recipe with parameters
PYESOREX_PLUGIN_DIR=/home/user/pycr2res/pyrecipes \
  uv run pyesorex cr2res_util_wavecorr \
  --cr2res_util_wavecorr.ref-order=5 \
  test.sof
```

**Why `pyrecipes/` only?** Pyesorex recursively scans directories for Python files. If pointed at repo root, it scans `.venv/`, which causes matplotlib import conflicts and segfaults. The `pyrecipes/` subdirectory solves this.

### Data example

You can get an example FITS file from https://neon.physics.uu.se/crires/examplespec.fits
and make a `test.sof` like this

```bash
echo examplespec.fits TAG > test.sof
```

### Code Quality
```bash
# Format and lint with Ruff 
uv run ruff format .
uv run ruff check .
uv run ruff check --fix .

# Run pre-commit hooks (runs automatically on commit, or manually)
uv run --with pre-commit pre-commit run --all-files

```

## Recipe Development

### PyCPL API Gotchas
When working with `cpl.core.Table` objects:
- Use `len(table)` NOT `table.size()` to get row count
- Use `table.column_names` property NOT `table.get_column_names()`
- Access columns with: `np.array(table["COLUMN_NAME"])`

### PyCPL docs
If needed, look at the PyCPL docs here: https://www.eso.org/sci/software/pycpl/pycpl-site/user/basics.html

### Recipe Structure
```python
from pycr2res.wavecorr import your_function  # Import from package
import cpl.core
import cpl.ui

class YourRecipe(cpl.ui.PyRecipe):
    _name = "recipe_name"
    # ... metadata ...

    def __init__(self):
        self.parameters = cpl.ui.ParameterList([
            cpl.ui.ParameterValue(
                name="recipe_name.param-name",
                context="recipe_name",
                description="...",
                default=value,
            ),
        ])

    def run(self, frameset: cpl.ui.FrameSet, settings: Dict[str, Any]):
        # Get parameters
        param = settings.get("recipe_name.param-name", default)

        # Process frames
        for frame in frameset:
            table = cpl.core.Table.load(frame.file, 1)  # ext number
            # ... do work ...

        return cpl.ui.FrameSet()  # Return products
```

### SOF Files
SOF (Set of Frames) files list input FITS files and their tags:
```
path/to/file.fits TAG_NAME
```

## Build System
- Uses **hatchling** (modern, not setuptools)
- Package config in `pyproject.toml` with `[tool.hatch.build.targets.wheel]`
- Set `package = true` in `[tool.uv]` to enable editable install

## Important Notes
- Do not commit changes without asking unless you are sure this is intended. NEVER push until asked explicitly.
- Before committing, always: `git fetch origin master && git rebase origin/master`
