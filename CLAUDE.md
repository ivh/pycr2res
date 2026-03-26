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
    newextract.py     # Spectrum extraction using CharSlit
    tellcorr.py       # Telluric correction (forward model fitting)
  pyrecipes/          # Recipe files for pyesorex discovery
    cr2res_util_hello.py
    cr2res_util_wavecorr.py
    cr2res_util_newextract.py
    cr2res_util_tellcorr.py
  CharSlit/           # Git submodule for CharSlit library
```

**Key Design:** Recipes in `pyrecipes/` are simple CPL interfaces. Business logic goes in the `pycr2res/` package.

## Nomenclature

- To _extract_ means to get the spectrum flux from pixel-space to a 1D array by collapsing one dimension (which does not need to align with pixel grid!). Don't use _extract_ or _extraction_ in other contexts to avoid confusiosn. Find other wordings or synonyms instead.
- A _spectrum_ has three arrays, flux, error and wavelength (wl). _Shifting_ a spectrum does not affect flux or error, only wl, i.e. the spectral bins stay the same, just get assigned new wl. Shift is not just just on offset by usually a low-order polynomial.
- _Resampling_ is when a spectrum's flux gets re-distributed from one wl-scale to another, i.e. the binning changes.
- CRIRES has three _detectors_ that each see a (single-digit) number of _spectral orders_. The spectra are saved separately for each _detector-order_ (aka _segment_) since they are non-contiguous.

## Development Setup

### Prerequisites
**For CharSlit submodule:** The CharSlit library (used by `cr2res_util_newextract`) requires:
- **CMake** (build system)
- **C++ compiler** with C++17 support
- nanobind and scikit-build-core (automatically installed by uv)

### Python Environment
This project uses **uv** for fast, modern Python package management.

```bash
# Initialize git submodules (IMPORTANT: do this first after cloning)
git submodule update --init --recursive

# Install dependencies (use uv, not pip)
# This will automatically build and install CharSlit from the submodule
uv sync

# Install with development dependencies
uv sync --all-extras

# Install pre-commit hooks (IMPORTANT: run this once after cloning)
uv run --with pre-commit pre-commit install
```


**Note on CharSlit submodule:**
- CharSlit is included as a git submodule in `CharSlit/`
- `uv sync` automatically builds and installs it in editable mode via the local path dependency in `pyproject.toml`
- The build process uses scikit-build-core with CMake, compiling C code with nanobind bindings
- If you update the CharSlit submodule, run `uv sync` again to rebuild

## Common Commands

### Using uv
**IMPORTANT: Always use `uv run` to execute Python commands.** This ensures the correct environment and dependencies.

```bash
# Run tests
uv run python -m pytest

# Run specific test
uv run python -m pytest tests/test_wavecorr.py -v
```

### Running Recipes with pyesorex
**CRITICAL:** Point `PYESOREX_PLUGIN_DIR` to the `pyrecipes/` subdirectory, NOT the repo root:

```bash
# Correct - avoids .venv scanning crash
PYESOREX_PLUGIN_DIR=/home/user/pycr2res/pyrecipes uv run pyesorex --recipes

# Run a recipe with parameters
PYESOREX_PLUGIN_DIR=/home/user/pycr2res/pyrecipes \
  uv run pyesorex cr2res_util_wavecorr \
  --ref-order=5 \
  test.sof
```

**Why `pyrecipes/` only?** Pyesorex recursively scans directories for Python files. If pointed at repo root, it scans `.venv/`, which causes matplotlib import conflicts and segfaults. The `pyrecipes/` subdirectory solves this.

```bash
# Run telluric correction (requires atmospheric model files)
PYESOREX_PLUGIN_DIR=/home/user/pycr2res/pyrecipes \
  uv run pyesorex cr2res_util_tellcorr \
  --atm-data-dir=$VIPER_ATMOS \
  --deg-norm=2 --deg-wave=2 \
  test.sof
```

### Running Recipes from Python
Instead of the shell, recipes can be run programmatically via the Pyesorex API. This is useful for scripting, plotting results, etc.

```python
from pyesorex.pyesorex import Pyesorex

p = Pyesorex()
p.recipe = "cr2res_util_tellcorr"
p.sof_location = "test.sof"
p.recipe_parameters.update({"atm-data-dir": "/path/to/atmos"})
products = p.run()  # returns cpl.ui.FrameSet

for frame in products:
    print(frame.file, frame.tag)  # path and category of each product
    table = cpl.core.Table.load(frame.file, 1)
```

`PYESOREX_PLUGIN_DIR` must still be set (env var or before running the script).

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
                name="param-name",
                context="recipe_name",
                description="...",
                default=value,
            ),
        ])

    def run(self, frameset: cpl.ui.FrameSet, settings: Dict[str, Any]):
        # Get parameters
        param = settings.get("param-name", default)

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

## Testing

### IDL Reference Data
The wavelength correction implementation has been validated against the original IDL code:
- Test data: `idl/wavecorr_test_data.sav` (input spectra)
- Reference results: `idl/wavecorr_result_data.sav` (expected output from IDL)
- Test compares Python output against IDL reference, achieving:
  - Detectors 1-2: >96% of pixels within 1% agreement
  - Detector 3: >85% of pixels within 1% agreement

Tests also generate diagnostic plots:
- `idl/wavecorr_ref_spectra.png` - Reference spectra with detected lines
- `idl/wavecorr_velocity.png` - Velocity corrections vs wavelength

The IDL `.sav` files can be read in Python using `scipy.io.readsav`. For IDL users,
test results are also saved as FITS (`wavecorr_python_result.fits`) which can be
read with `mrdfits` from IDL astrolib.

## Build System
- Uses **hatchling** (modern, not setuptools)
- Package config in `pyproject.toml` with `[tool.hatch.build.targets.wheel]`
- Set `package = true` in `[tool.uv]` to enable editable install

## Important Notes
- Do not commit changes without asking unless you are sure this is intended. NEVER push until asked explicitly.
- Before committing, always: `git fetch origin master && git rebase origin/master`
