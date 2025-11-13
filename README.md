# PyCr2res

Additional recipes for the [CRIRES+
pipeline](https://www.eso.org/sci/software/pipelines/cr2res/), built on ESO's
pyesorex framework.

## Algorithms

### cr2res_util_wavecorr

The purpose of this recipe is to correct for wavelength shifts in a time
sequence of observations, that can be either due to instrumental drift or the
fact that the PSF from adaptive optics not always fills the slit. Wavelengths
are made to match a reference order in a reference frame; absolute offset
remains.

**Input**
- SOF with a sequence of reduced spectra, standard table format as produced by
  the cr2res pipeline.
- Parameters
  - _ref-order_, which order number to use as the reference order. 

**Output**
- One FITS for each input spectrum with updated wavelength scales.
- Polynomials that translate old to new wavelength scales.

**Steps**
- Find usable lines in reference order
  - Goal is to use telluric lines only, not stellar lines in target
  - Fit lines in ref-order over the full sequence, reject lines that change in shape
- Use the fits of selected lines for their position in pixel
- Determine low order polynomial in **velocity** that matches the changes in
  line positions over the sequence, taking the reference order as zero point.
- Save polynomials and use them to save spectra with updated wavelength scales.

## Installation

### Prerequisites
- **CPL library** (required): `apt-get install libcpl-dev` on Ubuntu/Debian
- **Python 3.12** (3.13 has compatibility issues)

### Clone and run
```bash
git clone https://github.com/ivh/pycr2res
cd pycr2res
export PYESOREX_PLUGIN_DIR="$(pwd)/pyrecipes"
uv sync
uv run pyesorex --recipes
uv run pyesorex --man-page cr2res_util_wavecorr
```

## File Structure

- `pyrecipes/` - Recipe files for pyesorex discovery
- `pycr2res/` - Python package with shared utilities

See `CLAUDE.md` for development details.
