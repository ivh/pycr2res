# PyCr2res

Additional recipes for the [CRIRES+
pipeline](https://www.eso.org/sci/software/pipelines/cr2res/), built on ESO's
pyesorex framework.

## Recipes

### cr2res_util_wavecorr

The purpose of this recipe is to correct for wavelength shifts in a time
sequence of observations caused by adaptive optics PSF not consistently filling
the slit. Wavelengths are aligned to a reference order in a reference frame
by detecting and fitting telluric absorption lines.

**Input**
- SOF with a sequence of reduced spectra in standard CRIRES+ table format
  (CHIP1.INT1, CHIP2.INT1, CHIP3.INT1 extensions with columns like `05_01_SPEC`,
  `05_01_ERR`, `05_01_WL` for each order-trace)
- Parameters
  - `--ref-order` (default: 5): Order index (0-based) to use as reference
  - `--ref-phase` (default: 1): Frame index (0-based) to use as reference
  - `--poly-order` (default: 2): Polynomial order for wavelength correction
  - `--window` (default: 21): Line fitting window size in pixels
  - `--filter-width` (default: 60): Median filter width for continuum estimation
  - `--plot` (default: false): Generate diagnostic plots

**Output**
- One FITS file per input spectrum with corrected wavelength scales
  - Spectra are resampled onto the original wavelength grid after correction
  - Both flux (`_SPEC`) and uncertainties (`_ERR`) are corrected
- Optional diagnostic plots (if `--plot=true`):
  - Reference spectra with detected telluric lines marked
  - Velocity correction vs wavelength for all frames

**Algorithm**
1. Detect telluric absorption lines in reference spectrum (ref-order, ref-phase)
   across all three detectors
2. Fit each line with Gaussian/Lorentzian profiles to measure precise positions
3. For each frame, match the reference lines and measure position shifts
4. Fit polynomial correction as function of pixel position: `wl_new = poly(x) * wl`
5. Apply correction to all orders in all detectors via cubic spline resampling

The polynomial coefficients encode a velocity shift (multiplicative factor on
wavelength = Doppler shift). Typical corrections are 0.1-0.5 km/s.

### cr2res_util_tellcorr

Forward-models telluric absorption lines to simultaneously fit atmospheric
transmission, continuum shape, and wavelength calibration for each
detector-order segment. Based on viper (Koehler & Zechmeister).

Requires atmospheric model files (`stdAtmos_*.fits`), e.g. from a
[viper](https://github.com/mzechmeister/viper) installation. Set
`$VIPER_ATMOS` to point to the directory containing them.

**Input**
- SOF with one or more reduced spectra in standard CRIRES+ table format
- Parameters
  - `--atm-data-dir` (**required**): Path to atmospheric model files
  - `--deg-norm` (default: 3): Polynomial degree for continuum normalization
  - `--deg-wave` (default: 3): Polynomial degree for wavelength calibration
  - `--ip` (default: g): Instrumental profile model
  - `--kapsig` (default: 6): Kappa-sigma clipping threshold
  - `--telluric` (default: add): Telluric mode (`add` or `add2`)
  - `--plot` (default: 1): Generate diagnostic PNG

**Output**
- One FITS file per input with the same structure plus:
  - Updated `_WL` columns (fitted wavelength calibration)
  - New `_TELL` columns (telluric transmission model per segment)
  - New `_CONT` columns (continuum/normalization model per segment)
- Diagnostic PNG plot (if `--plot=1`)

**Example**
```bash
PYESOREX_PLUGIN_DIR=pyrecipes uv run pyesorex cr2res_util_tellcorr \
  --atm-data-dir=$VIPER_ATMOS \
  --deg-norm=2 --deg-wave=2 \
  test.sof
```

## Installation

### Prerequisites
- CMake and C++ compiler (only for CharSlit extraction submodule)
- Python packages get installed by uv (or pip) as below.

### Clone and run
```bash
git clone https://github.com/ivh/pycr2res
cd pycr2res

# Initialize submodules (CharSlit library)
git submodule update --init --recursive

# Install dependencies (this builds CharSlit automatically)
export PYESOREX_PLUGIN_DIR="$(pwd)/pyrecipes"
uv sync

# List available recipes
uv run pyesorex --recipes
uv run pyesorex --man-page cr2res_util_wavecorr
```

### Updating the CharSlit submodule

To update CharSlit to the latest version:

```bash
# Update the submodule to latest commit
cd CharSlit
git pull origin master
cd ..

# Rebuild with updated code
uv sync

# Commit the submodule update (updates the commit hash pointer)
git add CharSlit
git commit -m "Update CharSlit submodule"
```

**Why commit?** The pycr2res repo tracks which CharSlit commit to use (just a pointer/hash, not the actual files). When you update CharSlit, you need to commit the new pointer so others get the same version.

## File Structure

- `pyrecipes/` - Recipe files for pyesorex discovery
- `pycr2res/` - Python package with shared utilities

See `CLAUDE.md` for development details.
