# PyCr2res

Additional recipes for the [CRIRES+
pipeline](https://www.eso.org/sci/software/pipelines/cr2res/), built on ESO's
pyesorex framework.

## Algorithms

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
