"""Wavelength correction utilities for CRIRES+ spectroscopy.

This module implements wavelength correction by detecting telluric absorption lines
and computing polynomial shifts to align spectra to a reference.
"""

import numpy as np
from numpy.polynomial import polynomial as P
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit
from typing import Tuple, List
from astropy.io import fits


def median_smooth(arr: np.ndarray, width: float) -> np.ndarray:
    """Apply median filter (equivalent to IDL 'middle' function)."""
    size = int(width)
    if size % 2 == 0:
        size += 1
    # Ensure native byte order for scipy
    arr_native = np.ascontiguousarray(arr, dtype=np.float64)
    return median_filter(arr_native, size=size, mode='reflect')


def gauss_model(x: np.ndarray, amp: float, center: float, sigma: float) -> np.ndarray:
    """Gaussian profile model."""
    return amp * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def lorentz_model(x: np.ndarray, amp: float, center: float, gamma: float) -> np.ndarray:
    """Lorentzian profile model."""
    return amp * gamma**2 / ((x - center)**2 + gamma**2)


def gaussbox_model(x: np.ndarray, amp: float, center: float, sigma: float,
                   box_width: float) -> np.ndarray:
    """Gaussian convolved with box (simplified approximation)."""
    gauss = gauss_model(x, amp, center, sigma)
    return gauss


def fit_line_profile(x: np.ndarray, y: np.ndarray, model: str = 'gauss+lorentz'
                    ) -> Tuple[float, float, np.ndarray]:
    """
    Fit a line profile and return center position and residual.

    Parameters
    ----------
    x : np.ndarray
        Pixel coordinates
    y : np.ndarray
        Flux values (should be positive, i.e., inverted absorption)
    model : str
        Model type: 'gauss', 'lorentz', 'gaussbox', or 'gauss+lorentz' (auto-select)

    Returns
    -------
    center : float
        Fitted center position in pixel coordinates
    residual : float
        Sum of squared residuals
    fit : np.ndarray
        The fitted profile
    """
    if len(x) < 4 or np.all(y <= 0):
        return np.nan, np.inf, y * 0

    # Initial guesses
    amp0 = np.max(y)
    center0 = x[np.argmax(y)]
    sigma0 = len(x) / 4.0

    results = {}

    def try_fit(name, func, p0, bounds):
        try:
            popt, _ = curve_fit(func, x, y, p0=p0, bounds=bounds, maxfev=2000)
            fit = func(x, *popt)
            resid = np.sum((fit - y)**2)
            return popt, resid, fit
        except (RuntimeError, ValueError):
            return None, np.inf, None

    bounds_gauss = ([0, x[0], 0.1], [amp0*10, x[-1], len(x)])
    bounds_lorentz = ([0, x[0], 0.1], [amp0*10, x[-1], len(x)])

    if model == 'gauss' or model == 'gauss+lorentz':
        popt, resid, fit = try_fit('gauss', gauss_model, [amp0, center0, sigma0], bounds_gauss)
        if popt is not None:
            results['gauss'] = (popt[1], resid, fit)

    if model == 'lorentz' or model == 'gauss+lorentz':
        popt, resid, fit = try_fit('lorentz', lorentz_model, [amp0, center0, sigma0], bounds_lorentz)
        if popt is not None:
            results['lorentz'] = (popt[1], resid, fit)

    if model == 'gaussbox':
        bounds_gb = ([0, x[0], 0.1, 0.1], [amp0*10, x[-1], len(x), len(x)])
        popt, resid, fit = try_fit('gaussbox', gaussbox_model, [amp0, center0, sigma0, 2.0], bounds_gb)
        if popt is not None:
            results['gaussbox'] = (popt[1], resid, fit)

    if model == 'gauss+lorentz' and 'gauss' in results and 'lorentz' in results:
        # Also try average of gauss and lorentz (as in IDL code)
        g_center, g_resid, g_fit = results['gauss']
        l_center, l_resid, l_fit = results['lorentz']
        avg_fit = 0.5 * (g_fit + l_fit)
        avg_resid = np.sum((avg_fit - y)**2)
        avg_center = 0.5 * (g_center + l_center)
        results['avg'] = (avg_center, avg_resid, avg_fit)

    if not results:
        return np.nan, np.inf, y * 0

    # Pick best fit (lowest residual)
    best_key = min(results.keys(), key=lambda k: results[k][1])
    return results[best_key]


def spline_interp(x: np.ndarray, y: np.ndarray, x_new: float) -> float:
    """Interpolate using cubic spline (equivalent to IDL spline)."""
    if np.isscalar(x_new):
        x_new_arr = np.array([x_new])
    else:
        x_new_arr = np.asarray(x_new)

    # Use linear interpolation for robustness
    return np.interp(x_new_arr, x, y)[0] if np.isscalar(x_new) else np.interp(x_new_arr, x, y)


def detect_lines(spec: np.ndarray, wl: np.ndarray,
                 window: int = 21, filter_width: float = 60.0,
                 model: str = 'gauss+lorentz'
                ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Auto-detect absorption lines in a spectrum.

    This implements the first mode of wavecorr_ref where lines are automatically
    detected and their positions measured.

    Parameters
    ----------
    spec : np.ndarray
        Spectrum flux values
    wl : np.ndarray
        Wavelength array
    window : int
        Fitting window size in pixels
    filter_width : float
        Median filter width for continuum estimation
    model : str
        Line profile model for fitting

    Returns
    -------
    xref : np.ndarray
        Pixel positions of detected lines
    wref : np.ndarray
        Wavelength values of detected lines
    """
    nx = len(spec)

    # Estimate continuum and find absorption features
    continuum = median_smooth(spec, filter_width)

    # Invert spectrum to make absorption lines positive peaks
    o = spec.copy()
    threshold = continuum - np.std(spec - continuum)
    o[o > threshold] = np.max(spec)
    o = np.max(o) - o

    # Find regions where inverted spectrum is positive
    positive_mask = o > 0
    indices = np.where(positive_mask)[0]

    if len(indices) == 0:
        return np.array([]), np.array([])

    # Find contiguous regions (line candidates)
    gaps = np.where(np.diff(indices) > 1)[0]

    if len(gaps) == 0:
        i1 = [indices[0]]
        i2 = [indices[-1]]
    else:
        i1 = [indices[0]] + [indices[g + 1] for g in gaps]
        i2 = [indices[g] for g in gaps] + [indices[-1]]

    i1 = np.array(i1)
    i2 = np.array(i2)
    ii = (i1 + i2) // 2  # midpoints

    # Filter: require minimum width, not at edges, and significant absorption
    min_width = 8
    valid = (
        (i2 - i1 >= min_width) &
        (i1 > window // 2) &
        (i2 < nx - 1 - window // 2) &
        (spec[ii] < 0.9 * continuum[ii])
    )

    i1 = i1[valid]
    i2 = i2[valid]

    if len(i1) == 0:
        return np.array([]), np.array([])

    # Fit each line
    o_inv = np.max(spec) - spec  # inverted spectrum for fitting

    xref_list = []
    wref_list = []
    dev_list = []

    for j in range(len(i1)):
        xx = np.arange(i1[j], i2[j] + 1)
        yy = o_inv[xx]

        center, resid, _ = fit_line_profile(xx, yy, model=model)

        if np.isfinite(center) and i1[j] <= center <= i2[j]:
            wl_center = spline_interp(xx.astype(float), wl[xx], center)
            xref_list.append(center)
            wref_list.append(wl_center)
            dev_list.append(resid)

    if len(xref_list) == 0:
        return np.array([]), np.array([])

    xref = np.array(xref_list)
    wref = np.array(wref_list)
    dev = np.array(dev_list)

    # Remove outliers based on fit quality
    dev_smooth = median_smooth(dev, 10.0)
    dev_threshold = dev_smooth + np.std(dev - dev_smooth)
    good = dev < dev_threshold

    return xref[good], wref[good]


def match_lines(spec: np.ndarray, wl: np.ndarray,
                xref: np.ndarray, wref: np.ndarray,
                window: int = 21, model: str = 'gauss+lorentz'
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Match reference lines in a new spectrum.

    This implements the second mode of wavecorr_ref where known reference
    line positions are matched in a new spectrum.

    Parameters
    ----------
    spec : np.ndarray
        Spectrum flux values
    wl : np.ndarray
        Wavelength array
    xref : np.ndarray
        Reference pixel positions from detect_lines
    wref : np.ndarray
        Reference wavelength positions from detect_lines
    window : int
        Search/fitting window size in pixels
    model : str
        Line profile model for fitting

    Returns
    -------
    x0 : np.ndarray
        Measured pixel positions in this spectrum
    w0 : np.ndarray
        Measured wavelength positions in this spectrum
    wref_matched : np.ndarray
        Reference wavelengths for matched lines
    """
    o_inv = np.max(spec) - spec  # inverted spectrum
    nx = len(spec)

    x0_list = []
    w0_list = []
    wref_list = []

    for j, (xr, wr) in enumerate(zip(xref, wref)):
        j_center = int(round(xr))
        j1 = j_center - window // 2
        j2 = j_center + window // 2

        if j1 < 0 or j2 >= nx:
            continue

        # Find maximum in window to verify line is present
        jmax = j1 + np.argmax(o_inv[j1:j2+1])

        if abs(jmax - j_center) > window // 4:
            # Line not found at expected position
            continue

        # Fit the line
        i1 = j_center - window // 2
        i2 = j_center + window // 2
        xx = np.arange(i1, i2 + 1)
        yy = o_inv[xx]

        center, _, _ = fit_line_profile(xx, yy, model=model)

        if np.isfinite(center):
            wl_center = spline_interp(xx.astype(float), wl[xx], center)
            x0_list.append(center)
            w0_list.append(wl_center)
            wref_list.append(wr)

    return np.array(x0_list), np.array(w0_list), np.array(wref_list)


def compute_wavelength_correction(x0: np.ndarray, w0: np.ndarray,
                                  wref: np.ndarray, power: int = 1
                                 ) -> np.ndarray:
    """
    Compute polynomial coefficients for wavelength correction.

    The correction is: wl_new = poly(x) * wl

    Parameters
    ----------
    x0 : np.ndarray
        Measured pixel positions
    w0 : np.ndarray
        Measured wavelengths
    wref : np.ndarray
        Reference wavelengths
    power : int
        Polynomial order

    Returns
    -------
    coefficients : np.ndarray
        Polynomial coefficients (numpy polynomial convention: [c0, c1, c2, ...])
    """
    if len(x0) < power + 1:
        return np.array([1.0])  # No correction

    # Fit: wref/w0 as a function of x0
    ratio = wref / w0
    coeffs = np.polyfit(x0, ratio, power)
    # Convert to numpy polynomial convention (low to high degree)
    return coeffs[::-1]


def apply_correction(wl: np.ndarray, spec: np.ndarray, unc: np.ndarray,
                     coeffs: np.ndarray, x_offset: float = 0.0
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply wavelength correction by resampling spectrum.

    Parameters
    ----------
    wl : np.ndarray
        Original wavelength array
    spec : np.ndarray
        Spectrum flux values
    unc : np.ndarray
        Uncertainty values
    coeffs : np.ndarray
        Polynomial coefficients from compute_wavelength_correction
    x_offset : float
        Pixel offset for this detector segment

    Returns
    -------
    spec_new : np.ndarray
        Resampled spectrum on original wavelength grid
    unc_new : np.ndarray
        Resampled uncertainties on original wavelength grid
    """
    from scipy.interpolate import CubicSpline

    nx = len(wl)
    x = np.arange(nx) + x_offset

    # Compute corrected wavelength: wl_new = poly(x) * wl
    correction = P.polyval(x, coeffs)
    wl_corrected = correction * wl

    # Resample back to original wavelength grid using cubic spline
    # IDL spline is cubic spline
    cs_spec = CubicSpline(wl_corrected, spec, extrapolate=True)
    cs_unc = CubicSpline(wl_corrected, unc, extrapolate=True)

    spec_new = cs_spec(wl)
    unc_new = cs_unc(wl)

    return spec_new, unc_new


def load_fits_sequence(filenames: List[str]
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                 np.ndarray, np.ndarray, np.ndarray,
                                 np.ndarray, np.ndarray, np.ndarray,
                                 List[str]]:
    """
    Load a sequence of CRIRES+ FITS files into arrays for wavecorr.

    Parameters
    ----------
    filenames : List[str]
        List of FITS file paths

    Returns
    -------
    wave1, wave2, wave3 : np.ndarray
        Wavelength arrays for each detector, shape (n_orders, n_pixels)
    obs1, obs2, obs3 : np.ndarray
        Observation arrays, shape (n_phases, n_orders, n_pixels)
    unc1, unc2, unc3 : np.ndarray
        Uncertainty arrays, shape (n_phases, n_orders, n_pixels)
    order_names : List[str]
        Order identifiers (e.g., ['02_01', '03_01', ...])
    """
    if not filenames:
        raise ValueError("No input files provided")

    # Read first file to get structure
    with fits.open(filenames[0]) as hdul:
        # Get order names from column names
        chip1_cols = hdul['CHIP1.INT1'].data.dtype.names
        order_names = sorted(set(col.rsplit('_', 1)[0] for col in chip1_cols
                                  if col.endswith('_SPEC')))
        n_orders = len(order_names)
        n_pixels = len(hdul['CHIP1.INT1'].data)

    n_phases = len(filenames)

    # Initialize arrays
    wave1 = np.zeros((n_orders, n_pixels), dtype=np.float64)
    wave2 = np.zeros((n_orders, n_pixels), dtype=np.float64)
    wave3 = np.zeros((n_orders, n_pixels), dtype=np.float64)
    obs1 = np.zeros((n_phases, n_orders, n_pixels), dtype=np.float64)
    obs2 = np.zeros((n_phases, n_orders, n_pixels), dtype=np.float64)
    obs3 = np.zeros((n_phases, n_orders, n_pixels), dtype=np.float64)
    unc1 = np.zeros((n_phases, n_orders, n_pixels), dtype=np.float64)
    unc2 = np.zeros((n_phases, n_orders, n_pixels), dtype=np.float64)
    unc3 = np.zeros((n_phases, n_orders, n_pixels), dtype=np.float64)

    # Load wavelength from first file (assumed same for all)
    with fits.open(filenames[0]) as hdul:
        for i, order in enumerate(order_names):
            wave1[i, :] = hdul['CHIP1.INT1'].data[f'{order}_WL']
            wave2[i, :] = hdul['CHIP2.INT1'].data[f'{order}_WL']
            wave3[i, :] = hdul['CHIP3.INT1'].data[f'{order}_WL']

    # Load spectra from all files
    for phase, filename in enumerate(filenames):
        with fits.open(filename) as hdul:
            for i, order in enumerate(order_names):
                obs1[phase, i, :] = hdul['CHIP1.INT1'].data[f'{order}_SPEC']
                obs2[phase, i, :] = hdul['CHIP2.INT1'].data[f'{order}_SPEC']
                obs3[phase, i, :] = hdul['CHIP3.INT1'].data[f'{order}_SPEC']
                unc1[phase, i, :] = hdul['CHIP1.INT1'].data[f'{order}_ERR']
                unc2[phase, i, :] = hdul['CHIP2.INT1'].data[f'{order}_ERR']
                unc3[phase, i, :] = hdul['CHIP3.INT1'].data[f'{order}_ERR']

    return wave1, obs1, unc1, wave2, obs2, unc2, wave3, obs3, unc3, order_names


def save_fits_sequence(filenames: List[str],
                       obs1: np.ndarray, unc1: np.ndarray,
                       obs2: np.ndarray, unc2: np.ndarray,
                       obs3: np.ndarray, unc3: np.ndarray,
                       order_names: List[str],
                       output_dir: str = "."
                      ) -> List[str]:
    """
    Save corrected spectra back to FITS files.

    Parameters
    ----------
    filenames : List[str]
        Original input FITS file paths (used as templates)
    obs1, obs2, obs3 : np.ndarray
        Corrected observation arrays, shape (n_phases, n_orders, n_pixels)
    unc1, unc2, unc3 : np.ndarray
        Corrected uncertainty arrays, shape (n_phases, n_orders, n_pixels)
    order_names : List[str]
        Order identifiers
    output_dir : str
        Output directory for corrected files

    Returns
    -------
    output_files : List[str]
        List of output file paths
    """
    import os

    output_files = []

    for phase, filename in enumerate(filenames):
        # Generate output filename
        base = os.path.basename(filename)
        name, ext = os.path.splitext(base)
        output_file = os.path.join(output_dir, f"{name}_wavecorr{ext}")

        # Copy original file and modify
        with fits.open(filename) as hdul:
            # Update SPEC and ERR columns for each order
            for i, order in enumerate(order_names):
                hdul['CHIP1.INT1'].data[f'{order}_SPEC'] = obs1[phase, i, :]
                hdul['CHIP1.INT1'].data[f'{order}_ERR'] = unc1[phase, i, :]
                hdul['CHIP2.INT1'].data[f'{order}_SPEC'] = obs2[phase, i, :]
                hdul['CHIP2.INT1'].data[f'{order}_ERR'] = unc2[phase, i, :]
                hdul['CHIP3.INT1'].data[f'{order}_SPEC'] = obs3[phase, i, :]
                hdul['CHIP3.INT1'].data[f'{order}_ERR'] = unc3[phase, i, :]

            # Add history
            hdul[0].header['HISTORY'] = 'Wavelength correction applied by cr2res_util_wavecorr'

            hdul.writeto(output_file, overwrite=True)

        output_files.append(output_file)

    return output_files


def wavecorr(wave1: np.ndarray, obs1: np.ndarray, unc1: np.ndarray,
             wave2: np.ndarray, obs2: np.ndarray, unc2: np.ndarray,
             wave3: np.ndarray, obs3: np.ndarray, unc3: np.ndarray,
             ref_order: int = 1, ref_phase: int = 1, power: int = 2,
             window: int = 21, filter_width: float = 60.0,
             model: str = 'gauss+lorentz', return_diagnostics: bool = False
            ):
    """
    Main wavelength correction procedure for CRIRES+ data.

    Applies wavelength correction across all three detectors and all phases
    by detecting telluric lines in a reference spectrum and aligning all
    other spectra to match.

    Parameters
    ----------
    wave1, wave2, wave3 : np.ndarray
        Wavelength arrays for each detector, shape (n_orders, n_pixels)
    obs1, obs2, obs3 : np.ndarray
        Observation arrays for each detector, shape (n_phases, n_orders, n_pixels)
    unc1, unc2, unc3 : np.ndarray
        Uncertainty arrays for each detector, shape (n_phases, n_orders, n_pixels)
    ref_order : int
        Reference order index (0-based)
    ref_phase : int
        Reference phase index (0-based)
    power : int
        Polynomial order for wavelength correction
    window : int
        Line fitting window size
    filter_width : float
        Median filter width for continuum
    model : str
        Line profile model
    return_diagnostics : bool
        If True, also return diagnostic info for plotting

    Returns
    -------
    obs1_new, unc1_new : np.ndarray
        Corrected observation and uncertainty arrays for detector 1
    obs2_new, unc2_new : np.ndarray
        Corrected arrays for detector 2
    obs3_new, unc3_new : np.ndarray
        Corrected arrays for detector 3
    diagnostics : dict (only if return_diagnostics=True)
        Dictionary with diagnostic info for plotting
    """
    # Make copies to avoid modifying input
    obs1_new = obs1.copy()
    obs2_new = obs2.copy()
    obs3_new = obs3.copy()
    unc1_new = unc1.copy()
    unc2_new = unc2.copy()
    unc3_new = unc3.copy()

    nphase = obs1.shape[0]
    nord = wave1.shape[0]
    nwl1 = wave1.shape[1]
    nwl2 = wave2.shape[1]

    # Pixel coordinate offsets for combined detector array
    x1_offset = 0.0
    x2_offset = float(nwl1)
    x3_offset = float(nwl1 + nwl2)

    # Step 1: Detect reference lines in each detector for the reference phase
    ww1 = wave1[ref_order, :]
    oo1 = obs1[ref_phase, ref_order, :]
    xref1, wref1 = detect_lines(oo1, ww1, window=window, filter_width=filter_width, model=model)

    ww2 = wave2[ref_order, :]
    oo2 = obs2[ref_phase, ref_order, :]
    xref2, wref2 = detect_lines(oo2, ww2, window=window, filter_width=filter_width, model=model)

    ww3 = wave3[ref_order, :]
    oo3 = obs3[ref_phase, ref_order, :]
    xref3, wref3 = detect_lines(oo3, ww3, window=window, filter_width=filter_width, model=model)

    # Step 2: Concatenate reference data across all three detectors
    # Adjust pixel positions for detector offsets
    xref = np.concatenate([xref1, xref2 + nwl1, xref3 + nwl1 + nwl2])
    wref = np.concatenate([wref1, wref2, wref3])

    # Combined wavelength and spectrum arrays for reference order
    w_combined = np.concatenate([wave1[ref_order, :], wave2[ref_order, :], wave3[ref_order, :]])

    # Store coefficients for each phase if diagnostics requested
    all_coeffs = []

    # Step 3: Process each phase
    for jphase in range(nphase):
        # Combined spectrum for this phase
        o_combined = np.concatenate([
            obs1[jphase, ref_order, :],
            obs2[jphase, ref_order, :],
            obs3[jphase, ref_order, :]
        ])

        # Match lines in this spectrum
        x0, w0, matched_wref = match_lines(o_combined, w_combined, xref, wref, window=window, model=model)

        if len(x0) < power + 1:
            # Not enough lines matched, skip this phase
            continue

        # Compute polynomial coefficients
        coeffs = compute_wavelength_correction(x0, w0, matched_wref, power=power)
        all_coeffs.append(coeffs)

        # Step 4: Apply correction to all orders in all detectors
        for jord in range(nord):
            # Detector 1
            spec_new, unc_new_arr = apply_correction(
                wave1[jord, :], obs1[jphase, jord, :], unc1[jphase, jord, :],
                coeffs, x_offset=x1_offset
            )
            obs1_new[jphase, jord, :] = spec_new
            unc1_new[jphase, jord, :] = unc_new_arr

            # Detector 2
            spec_new, unc_new_arr = apply_correction(
                wave2[jord, :], obs2[jphase, jord, :], unc2[jphase, jord, :],
                coeffs, x_offset=x2_offset
            )
            obs2_new[jphase, jord, :] = spec_new
            unc2_new[jphase, jord, :] = unc_new_arr

            # Detector 3
            spec_new, unc_new_arr = apply_correction(
                wave3[jord, :], obs3[jphase, jord, :], unc3[jphase, jord, :],
                coeffs, x_offset=x3_offset
            )
            obs3_new[jphase, jord, :] = spec_new
            unc3_new[jphase, jord, :] = unc_new_arr

    if return_diagnostics:
        diagnostics = {
            'ref_spec1': oo1,
            'ref_spec2': oo2,
            'ref_spec3': oo3,
            'ref_wl1': ww1,
            'ref_wl2': ww2,
            'ref_wl3': ww3,
            'xref1': xref1,
            'xref2': xref2,
            'xref3': xref3,
            'wref1': wref1,
            'wref2': wref2,
            'wref3': wref3,
            'coeffs': all_coeffs,
            'w_combined': w_combined,
        }
        return obs1_new, unc1_new, obs2_new, unc2_new, obs3_new, unc3_new, diagnostics

    return obs1_new, unc1_new, obs2_new, unc2_new, obs3_new, unc3_new


# Speed of light in km/s
C_KMS = 299792.458


def plot_reference_spectra(diagnostics: dict, output_file: str = "wavecorr_ref_spectra.png"):
    """
    Plot reference spectra with detected lines marked.

    Parameters
    ----------
    diagnostics : dict
        Diagnostics dictionary from wavecorr with return_diagnostics=True
    output_file : str
        Output filename for the plot
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=False)

    for i, (ax, det) in enumerate(zip(axes, [1, 2, 3])):
        wl = diagnostics[f'ref_wl{det}']
        spec = diagnostics[f'ref_spec{det}']
        wref = diagnostics[f'wref{det}']

        ax.plot(wl, spec, 'k-', lw=0.5, label=f'Detector {det}')

        # Mark detected lines
        for w in wref:
            ax.axvline(w, color='r', lw=0.8, alpha=0.7)

        ax.set_ylabel('Flux')
        ax.set_title(f'Detector {det}: {len(wref)} lines detected')
        ax.legend(loc='upper right')

    axes[-1].set_xlabel('Wavelength')
    fig.suptitle('Reference Spectra with Detected Lines', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Saved reference spectra plot to {output_file}")


def plot_velocity_correction(diagnostics: dict, output_file: str = "wavecorr_velocity.png"):
    """
    Plot the velocity correction as a function of wavelength.

    Parameters
    ----------
    diagnostics : dict
        Diagnostics dictionary from wavecorr with return_diagnostics=True
    output_file : str
        Output filename for the plot
    """
    import matplotlib.pyplot as plt

    coeffs_list = diagnostics['coeffs']
    w_combined = diagnostics['w_combined']

    if not coeffs_list:
        print("No correction coefficients to plot")
        return

    # Compute velocity correction for each phase
    n_pixels = len(w_combined)
    x = np.arange(n_pixels)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))

    # Plot 1: All phases as light lines, median as bold
    velocities = []
    for coeffs in coeffs_list:
        correction = P.polyval(x, coeffs)
        vel = C_KMS * (correction - 1.0)  # Convert to velocity in km/s
        velocities.append(vel)
        ax1.plot(w_combined, vel, 'b-', alpha=0.1, lw=0.5)

    velocities = np.array(velocities)
    median_vel = np.median(velocities, axis=0)
    ax1.plot(w_combined, median_vel, 'r-', lw=2, label='Median')

    ax1.set_xlabel('Wavelength')
    ax1.set_ylabel('Velocity correction (km/s)')
    ax1.set_title('Velocity Correction vs Wavelength (all phases)')
    ax1.legend()
    ax1.axhline(0, color='gray', ls='--', lw=0.5)

    # Plot 2: Velocity at detector centers vs phase number
    # Get representative wavelengths for each detector
    nwl1 = len(diagnostics['ref_wl1'])
    nwl2 = len(diagnostics['ref_wl2'])
    det_centers = [nwl1 // 2, nwl1 + nwl2 // 2, nwl1 + nwl2 + len(diagnostics['ref_wl3']) // 2]
    det_wls = [w_combined[c] for c in det_centers]

    phases = np.arange(len(coeffs_list))
    for det_idx, (center, wl_center) in enumerate(zip(det_centers, det_wls), 1):
        vel_at_center = [C_KMS * (P.polyval(center, c) - 1.0) for c in coeffs_list]
        ax2.plot(phases, vel_at_center, 'o-', ms=2, lw=1,
                 label=f'Det {det_idx} ({wl_center:.1f} nm)')

    ax2.set_xlabel('Phase number')
    ax2.set_ylabel('Velocity correction (km/s)')
    ax2.set_title('Velocity Correction vs Phase')
    ax2.legend()
    ax2.axhline(0, color='gray', ls='--', lw=0.5)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Saved velocity correction plot to {output_file}")
