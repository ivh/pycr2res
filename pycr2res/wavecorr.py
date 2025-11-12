"""Wavelength correction utilities for CRIRES+ spectroscopy."""

import numpy as np
from typing import Dict, List, Tuple


def extract_order_data(
    table, order_num: int, detector: int = 1
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract wavelength, spectrum, and error for a given order.

    Parameters
    ----------
    table : cpl.core.Table
        FITS table containing spectral data
    order_num : int
        Order number (e.g., 2-9 for CRIRES+)
    detector : int
        Detector number (1-3 for CRIRES+)

    Returns
    -------
    wavelength : np.ndarray
        Wavelength array
    spectrum : np.ndarray
        Spectrum flux values
    error : np.ndarray
        Error values
    """
    prefix = f"{order_num:02d}_{detector:02d}"
    wl = np.array(table[f"{prefix}_WL"])
    spec = np.array(table[f"{prefix}_SPEC"])
    err = np.array(table[f"{prefix}_ERR"])
    return wl, spec, err


def cross_correlate_orders(
    ref_wl: np.ndarray,
    ref_spec: np.ndarray,
    target_wl: np.ndarray,
    target_spec: np.ndarray,
) -> float:
    """
    Cross-correlate two spectra to find wavelength shift.

    Parameters
    ----------
    ref_wl : np.ndarray
        Reference wavelength array
    ref_spec : np.ndarray
        Reference spectrum
    target_wl : np.ndarray
        Target wavelength array
    target_spec : np.ndarray
        Target spectrum to align

    Returns
    -------
    shift : float
        Wavelength shift in the same units as input wavelengths
    """
    # Placeholder for actual cross-correlation implementation
    # This would use scipy.signal.correlate or similar
    return 0.0


def select_lines(spec: np.ndarray, threshold: float = 3.0) -> np.ndarray:
    """
    Select spectral lines above a certain threshold.

    Parameters
    ----------
    spec : np.ndarray
        Spectrum array
    threshold : float
        Threshold in units of median absolute deviation

    Returns
    -------
    indices : np.ndarray
        Indices of pixels containing spectral lines
    """
    # Calculate median and MAD for robust statistics
    median = np.median(spec)
    mad = np.median(np.abs(spec - median))

    # Find pixels above threshold
    indices = np.where(np.abs(spec - median) > threshold * mad)[0]

    return indices


def shift_spectrum(
    wavelength: np.ndarray, spectrum: np.ndarray, error: np.ndarray, shift: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Apply wavelength shift to spectrum.

    Parameters
    ----------
    wavelength : np.ndarray
        Original wavelength array
    spectrum : np.ndarray
        Spectrum to shift
    error : np.ndarray
        Error array
    shift : float
        Wavelength shift to apply

    Returns
    -------
    shifted_wl : np.ndarray
        Shifted wavelength array
    shifted_spec : np.ndarray
        Interpolated spectrum on new wavelength grid
    shifted_err : np.ndarray
        Interpolated errors
    """
    shifted_wl = wavelength + shift
    # Placeholder - would interpolate spectrum onto new grid
    return shifted_wl, spectrum, error
