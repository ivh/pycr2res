"""Wavelength correction utilities for CRIRES+ spectroscopy."""

import numpy as np
from typing import Dict, List, Tuple


def get_order_data(
    table, order_num: int, detector: int = 1
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get wavelength, spectrum, and error for a given order.

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



def select_lines(spec: np.ndarray, err: np.ndarray) -> np.ndarray:
    """
    Select spectral lines for wavelength correction.
    """
    # Placeholder - would implement line selection logic
    line_positions = np.array([]) # positoins of selected lines in pixel space
    return line_positions

def shift_wl(wavelength: np.ndarray, polynomial: np.ndarray) -> np.ndarray:
    """
    Shift wavelength array using polynomial coefficients.
    """
    shifted_wl = wavelength + np.polyval(polynomial, wavelength)
    return shifted_wl

def wavecorr_main(table, ref_order: int) -> Tuple[cpl.core.Table, Dict[int, float]]:
    """
    Main function to perform wavelength correction.

    Parameters
    ----------
    table : cpl.core.Table
        FITS table with spectral data
    ref_order : int
        Reference order number for alignment

    Returns
    -------
    outtable : cpl.core.Table
        Table with corrected wavelengths
    polys : Dict[int, float]
        Polynomial coefficients for shifts per order
    """
    outtable = table.copy()
    polys = {}

    # Get reference wavelength and spectrum
    ref_wl, ref_spec, ref_err = get_order_data(table, ref_order)

    # Select lines in reference spectrum
    ref_lines = select_lines(ref_spec, ref_err)

    # Loop over orders to compute shifts and apply corrections
    for order_num in range(2, 10):  # Example for CRIRES+ orders 2-9
        wl, spec, err = get_order_data(table, order_num)

        # Select lines in current spectrum
        curr_lines = select_lines(spec, err)

        # Compute shift (placeholder logic)
        shift = np.median(ref_lines - curr_lines) if len(curr_lines) > 0 else 0.0
        polys[order_num] = shift

        # Apply shift to current spectrum
        shifted_wl, shifted_spec, shifted_err = shift_spectrum(wl, spec, err, shift)

        # Update output table with shifted data
        prefix = f"{order_num:02d}_01"  # Assuming detector 1 for simplicity
        outtable[f"{prefix}_WL"] = shifted_wl
        outtable[f"{prefix}_SPEC"] = shifted_spec
        outtable[f"{prefix}_ERR"] = shifted_err

    return outtable, polys