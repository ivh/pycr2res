"""Wavelength correction utilities for CRIRES+ spectroscopy."""

import numpy as np
from typing import Dict, List, Tuple
import cpl.core
import cpl.ui


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

def wavecorr_main(
    frameset: cpl.ui.FrameSet, ref_order: int
) -> Tuple[cpl.ui.FrameSet, Dict[str, np.ndarray]]:
    """
    Main function to perform wavelength correction on a sequence of frames.

    Parameters
    ----------
    frameset : cpl.ui.FrameSet
        Input frames with spectral data
    ref_order : int
        Reference order number for alignment

    Returns
    -------
    output_frameset : cpl.ui.FrameSet
        Frameset with corrected wavelengths
    polys : Dict[str, np.ndarray]
        Polynomial coefficients for shifts. Keys are "frame_index_order_detector"
        format, values are polynomial coefficient arrays.
    """
    output_frameset = cpl.ui.FrameSet()
    polys = {}

    # Load all tables from frameset
    tables = []
    frames_list = list(frameset)
    for frame in frames_list:
        table = cpl.core.Table.load(frame.file, 1)
        tables.append((frame, table))

    # Validate that reference order exists in the data
    if tables:
        first_table = tables[0][1]
        available_orders = set()
        for col_name in first_table.column_names:
            if "_WL" in col_name:
                order_num = int(col_name.split("_")[0])
                available_orders.add(order_num)

        if ref_order not in available_orders:
            raise ValueError(
                f"Reference order {ref_order} not found in data. "
                f"Available orders: {sorted(available_orders)}"
            )

    # TODO: Implement cross-frame line selection and polynomial fitting
    # For now, process each frame independently with placeholder logic

    for frame_idx, (frame, table) in enumerate(tables):
        # We modify the table in place since each frame has its own loaded table
        # Get reference wavelength and spectrum from this frame
        ref_wl, ref_spec, ref_err = get_order_data(table, ref_order)

        # Select lines in reference spectrum
        ref_lines = select_lines(ref_spec, ref_err)

        # Loop over orders to compute shifts and apply corrections
        for order_num in range(2, 10):  # Example for CRIRES+ orders 2-9
            for detector in [1]:  # TODO: Extend to detectors 1-3
                try:
                    wl, spec, err = get_order_data(table, order_num, detector)

                    # Select lines in current spectrum
                    curr_lines = select_lines(spec, err)

                    # Compute shift (placeholder logic)
                    shift = (
                        np.median(ref_lines - curr_lines)
                        if len(curr_lines) > 0
                        else 0.0
                    )

                    # Store polynomial (currently just a constant shift)
                    poly_key = f"{frame_idx:03d}_{order_num:02d}_{detector:02d}"
                    polys[poly_key] = np.array([shift, 0.0])  # [constant, linear]

                    # Apply shift to wavelength only (not flux or error)
                    shifted_wl = shift_wl(wl, np.array([shift, 0.0]))

                    # Update table with shifted wavelength (in place)
                    prefix = f"{order_num:02d}_{detector:02d}"
                    table[f"{prefix}_WL"] = shifted_wl
                    # SPEC and ERR remain unchanged

                except Exception:
                    # Skip orders/detectors that don't exist in this frame
                    continue

        # Table has been modified in place, keep the reference
        tables[frame_idx] = (frame, table)

    return tables, polys