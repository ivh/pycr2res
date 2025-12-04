"""Spectrum extraction using CharSlit slit decomposition for CRIRES+ spectroscopy."""

import numpy as np
from typing import Tuple, List, Optional
from astropy.io import fits
import charslit


def load_fits_for_extraction(filename: str, detector: int = 1, order_idx: int = 0
                             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                       np.ndarray, np.ndarray, str]:
    """
    Load a CRIRES+ FITS file and prepare data for slit decomposition.

    Parameters
    ----------
    filename : str
        Path to CRIRES+ FITS file
    detector : int
        Detector number (1, 2, or 3)
    order_idx : int
        Order index (0-based)

    Returns
    -------
    im : np.ndarray
        2D spectral image (nrows, ncols)
    pix_unc : np.ndarray
        Pixel uncertainties (nrows, ncols)
    mask : np.ndarray
        Pixel mask (nrows, ncols), uint8
    wavelength : np.ndarray
        Wavelength array (ncols,)
    order_name : str
        Order identifier (e.g., '05_01')
    """
    with fits.open(filename) as hdul:
        ext_name = f'CHIP{detector}.INT1'

        chip_cols = hdul[ext_name].data.dtype.names
        order_names = sorted(set(col.rsplit('_', 1)[0] for col in chip_cols
                                  if col.endswith('_SPEC')))

        if order_idx >= len(order_names):
            raise ValueError(f"Order index {order_idx} out of range (max {len(order_names)-1})")

        order_name = order_names[order_idx]

        spec_col = f'{order_name}_SPEC'
        err_col = f'{order_name}_ERR'
        wl_col = f'{order_name}_WL'

        wavelength = np.array(hdul[ext_name].data[wl_col], dtype=np.float64)
        im = np.array(hdul[ext_name].data[spec_col], dtype=np.float64)
        pix_unc = np.array(hdul[ext_name].data[err_col], dtype=np.float64)

        # Create mask (0 = good, 1 = bad)
        mask = np.zeros_like(im, dtype=np.uint8)
        mask[~np.isfinite(im)] = 1
        mask[~np.isfinite(pix_unc)] = 1
        mask[pix_unc <= 0] = 1

        # Handle NaNs
        im = np.nan_to_num(im, nan=0.0, posinf=0.0, neginf=0.0)
        pix_unc = np.nan_to_num(pix_unc, nan=1e10, posinf=1e10, neginf=1e10)

    return im, pix_unc, mask, wavelength, order_name


def run_slitdec(im: np.ndarray,
                pix_unc: np.ndarray,
                mask: np.ndarray,
                ycen: Optional[np.ndarray] = None,
                slitcurve: Optional[np.ndarray] = None,
                slitdeltas: Optional[np.ndarray] = None,
                osample: int = 6,
                lambda_sP: float = 0.0,
                lambda_sL: float = 0.1,
                maxiter: int = 20
               ) -> dict:
    """
    Run CharSlit slit decomposition.

    Parameters
    ----------
    im : np.ndarray
        2D spectral image (nrows, ncols)
    pix_unc : np.ndarray
        Pixel uncertainties (nrows, ncols)
    mask : np.ndarray
        Pixel mask (nrows, ncols), uint8
    ycen : np.ndarray, optional
        Order center offsets (ncols,). Default: zeros
    slitcurve : np.ndarray, optional
        Curvature coefficients (ncols, 3). Default: zeros
    slitdeltas : np.ndarray, optional
        Horizontal offsets. Default: uniform spacing
    osample : int
        Oversampling factor (default: 6)
    lambda_sP : float
        Spectrum smoothing (default: 0.0)
    lambda_sL : float
        Slit function smoothing (default: 0.1)
    maxiter : int
        Maximum iterations (default: 20)

    Returns
    -------
    result : dict
        Result dictionary from slitchar.slitdec() containing:
        - spectrum: extracted 1D spectrum (ncols,)
        - slitfunction: slit illumination function (ny,)
        - model: reconstructed 2D image (nrows, ncols)
        - uncertainty: spectrum uncertainties (ncols,)
        - mask: updated mask (nrows, ncols)
        - info: status info array
        - return_code: return code
    """
    nrows, ncols = im.shape

    if ycen is None:
        ycen = np.zeros(ncols, dtype=np.float64)

    if slitcurve is None:
        slitcurve = np.zeros((ncols, 3), dtype=np.float64)

    if slitdeltas is None:
        slitdeltas = np.arange(nrows, dtype=np.float64) - nrows / 2.0

    result = charslit.slitdec(
        im=im,
        pix_unc=pix_unc,
        mask=mask,
        ycen=ycen,
        slitcurve=slitcurve,
        slitdeltas=slitdeltas,
        osample=osample,
        lambda_sP=lambda_sP,
        lambda_sL=lambda_sL,
        maxiter=maxiter
    )

    return result


def save_extraction_results(filename: str,
                            wavelength: np.ndarray,
                            spectrum: np.ndarray,
                            uncertainty: np.ndarray,
                            slitfunction: np.ndarray,
                            model: np.ndarray,
                            order_name: str,
                            detector: int,
                            output_dir: str = "."
                           ) -> str:
    """
    Save extraction results to FITS file.

    Parameters
    ----------
    filename : str
        Original input FITS file path (used as template)
    wavelength : np.ndarray
        Wavelength array (ncols,)
    spectrum : np.ndarray
        Extracted 1D spectrum (ncols,)
    uncertainty : np.ndarray
        Spectrum uncertainties (ncols,)
    slitfunction : np.ndarray
        Slit illumination function (ny,)
    model : np.ndarray
        Reconstructed 2D image (nrows, ncols)
    order_name : str
        Order identifier
    detector : int
        Detector number (1, 2, or 3)
    output_dir : str
        Output directory

    Returns
    -------
    output_file : str
        Path to output FITS file
    """
    import os

    base = os.path.basename(filename)
    name, ext = os.path.splitext(base)
    output_file = os.path.join(output_dir,
                               f"{name}_extracted_chip{detector}_ord{order_name}{ext}")

    hdu_list = []

    primary = fits.PrimaryHDU()
    primary.header['HISTORY'] = 'Spectrum extraction by cr2res_util_newextract'
    primary.header['DETECTOR'] = detector
    primary.header['ORDER'] = order_name
    hdu_list.append(primary)

    cols = [
        fits.Column(name='WAVELENGTH', format='D', array=wavelength),
        fits.Column(name='SPECTRUM', format='D', array=spectrum),
        fits.Column(name='UNCERTAINTY', format='D', array=uncertainty),
    ]
    spec_table = fits.BinTableHDU.from_columns(cols, name='SPECTRUM')
    hdu_list.append(spec_table)

    slit_col = fits.Column(name='SLITFUNCTION', format='D', array=slitfunction)
    slit_table = fits.BinTableHDU.from_columns([slit_col], name='SLITFUNCTION')
    hdu_list.append(slit_table)

    model_hdu = fits.ImageHDU(data=model, name='MODEL')
    hdu_list.append(model_hdu)

    hdul = fits.HDUList(hdu_list)
    hdul.writeto(output_file, overwrite=True)

    return output_file


def extract_spectrum(filename: str,
                    detector: int = 1,
                    order_idx: int = 0,
                    osample: int = 6,
                    lambda_sP: float = 0.0,
                    lambda_sL: float = 0.1,
                    maxiter: int = 20,
                    output_dir: str = "."
                   ) -> str:
    """
    Extract spectrum from CRIRES+ FITS file using CharSlit.

    Parameters
    ----------
    filename : str
        Path to CRIRES+ FITS file
    detector : int
        Detector number (1, 2, or 3)
    order_idx : int
        Order index (0-based)
    osample : int
        Oversampling factor
    lambda_sP : float
        Spectrum smoothing parameter
    lambda_sL : float
        Slit function smoothing parameter
    maxiter : int
        Maximum iterations
    output_dir : str
        Output directory

    Returns
    -------
    output_file : str
        Path to output FITS file
    """
    im, pix_unc, mask, wavelength, order_name = load_fits_for_extraction(
        filename, detector, order_idx
    )

    result = run_slitdec(
        im=im,
        pix_unc=pix_unc,
        mask=mask,
        osample=osample,
        lambda_sP=lambda_sP,
        lambda_sL=lambda_sL,
        maxiter=maxiter
    )

    if result['return_code'] != 0:
        print(f"Warning: slitdec returned code {result['return_code']}")
        print(f"Info: {result['info']}")

    output_file = save_extraction_results(
        filename=filename,
        wavelength=wavelength,
        spectrum=result['spectrum'],
        uncertainty=result['uncertainty'],
        slitfunction=result['slitfunction'],
        model=result['model'],
        order_name=order_name,
        detector=detector,
        output_dir=output_dir
    )

    return output_file
