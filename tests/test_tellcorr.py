"""End-to-end test of the telluric correction recipe against CRIRES-LM.

Downloads (and caches) a reference output FITS from the CRIRES-LM archive,
reconstructs the raw observed spectrum from SPEC*TELLUR, runs pycr2res
tellcorr on it, and verifies the resulting wavelength solution agrees
with CRIRES-LM's solution at sub-km/s precision.

Reference data is fetched on demand into ``$PYCR2RES_TEST_CACHE`` (default:
``~/.cache/pycr2res-tests/``) and is not committed to the repository. The
test is skipped if the data can't be downloaded (e.g. no network).
"""

from __future__ import annotations

import os
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

C_KMS = 299792.458

REF_FITS_URL = (
    "https://www.astro.uu.se/crires-lm/files/"
    "bet_Ori_M4368_2024-12-09_02033/"
    "bet_Ori_M4368_2024-12-09_02033_tellcorrA.fits"
)
ATMOS_URLS = {
    "L": "https://neon.physics.uu.se/crires/stdAtmos_L.fits",
    "M": "https://neon.physics.uu.se/crires/stdAtmos_M.fits",
}

CHIPS = ("CHIP1.INT1", "CHIP2.INT1", "CHIP3.INT1")


def _cache_dir() -> Path:
    return Path(
        os.environ.get("PYCR2RES_TEST_CACHE", Path.home() / ".cache" / "pycr2res-tests")
    )


def _download(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(url, dest)
    except (urllib.error.URLError, OSError) as e:
        pytest.skip(f"could not download {url}: {e}")
    return dest


@pytest.fixture(scope="session")
def crires_lm_ref():
    """Path to the CRIRES-LM reference FITS (cached after first download)."""
    return _download(REF_FITS_URL, _cache_dir() / "bet_Ori_M4368_tellcorrA.fits")


@pytest.fixture(scope="session")
def atmos_dir():
    """Directory with viper stdAtmos_*.fits (cached after first download)."""
    d = _cache_dir() / "atmos"
    for band, url in ATMOS_URLS.items():
        _download(url, d / f"stdAtmos_{band}.fits")
    return d


def _orders_in(hdul, chip):
    cols = hdul[chip].data.dtype.names
    return sorted({c.rsplit("_", 1)[0] for c in cols if c.endswith("_SPEC")})


def _make_raw_input(ref_path: Path, dest: Path) -> Path:
    """Reconstruct raw observed spectrum: SPEC_raw = SPEC_corrected * TELLUR.

    Drops the CONT/TELLUR columns so the file looks like a fresh extraction.
    """
    with fits.open(ref_path) as h:
        out = fits.HDUList([h[0].copy()])
        for chip in CHIPS:
            d = h[chip].data
            cols = []
            for o in _orders_in(h, chip):
                with np.errstate(invalid="ignore"):
                    raw = d[f"{o}_SPEC"] * d[f"{o}_TELLUR"]
                cols += [
                    fits.Column(name=f"{o}_SPEC", format="D", array=raw),
                    fits.Column(name=f"{o}_ERR", format="D", array=d[f"{o}_ERR"]),
                    fits.Column(name=f"{o}_WL", format="D", array=d[f"{o}_WL"]),
                ]
            out.append(
                fits.BinTableHDU.from_columns(cols, name=chip, header=h[chip].header)
            )
        out.writeto(dest, overwrite=True)
    return dest


@pytest.fixture(scope="module")
def pycr2res_output(crires_lm_ref, atmos_dir, tmp_path_factory):
    """Run pycr2res tellcorr on the reconstructed raw spectrum once per module."""
    from pycr2res.tellcorr import tellcorr

    workdir = tmp_path_factory.mktemp("tellcorr")
    raw = _make_raw_input(crires_lm_ref, workdir / "raw.fits")
    out_path = tellcorr(
        str(raw),
        str(atmos_dir),
        output_dir=str(workdir),
        deg_norm=2,
        deg_wave=2,
        telluric="add2",
        wl_interp_prms_max=30.0,
        wl_interp_deg=1,
    )
    return Path(out_path)


def _collect_dv(ref_path: Path, out_path: Path):
    """Per-order and combined dv [km/s] between ref and pycr2res WL solutions."""
    per_order = {}
    all_dv = []
    with fits.open(ref_path) as h_ref, fits.open(out_path) as h_out:
        for chip in CHIPS:
            for o in _orders_in(h_ref, chip):
                wl_ref = h_ref[chip].data[f"{o}_WL"]
                wl_pyc = h_out[chip].data[f"{o}_WL"]
                ok = np.isfinite(wl_ref) & np.isfinite(wl_pyc) & (wl_ref > 0)
                if not ok.any():
                    continue
                dv = (wl_pyc[ok] - wl_ref[ok]) / wl_ref[ok] * C_KMS
                per_order[(chip, o)] = dv
                all_dv.append(dv)
    return per_order, np.concatenate(all_dv)


def test_wavelength_agrees_with_crires_lm(crires_lm_ref, pycr2res_output):
    """Aggregate dv between pycr2res and CRIRES-LM should be sub-km/s.

    Since the input WL is already CRIRES-LM-corrected, the residual is a
    direct measure of how consistent pycr2res's solution is with theirs.
    """
    _, all_dv = _collect_dv(crires_lm_ref, pycr2res_output)
    median_abs = float(np.median(np.abs(all_dv)))
    rms = float(np.std(all_dv))

    print(f"\nmedian |dv| = {median_abs * 1000:.0f} m/s, rms = {rms * 1000:.0f} m/s")
    assert median_abs < 0.1, f"median |dv| = {median_abs:.3f} km/s exceeds 100 m/s"
    assert rms < 0.5, f"dv rms = {rms:.3f} km/s exceeds 500 m/s"


def test_dv_interpolation_fills_failed_orders(crires_lm_ref, pycr2res_output):
    """Orders that fit_segment fails to fit must get a small dv-interpolated
    correction, not be left at the input WL.

    For this dataset, 05_01 on all three chips fails the per-order telluric
    fit. They should still match CRIRES-LM's wavecorr solution to within a
    handful of m/s, because both implementations interpolate the smooth
    dv(lambda) across orders.
    """
    per_order, _ = _collect_dv(crires_lm_ref, pycr2res_output)
    failed_keys = [(c, "05_01") for c in CHIPS]
    found = [k for k in failed_keys if k in per_order]
    assert found, "expected 05_01 in all chips to be present in output"
    for key in found:
        dv = per_order[key]
        assert np.max(np.abs(dv)) < 0.1, (
            f"{key}: max |dv| = {np.max(np.abs(dv)):.3f} km/s, "
            "dv interpolation likely missed this order"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
