"""End-to-end test of the telluric correction recipe against CRIRES-LM.

Downloads (and caches) the original pre-tellcorr extracted spectrum and the
CRIRES-LM tellcorr reference output. Runs pycr2res tellcorr on the raw
input and verifies the resulting wavelength solution agrees with CRIRES-LM's
at sub-km/s precision. The input WL is uncorrected, so this exercises the
full correction (typical pipeline shifts here are 0.6-5 km/s per order).

Reference data is fetched on demand into ``$PYCR2RES_TEST_CACHE`` (default:
``~/.cache/pycr2res-tests/``) and is not committed to the repository. The
test is skipped if the data can't be downloaded (e.g. no network).
"""

from __future__ import annotations

import os
import shutil
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

C_KMS = 299792.458

INPUT_FITS_URL = (
    "https://www.astro.uu.se/crires-lm/files/"
    "bet_Ori_M4368_2024-12-09_02033/"
    "cr2res_obs_nodding_extractedA.fits"
)
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
def crires_lm_input():
    """Pre-tellcorr extracted spectrum with uncorrected WL (cached)."""
    return _download(
        INPUT_FITS_URL, _cache_dir() / "cr2res_obs_nodding_extractedA.fits"
    )


@pytest.fixture(scope="session")
def crires_lm_ref():
    """CRIRES-LM tellcorr output, used as the WL reference (cached)."""
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


@pytest.fixture(scope="module")
def pycr2res_output(crires_lm_input, atmos_dir, tmp_path_factory):
    """Run pycr2res tellcorr on the original raw extraction once per module."""
    from pycr2res.tellcorr import tellcorr

    workdir = tmp_path_factory.mktemp("tellcorr")
    raw = workdir / crires_lm_input.name
    shutil.copy(crires_lm_input, raw)
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

    Input WL is the raw uncorrected scale (typical pipeline shift of
    ~0.6-5 km/s per order), so this measures whether pycr2res
    independently recovers the same solution CRIRES-LM did.
    """
    _, all_dv = _collect_dv(crires_lm_ref, pycr2res_output)
    median_abs = float(np.median(np.abs(all_dv)))
    rms = float(np.std(all_dv))

    print(f"\nmedian |dv| = {median_abs * 1000:.0f} m/s, rms = {rms * 1000:.0f} m/s")
    assert median_abs < 0.1, f"median |dv| = {median_abs:.3f} km/s exceeds 100 m/s"
    assert rms < 0.5, f"dv rms = {rms:.3f} km/s exceeds 500 m/s"


def test_dv_interpolation_fills_failed_orders(crires_lm_ref, pycr2res_output):
    """Orders that fit_segment fails to fit must get a dv-interpolated
    correction from neighboring orders, not be left at the raw input WL.

    For this dataset, 05_01 on all three chips fails the per-order telluric
    fit (BIC prefers no-telluric). The raw input WL on those orders is off
    from CRIRES-LM by ~2.8 km/s; interpolation should walk it back to
    within a couple hundred m/s.
    """
    per_order, _ = _collect_dv(crires_lm_ref, pycr2res_output)
    failed_keys = [(c, "05_01") for c in CHIPS]
    found = [k for k in failed_keys if k in per_order]
    assert found, "expected 05_01 in all chips to be present in output"
    for key in found:
        dv = per_order[key]
        assert np.max(np.abs(dv)) < 0.2, (
            f"{key}: max |dv| = {np.max(np.abs(dv)):.3f} km/s, "
            "dv interpolation likely missed this order"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
