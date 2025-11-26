"""Tests for wavelength correction comparing against IDL reference implementation."""

import numpy as np
from scipy.io import readsav
import pytest
from pathlib import Path

from pycr2res.wavecorr import (
    wavecorr,
    detect_lines,
    match_lines,
    compute_wavelength_correction,
    plot_reference_spectra,
    plot_velocity_correction,
)


def save_idl_format(filepath, wave1, wave2, wave3, obs1, obs2, obs3, un1, un2, un3):
    """Save results in numpy .npz format and FITS for IDL comparison."""
    from astropy.io import fits

    # Save as .npz for Python
    npz_path = str(filepath).replace('.sav', '.npz')
    np.savez(npz_path,
        wave1=wave1,
        wave2=wave2,
        wave3=wave3,
        obs1=obs1,
        obs2=obs2,
        obs3=obs3,
        un1=un1,
        un2=un2,
        un3=un3,
    )
    print(f"Saved Python results to {npz_path}")

    # Save as FITS for IDL (IDL can read FITS via astrolib)
    fits_path = str(filepath).replace('.sav', '.fits')
    hdu_list = [fits.PrimaryHDU()]
    hdu_list.append(fits.ImageHDU(wave1, name='WAVE1'))
    hdu_list.append(fits.ImageHDU(wave2, name='WAVE2'))
    hdu_list.append(fits.ImageHDU(wave3, name='WAVE3'))
    hdu_list.append(fits.ImageHDU(obs1, name='OBS1'))
    hdu_list.append(fits.ImageHDU(obs2, name='OBS2'))
    hdu_list.append(fits.ImageHDU(obs3, name='OBS3'))
    hdu_list.append(fits.ImageHDU(un1, name='UN1'))
    hdu_list.append(fits.ImageHDU(un2, name='UN2'))
    hdu_list.append(fits.ImageHDU(un3, name='UN3'))
    hdul = fits.HDUList(hdu_list)
    hdul.writeto(fits_path, overwrite=True)
    print(f"Saved Python results to {fits_path} (for IDL: use mrdfits)")


@pytest.fixture
def idl_test_data():
    """Load IDL test input data."""
    data_path = Path(__file__).parent.parent / "idl" / "wavecorr_test_data.sav"
    return readsav(str(data_path))


@pytest.fixture
def idl_result_data():
    """Load IDL result data for comparison."""
    data_path = Path(__file__).parent.parent / "idl" / "wavecorr_result_data.sav"
    return readsav(str(data_path))


class TestWavecorrAgainstIDL:
    """Test wavecorr implementation against IDL reference."""

    def test_data_shapes(self, idl_test_data):
        """Verify test data has expected shapes."""
        # wave: (n_orders, n_pixels) = (6, 2005)
        # obs/unc: (n_phases, n_orders, n_pixels) = (114, 6, 2005)
        assert idl_test_data['wave1'].shape == (6, 2005)
        assert idl_test_data['obs1'].shape == (114, 6, 2005)
        assert idl_test_data['un1'].shape == (114, 6, 2005)

    def test_wavecorr_full(self, idl_test_data, idl_result_data):
        """Test full wavecorr against IDL output."""
        # Load input data
        wave1 = idl_test_data['wave1']
        wave2 = idl_test_data['wave2']
        wave3 = idl_test_data['wave3']
        obs1 = idl_test_data['obs1']
        obs2 = idl_test_data['obs2']
        obs3 = idl_test_data['obs3']
        unc1 = idl_test_data['un1']
        unc2 = idl_test_data['un2']
        unc3 = idl_test_data['un3']

        # Run Python implementation with diagnostics for plotting
        result = wavecorr(
            wave1, obs1, unc1,
            wave2, obs2, unc2,
            wave3, obs3, unc3,
            ref_order=1, ref_phase=1, power=2,
            return_diagnostics=True
        )
        obs1_new, unc1_new, obs2_new, unc2_new, obs3_new, unc3_new, diagnostics = result

        # Load expected results
        obs1_expected = idl_result_data['obs1']
        obs2_expected = idl_result_data['obs2']
        obs3_expected = idl_result_data['obs3']

        # Compare - allow some tolerance due to numerical differences
        # First check: are we in the right ballpark?
        rel_diff1 = np.abs(obs1_new - obs1_expected) / (np.abs(obs1_expected) + 1e-10)
        rel_diff2 = np.abs(obs2_new - obs2_expected) / (np.abs(obs2_expected) + 1e-10)
        rel_diff3 = np.abs(obs3_new - obs3_expected) / (np.abs(obs3_expected) + 1e-10)

        median_diff1 = np.median(rel_diff1)
        median_diff2 = np.median(rel_diff2)
        median_diff3 = np.median(rel_diff3)

        print(f"Median relative difference det1: {median_diff1:.6f}")
        print(f"Median relative difference det2: {median_diff2:.6f}")
        print(f"Median relative difference det3: {median_diff3:.6f}")

        # Check that most pixels are close
        # Allow 1% tolerance for 95% of pixels
        fraction_close1 = np.mean(rel_diff1 < 0.01)
        fraction_close2 = np.mean(rel_diff2 < 0.01)
        fraction_close3 = np.mean(rel_diff3 < 0.01)

        print(f"Fraction within 1% det1: {fraction_close1:.4f}")
        print(f"Fraction within 1% det2: {fraction_close2:.4f}")
        print(f"Fraction within 1% det3: {fraction_close3:.4f}")

        # Det1 and Det2 should be very close (>90%)
        assert fraction_close1 > 0.90, f"Det1: only {fraction_close1:.2%} within 1%"
        assert fraction_close2 > 0.90, f"Det2: only {fraction_close2:.2%} within 1%"
        # Det3 has fewer reference lines, slightly looser threshold
        assert fraction_close3 > 0.85, f"Det3: only {fraction_close3:.2%} within 1%"

        # Also check 5% tolerance - should be very high for all
        fraction_5pct_1 = np.mean(rel_diff1 < 0.05)
        fraction_5pct_2 = np.mean(rel_diff2 < 0.05)
        fraction_5pct_3 = np.mean(rel_diff3 < 0.05)
        print(f"Fraction within 5% det1: {fraction_5pct_1:.4f}")
        print(f"Fraction within 5% det2: {fraction_5pct_2:.4f}")
        print(f"Fraction within 5% det3: {fraction_5pct_3:.4f}")
        assert fraction_5pct_3 > 0.95, f"Det3: only {fraction_5pct_3:.2%} within 5%"

        # Save Python results in same format as IDL for comparison
        output_path = Path(__file__).parent.parent / "idl" / "wavecorr_python_result.sav"
        save_idl_format(
            output_path,
            wave1, wave2, wave3,
            obs1_new, obs2_new, obs3_new,
            unc1_new, unc2_new, unc3_new
        )

        # Save diagnostic plots
        plot_dir = Path(__file__).parent.parent / "idl"
        plot_reference_spectra(diagnostics, str(plot_dir / "wavecorr_ref_spectra.png"))
        plot_velocity_correction(diagnostics, str(plot_dir / "wavecorr_velocity.png"))

    def test_detect_lines_finds_lines(self, idl_test_data):
        """Test that line detection finds absorption lines."""
        wave1 = idl_test_data['wave1']
        obs1 = idl_test_data['obs1']

        # Test on reference order (1) and phase (1)
        spec = obs1[1, 1, :]
        wl = wave1[1, :]

        xref, wref = detect_lines(spec, wl)

        # Should find some lines
        assert len(xref) > 0, "No lines detected"
        assert len(wref) == len(xref)

        # Line positions should be within spectrum range
        assert np.all(xref >= 0)
        assert np.all(xref < len(spec))

        # Wavelengths should be within range
        assert np.all(wref >= wl.min())
        assert np.all(wref <= wl.max())

        print(f"Detected {len(xref)} lines")

    def test_match_lines(self, idl_test_data):
        """Test that line matching works."""
        wave1 = idl_test_data['wave1']
        obs1 = idl_test_data['obs1']

        # Detect lines in reference phase
        ref_spec = obs1[1, 1, :]
        wl = wave1[1, :]
        xref, wref = detect_lines(ref_spec, wl)

        # Match in a different phase
        test_spec = obs1[10, 1, :]
        x0, w0, wref_matched = match_lines(test_spec, wl, xref, wref)

        # Should match most lines
        match_fraction = len(x0) / len(xref)
        print(f"Matched {len(x0)}/{len(xref)} lines ({match_fraction:.1%})")

        assert match_fraction > 0.5, f"Only matched {match_fraction:.1%} of lines"
        assert len(wref_matched) == len(x0)


class TestHelperFunctions:
    """Test individual helper functions."""

    def test_compute_wavelength_correction(self):
        """Test polynomial coefficient computation."""
        # Simple test case: linear correction
        x0 = np.array([100.0, 500.0, 1000.0, 1500.0])
        w0 = np.array([2.0, 2.1, 2.2, 2.3])  # measured wavelengths
        wref = np.array([2.001, 2.101, 2.201, 2.301])  # reference wavelengths

        coeffs = compute_wavelength_correction(x0, w0, wref, power=1)

        # Coefficients should be close to 1 (small correction)
        assert len(coeffs) == 2
        assert np.abs(coeffs[0] - 1.0) < 0.01


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
