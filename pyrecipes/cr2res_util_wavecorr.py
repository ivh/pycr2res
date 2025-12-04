from typing import Any, Dict

import cpl.core
import cpl.ui

from pycr2res.wavecorr import (
    load_fits_sequence,
    save_fits_sequence,
    wavecorr,
    plot_reference_spectra,
    plot_velocity_correction,
)


class WaveCorr(cpl.ui.PyRecipe):
    _name = "cr2res_util_wavecorr"
    _version = "0.2"
    _author = "Thomas Marquart"
    _email = "thomas.marquart@astro.uu.se"
    _copyright = "GPL-3.0-or-later"
    _synopsis = "Shift the wavelength scale to a common reference"
    _description = (
        "In good seeing, an observing sequence with CRIRES will have \n"
        + "wavelength shifts between frames because the AO PSF does not \n"
        + "fill the slit. This recipe shifts all frames and orders in \n"
        + "the sequence to a common reference frame by aligning telluric lines."
    )

    def __init__(self):
        self.parameters = cpl.ui.ParameterList(
            [
                cpl.ui.ParameterValue(
                    name="ref-order",
                    context="cr2res_util_wavecorr",
                    description="Reference order index (0-based) for wavelength alignment",
                    default=5,
                ),
                cpl.ui.ParameterValue(
                    name="ref-phase",
                    context="cr2res_util_wavecorr",
                    description="Reference phase/frame index (0-based)",
                    default=1,
                ),
                cpl.ui.ParameterValue(
                    name="poly-order",
                    context="cr2res_util_wavecorr",
                    description="Polynomial order for wavelength correction",
                    default=2,
                ),
                cpl.ui.ParameterValue(
                    name="window",
                    context="cr2res_util_wavecorr",
                    description="Line fitting window size in pixels",
                    default=21,
                ),
                cpl.ui.ParameterValue(
                    name="filter-width",
                    context="cr2res_util_wavecorr",
                    description="Median filter width for continuum estimation",
                    default=60.0,
                ),
                cpl.ui.ParameterValue(
                    name="plot",
                    context="cr2res_util_wavecorr",
                    description="Generate diagnostic plots",
                    default=False,
                ),
            ]
        )

    def run(
        self, frameset: cpl.ui.FrameSet, settings: Dict[str, Any]
    ) -> cpl.ui.FrameSet:
        # Get parameters
        ref_order = settings.get("ref-order", 5)
        ref_phase = settings.get("ref-phase", 1)
        poly_order = settings.get("poly-order", 2)
        window = settings.get("window", 21)
        filter_width = settings.get("filter-width", 60.0)
        do_plot = settings.get("plot", False)

        # Get input filenames
        filenames = [frame.file for frame in frameset]
        n_files = len(filenames)

        print(f"Processing {n_files} input frames")
        print(f"Reference order: {ref_order}, reference phase: {ref_phase}")
        print(f"Polynomial order: {poly_order}")
        print(f"Line fitting window: {window} pixels, continuum filter width: {filter_width} pixels")
        if do_plot:
            print("Diagnostic plots will be generated")

        if n_files == 0:
            raise ValueError("No input files provided")

        # Load all FITS files into arrays
        print("Loading FITS files...")
        (wave1, obs1, unc1,
         wave2, obs2, unc2,
         wave3, obs3, unc3,
         order_names) = load_fits_sequence(filenames)

        print(f"Loaded {obs1.shape[0]} phases, {obs1.shape[1]} orders, {obs1.shape[2]} pixels")
        print(f"Orders: {order_names}")

        # Validate ref_order
        if ref_order < 0 or ref_order >= len(order_names):
            raise ValueError(
                f"ref-order {ref_order} out of range. "
                f"Valid range: 0-{len(order_names)-1}"
            )

        # Validate ref_phase
        if ref_phase < 0 or ref_phase >= n_files:
            raise ValueError(
                f"ref-phase {ref_phase} out of range. "
                f"Valid range: 0-{n_files-1}"
            )

        # Run wavelength correction
        print("Running wavelength correction...")
        result = wavecorr(
            wave1, obs1, unc1,
            wave2, obs2, unc2,
            wave3, obs3, unc3,
            ref_order=ref_order,
            ref_phase=ref_phase,
            power=poly_order,
            window=window,
            filter_width=filter_width,
            return_diagnostics=do_plot
        )

        if do_plot:
            (obs1_new, unc1_new,
             obs2_new, unc2_new,
             obs3_new, unc3_new,
             diagnostics) = result

            # Generate plots
            print("Generating diagnostic plots...")
            plot_reference_spectra(diagnostics, "wavecorr_ref_spectra.png")
            plot_velocity_correction(diagnostics, "wavecorr_velocity.png")
        else:
            (obs1_new, unc1_new,
             obs2_new, unc2_new,
             obs3_new, unc3_new) = result

        # Save corrected spectra
        print("Saving corrected spectra...")
        output_files = save_fits_sequence(
            filenames,
            obs1_new, unc1_new,
            obs2_new, unc2_new,
            obs3_new, unc3_new,
            order_names,
            output_dir="."
        )

        # Create output frameset
        output_frameset = cpl.ui.FrameSet()
        for output_file in output_files:
            print(f"  {output_file}")
            frame = cpl.ui.Frame(
                file=output_file,
                tag="WAVECORR_SPEC",
                group=cpl.ui.Frame.FrameGroup.PRODUCT
            )
            output_frameset.append(frame)

        print(f"Done. Produced {len(output_files)} output files.")
        return output_frameset
