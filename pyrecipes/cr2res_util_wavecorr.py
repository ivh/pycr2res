from typing import Any, Dict
import os

import numpy as np

import cpl.core
import cpl.ui
from pycr2res.wavecorr import wavecorr_main


class WaveCorr(cpl.ui.PyRecipe):
    _name = "cr2res_util_wavecorr"
    _version = "0.1"
    _author = "Thomas Marquart"
    _email = "thomas.marquart@astro.uu.se"
    _copyright = "GPL-3.0-or-later"
    _synopsis = "Shift the wavelenth scale to a common reference"
    _description = (
        "In good seeing, an observing sequence with CRIRES will have \n"
        + "wavelength shifts between frames because the AO PSF does not \n"
        + "fill the slit. This recipe shifts all frames and orders in \n"
        + "the sequence to a common reference frame."
    )

    def __init__(self):
        self.parameters = cpl.ui.ParameterList(
            [
                cpl.ui.ParameterValue(
                    name="cr2res_util_wavecorr.ref-order",
                    context="cr2res_util_wavecorr",
                    description="Reference order number for wavelength alignment",
                    default=1,
                ),
            ]
        )

    def run(
        self, frameset: cpl.ui.FrameSet, settings: Dict[str, Any]
    ) -> cpl.ui.FrameSet:
        # Get the reference order parameter
        ref_order = settings.get("cr2res_util_wavecorr.ref-order", 1)
        print(f"Reference order: {ref_order}")
        print(f"Processing {len(list(frameset))} input frames")

        # Process the full frameset
        corrected_tables, polys = wavecorr_main(frameset, ref_order)

        # Create output frameset
        output_frameset = cpl.ui.FrameSet()

        # Save corrected spectra to FITS files
        for frame_idx, (original_frame, corrected_table) in enumerate(corrected_tables):
            # Generate output filename based on input
            base_name = os.path.splitext(os.path.basename(original_frame.file))[0]
            output_file = f"{base_name}_wavecorr.fits"

            print(f"Saving corrected spectrum to: {output_file}")

            # Save the corrected table as a FITS file
            corrected_table.save(output_file, mode=cpl.core.Table.SaveMode.CREATE)

            # Add to output frameset
            output_frame = cpl.ui.Frame(
                file=output_file, tag="WAVECORR_SPEC", group=cpl.ui.Frame.FrameGroup.PRODUCT
            )
            output_frameset.append(output_frame)

        # Save polynomial coefficients to a FITS table
        if polys:
            poly_table = self._create_polynomial_table(polys)
            poly_output_file = "wavecorr_polynomials.fits"
            print(f"Saving polynomial coefficients to: {poly_output_file}")

            poly_table.save(poly_output_file, mode=cpl.core.Table.SaveMode.CREATE)

            # Add to output frameset
            poly_frame = cpl.ui.Frame(
                file=poly_output_file,
                tag="WAVECORR_POLY",
                group=cpl.ui.Frame.FrameGroup.PRODUCT,
            )
            output_frameset.append(poly_frame)

        return output_frameset

    def _create_polynomial_table(
        self, polys: Dict[str, np.ndarray]
    ) -> cpl.core.Table:
        """
        Create a CPL table to store polynomial coefficients.

        Parameters
        ----------
        polys : Dict[str, np.ndarray]
            Dictionary mapping frame_order_detector keys to polynomial coefficients

        Returns
        -------
        table : cpl.core.Table
            Table with polynomial coefficients
        """
        # Determine the maximum polynomial degree
        max_degree = max(len(coefs) for coefs in polys.values())

        # Create table with appropriate columns
        table = cpl.core.Table()

        # Create lists for each column
        keys_list = []
        frame_ids = []
        order_nums = []
        detector_nums = []
        coef_arrays = [[] for _ in range(max_degree)]

        for key, coefs in polys.items():
            # Parse key format: "frame_index_order_detector"
            parts = key.split("_")
            frame_id = int(parts[0])
            order_num = int(parts[1])
            detector = int(parts[2])

            keys_list.append(key)
            frame_ids.append(frame_id)
            order_nums.append(order_num)
            detector_nums.append(detector)

            # Store coefficients (pad with zeros if needed)
            for i in range(max_degree):
                coef_arrays[i].append(coefs[i] if i < len(coefs) else 0.0)

        # Add columns to table
        table.set_column_array("KEY", np.array(keys_list, dtype=str))
        table.set_column_array("FRAME_ID", np.array(frame_ids, dtype=int))
        table.set_column_array("ORDER", np.array(order_nums, dtype=int))
        table.set_column_array("DETECTOR", np.array(detector_nums, dtype=int))

        for i in range(max_degree):
            table.set_column_array(f"COEF_{i}", np.array(coef_arrays[i], dtype=float))

        return table
