from typing import Any, Dict

import cpl.core
import cpl.ui


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

        # Process each input frame
        for frame in frameset:
            print(f"Processing frame: {frame.file}")
            print(f"  Tag: {frame.tag}")
            print(f"  Group: {frame.group}")

            # Read the FITS table from the frame
            try:
                table = cpl.core.Table.load(frame.file, 1)
                print(f"  Loaded table with {len(table)} rows")
                print(f"  Table columns: {table.column_names}")
            except Exception as e:
                print(f"  Error reading table: {e}")

        return cpl.ui.FrameSet()
