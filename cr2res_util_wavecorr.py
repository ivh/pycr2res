from typing import Any, Dict

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

    def run(
        self, frameset: cpl.ui.FrameSet, settings: Dict[str, Any]
    ) -> cpl.ui.FrameSet:
        for frame in frameset:
            print(f"Hello, {frame.file}!")
        return cpl.ui.FrameSet()
