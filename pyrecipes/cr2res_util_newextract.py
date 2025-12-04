from typing import Any, Dict

import cpl.core
import cpl.ui

from pycr2res.newextract import extract_spectrum


class NewExtract(cpl.ui.PyRecipe):
    _name = "cr2res_util_newextract"
    _version = "0.1"
    _author = "Thomas Marquart"
    _email = "thomas.marquart@astro.uu.se"
    _copyright = "GPL-3.0-or-later"
    _synopsis = "Extract spectrum using CharSlit slit decomposition"
    _description = (
        "Extract a 1D spectrum from a 2D spectral image using the CharSlit\n"
        + "slit decomposition algorithm. This simultaneously models the spectrum\n"
        + "and the slit illumination function, providing optimal extraction\n"
        + "with built-in outlier rejection."
    )

    def __init__(self):
        self.parameters = cpl.ui.ParameterList(
            [
                cpl.ui.ParameterValue(
                    name="cr2res_util_newextract.detector",
                    context="cr2res_util_newextract",
                    description="Detector number (1, 2, or 3)",
                    default=1,
                ),
                cpl.ui.ParameterValue(
                    name="cr2res_util_newextract.order-idx",
                    context="cr2res_util_newextract",
                    description="Order index (0-based) to process",
                    default=0,
                ),
                cpl.ui.ParameterValue(
                    name="cr2res_util_newextract.osample",
                    context="cr2res_util_newextract",
                    description="Oversampling factor for slit function",
                    default=6,
                ),
                cpl.ui.ParameterValue(
                    name="cr2res_util_newextract.lambda-sp",
                    context="cr2res_util_newextract",
                    description="Spectrum smoothing parameter (0.0 = no smoothing)",
                    default=0.0,
                ),
                cpl.ui.ParameterValue(
                    name="cr2res_util_newextract.lambda-sl",
                    context="cr2res_util_newextract",
                    description="Slit function smoothing parameter",
                    default=0.1,
                ),
                cpl.ui.ParameterValue(
                    name="cr2res_util_newextract.maxiter",
                    context="cr2res_util_newextract",
                    description="Maximum number of iterations",
                    default=20,
                ),
            ]
        )

    def run(
        self, frameset: cpl.ui.FrameSet, settings: Dict[str, Any]
    ) -> cpl.ui.FrameSet:
        detector = settings.get("cr2res_util_newextract.detector", 1)
        order_idx = settings.get("cr2res_util_newextract.order-idx", 0)
        osample = settings.get("cr2res_util_newextract.osample", 6)
        lambda_sP = settings.get("cr2res_util_newextract.lambda-sp", 0.0)
        lambda_sL = settings.get("cr2res_util_newextract.lambda-sl", 0.1)
        maxiter = settings.get("cr2res_util_newextract.maxiter", 20)

        filenames = [frame.file for frame in frameset]
        n_files = len(filenames)

        print(f"Processing {n_files} input frames")
        print(f"Detector: {detector}, Order index: {order_idx}")
        print(f"Oversampling: {osample}, Max iterations: {maxiter}")
        print(f"Smoothing - spectrum: {lambda_sP}, slit function: {lambda_sL}")

        if n_files == 0:
            raise ValueError("No input files provided")

        if detector not in [1, 2, 3]:
            raise ValueError(
                f"detector must be 1, 2, or 3. Got: {detector}"
            )

        if order_idx < 0:
            raise ValueError(
                f"order-idx must be >= 0. Got: {order_idx}"
            )

        output_frameset = cpl.ui.FrameSet()

        for filename in filenames:
            print(f"Processing {filename}...")

            output_file = extract_spectrum(
                filename=filename,
                detector=detector,
                order_idx=order_idx,
                osample=osample,
                lambda_sP=lambda_sP,
                lambda_sL=lambda_sL,
                maxiter=maxiter,
                output_dir="."
            )

            print(f"  Output: {output_file}")

            frame = cpl.ui.Frame(
                file=output_file,
                tag="EXTRACTED_SPEC",
                group=cpl.ui.Frame.FrameGroup.PRODUCT
            )
            output_frameset.append(frame)

        print(f"Done. Produced {len(output_frameset)} output files.")
        return output_frameset
