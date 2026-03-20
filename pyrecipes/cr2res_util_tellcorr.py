from typing import Any, Dict

import cpl.core
import cpl.ui

from pycr2res.tellcorr import plot_tellcorr, tellcorr


class TellCorr(cpl.ui.PyRecipe):
    _name = "cr2res_util_tellcorr"
    _version = "0.1"
    _author = "Thomas Marquart"
    _email = "thomas.marquart@astro.uu.se"
    _copyright = "GPL-3.0-or-later"
    _synopsis = "Fit and remove telluric absorption features"
    _description = (
        "Forward-models telluric absorption lines to simultaneously fit \n"
        "atmospheric transmission, continuum shape, and wavelength calibration. \n"
        "Outputs a file with updated wavelength scale and additional columns \n"
        "for the telluric and continuum models. \n"
        "Based on viper (Koehler & Zechmeister)."
    )

    def __init__(self):
        self.parameters = cpl.ui.ParameterList(
            [
                cpl.ui.ParameterValue(
                    name="atm-data-dir",
                    context="cr2res_util_tellcorr",
                    description="Path to atmospheric model directory (stdAtmos_*.fits)",
                    default="",
                ),
                cpl.ui.ParameterValue(
                    name="deg-norm",
                    context="cr2res_util_tellcorr",
                    description="Polynomial degree for flux normalization",
                    default=3,
                ),
                cpl.ui.ParameterValue(
                    name="deg-wave",
                    context="cr2res_util_tellcorr",
                    description="Polynomial degree for wavelength calibration",
                    default=3,
                ),
                cpl.ui.ParameterValue(
                    name="deg-bkg",
                    context="cr2res_util_tellcorr",
                    description="Degree of additive background model (0 for none)",
                    default=0,
                ),
                cpl.ui.ParameterValue(
                    name="ip",
                    context="cr2res_util_tellcorr",
                    description="IP model (g: Gaussian, sg: super-Gaussian, ag: asymmetric, bg: biGaussian, mg: multi-Gaussian, mcg: multi-central)",
                    default="g",
                ),
                cpl.ui.ParameterValue(
                    name="kapsig",
                    context="cr2res_util_tellcorr",
                    description="Kappa-sigma clipping threshold (0 to disable)",
                    default=6.0,
                ),
                cpl.ui.ParameterValue(
                    name="tell-bic",
                    context="cr2res_util_tellcorr",
                    description="BIC threshold for telluric model selection",
                    default=10.0,
                ),
                cpl.ui.ParameterValue(
                    name="vcut",
                    context="cr2res_util_tellcorr",
                    description="Velocity cutoff for edge trimming [km/s]",
                    default=100.0,
                ),
                cpl.ui.ParameterValue(
                    name="telluric",
                    context="cr2res_util_tellcorr",
                    description="Telluric mode (add: per-molecule coefficients, add2: combined non-water)",
                    default="add",
                ),
                cpl.ui.ParameterValue(
                    name="molec",
                    context="cr2res_util_tellcorr",
                    description="Molecules to fit (all: auto-select, or comma-separated list)",
                    default="all",
                ),
                cpl.ui.ParameterValue(
                    name="iphs",
                    context="cr2res_util_tellcorr",
                    description="IP half-size in model pixels",
                    default=50,
                ),
                cpl.ui.ParameterValue(
                    name="plot",
                    context="cr2res_util_tellcorr",
                    description="Generate diagnostic plot (0=off, 1=on)",
                    default=1,
                ),
            ]
        )

    def run(
        self, frameset: cpl.ui.FrameSet, settings: Dict[str, Any]
    ) -> cpl.ui.FrameSet:
        atm_data_dir = settings.get("atm-data-dir", "")
        if not atm_data_dir:
            raise ValueError(
                "atm-data-dir must be set to the directory containing "
                "stdAtmos_*.fits atmospheric model files"
            )

        deg_norm = settings.get("deg-norm", 3)
        deg_wave = settings.get("deg-wave", 3)
        deg_bkg = settings.get("deg-bkg", 0)
        ip_type = settings.get("ip", "g")
        kapsig = settings.get("kapsig", 6.0)
        tell_bic = settings.get("tell-bic", 10.0)
        vcut = settings.get("vcut", 100.0)
        telluric_mode = settings.get("telluric", "add")
        molecules = settings.get("molec", "all")
        iphs = settings.get("iphs", 50)
        do_plot = settings.get("plot", 1)

        filenames = [frame.file for frame in frameset]
        n_files = len(filenames)

        print(f"Processing {n_files} input file(s)")
        print(f"Atmospheric models: {atm_data_dir}")
        print(
            f"deg_norm={deg_norm}, deg_wave={deg_wave}, ip={ip_type}, "
            f"kapsig={kapsig}, telluric={telluric_mode}"
        )

        if n_files == 0:
            raise ValueError("No input files provided")

        output_frameset = cpl.ui.FrameSet()

        for filename in filenames:
            print(f"\n--- {filename} ---")
            output_file = tellcorr(
                filename,
                atm_data_dir,
                output_dir=".",
                deg_norm=deg_norm,
                deg_wave=deg_wave,
                ip_type=ip_type,
                kapsig=kapsig,
                vcut=vcut,
                deg_bkg=deg_bkg,
                tell_bic=tell_bic,
                iphs=iphs,
                telluric=telluric_mode,
                molecules=molecules,
            )
            print(f"  -> {output_file}")
            if do_plot:
                plot_tellcorr(output_file)
            frame = cpl.ui.Frame(
                file=output_file,
                tag="TELLCORR_SPEC",
                group=cpl.ui.Frame.FrameGroup.PRODUCT,
            )
            output_frameset.append(frame)

        print(f"\nDone. Produced {len(filenames)} output file(s).")
        return output_frameset
