"""Telluric correction for CRIRES+ spectra.

Fits a forward model of telluric absorption to determine atmospheric
transmission, continuum shape, and wavelength calibration simultaneously.
Based on viper (Koehler & Zechmeister).
"""

import os
from collections import defaultdict

import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.special import erf

c = 299792.458  # speed of light [km/s]


# ---------------------------------------------------------------------------
# Parameter handling (from vipere)
# ---------------------------------------------------------------------------


class param(float):
    def __new__(cls, value, unc=None):
        instance = super().__new__(cls, value)
        instance.unc = unc
        instance.value = value
        return instance

    def __repr__(self):
        return f"{self.value}" + ("" if self.unc is None else f" +/- {self.unc}")


class nesteddict(dict):
    __getattr__ = dict.__getitem__

    def values(self):
        return [*super().values()]

    def keys(self):
        return [*super().keys()]

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self[key[0]][key[1]]
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            self[key[0]][key[1]] = value
        else:
            super().__setitem__(key, value)

    def __setattr__(self, key, value):
        self[key] = value

    def update(self, *args, **kwargs):
        for d in args + (kwargs,):
            for k in d:
                self[k] = d[k]

    def flat(self):
        d = {}
        for key, values in self.items():
            if isinstance(values, list):
                for i, val in enumerate(values):
                    d[(key, i)] = val
            elif isinstance(values, dict):
                for k, val in values.items():
                    d[(key, k)] = val
            else:
                d[key] = values
        return d

    def __add__(self, d):
        return self.__class__(self, d)

    def __repr__(self):
        return "\n".join(
            f"{k}: " + repr(v).replace("\n", ", ") for k, v in self.items()
        )


class Params(nesteddict):
    def __setitem__(self, key, value):
        super().__setitem__(key, self._as_param(value))

    def _as_param(self, value):
        if isinstance(value, (param, Params)):
            return value
        elif isinstance(value, (float, int)):
            return param(value)
        elif isinstance(value, tuple):
            return param(*value)
        elif type(value).__name__ in ("list", "ndarray"):
            return [self._as_param(val) for val in value]
        elif isinstance(value, dict):
            return Params(value)
        else:
            raise TypeError(f"{type(value)} not supported for param")

    def vary(self):
        return {k: v for k, v in self.flat().items() if v.unc != 0}


# ---------------------------------------------------------------------------
# Instrumental profiles (from vipere)
# ---------------------------------------------------------------------------


def IP_g(vk, s=2.2):
    """Gaussian IP."""
    ip_k = np.exp(-((vk / s) ** 2) / 2)
    ip_k /= ip_k.sum()
    return ip_k


def IP_sg(vk, s=2.2, e=2.0):
    """Super Gaussian."""
    ip_k = np.exp(-(abs(vk / s) ** e))
    ip_k /= ip_k.sum()
    return ip_k


def IP_ag(vk, s=2.2, a=0):
    """Asymmetric (skewed) Gaussian."""
    b = a / np.sqrt(1 + a**2) * np.sqrt(2 / np.pi)
    ss = s / np.sqrt(1 - b**2)
    vk = (vk + ss * b) / ss
    ip_k = np.exp(-(vk**2) / 2) * (1 + erf(a / np.sqrt(2) * vk))
    ip_k /= ip_k.sum()
    return ip_k


def IP_agr(vk, s, a=0):
    a = 10 * np.tanh(a / 10)
    return IP_ag(vk, s, a=a)


def IP_bg(vk, s1=2.0, s2=2.0):
    """BiGaussian."""
    xc = np.sqrt(2 / np.pi) * (-(s1**2) + s2**2) / (s1 + s2)
    vck = vk + xc
    ip_k = np.exp(-0.5 * (vck / np.where(vck < 0, s1, s2)) ** 2)
    ip_k /= ip_k.sum()
    return ip_k


def IP_mcg(vk, s0=2, a1=0.1):
    """Multiple central Gaussians."""
    s1 = 4 * s0
    a1 = a1 / 10
    ip_k = np.exp(-((vk / s0) ** 2))
    ip_k += a1 * np.exp(-((vk / s1) ** 2))
    ip_k = ip_k.clip(0, None)
    ip_k /= ip_k.sum()
    return ip_k


def IP_mg(vk, *a):
    """Multiple uniformly spaced Gaussians."""
    s = 0.9
    dx = s
    na = len(a) + 1
    mid = len(a) // 2
    a = np.tanh(a)
    a = [*a[:mid], 1, *a[mid:]]
    xl = np.arange(na)
    xm = np.dot(xl, a) / sum(a)
    xc = (dx * (xl - xm))[:, np.newaxis]
    ip_k = np.exp(-(((vk - xc) / s) ** 2))
    ip_k = np.dot(a, ip_k)
    ip_k /= ip_k.sum()
    return ip_k


IPs = {
    "g": IP_g,
    "sg": IP_sg,
    "ag": IP_ag,
    "agr": IP_agr,
    "bg": IP_bg,
    "mg": IP_mg,
    "mcg": IP_mcg,
}


# ---------------------------------------------------------------------------
# Forward model (from vipere)
# ---------------------------------------------------------------------------


def poly(x, a):
    return np.polyval(a[::-1], x)


class model:
    """Forward model: norm(x) * conv(IP, stellar * atm)(wavelength(x))."""

    def __init__(self, S_star, lnwave_j, fluxes_molec, IP, IP_hs=50, xcen=0):
        self.xcen = xcen
        self.S_star = S_star
        self.lnwave_j = lnwave_j
        self.fluxes_molec = fluxes_molec
        self.IP = IP
        self.dx = lnwave_j[1] - lnwave_j[0]
        self.IP_hs = IP_hs
        self.vk = np.arange(-IP_hs, IP_hs + 1) * self.dx * c
        self.lnwave_j_eff = lnwave_j[IP_hs:-IP_hs]

    def __call__(self, pixel, rv=0, norm=[1], wave=[], ip=[], atm=[], bkg=[0], **_kw):
        spec_gas = 1
        if len(self.fluxes_molec):
            flux_atm = np.nanprod(
                np.power(
                    self.fluxes_molec,
                    np.abs(atm[: len(self.fluxes_molec)])[:, np.newaxis],
                ),
                axis=0,
            )
            spec_gas = flux_atm

        Sj_eff = np.convolve(
            self.IP(self.vk, *ip),
            self.S_star(self.lnwave_j - rv / c) * (spec_gas + bkg[0]),
            mode="valid",
        )

        lnwave_obs = np.log(poly(pixel - self.xcen, wave))
        Si_eff = np.interp(lnwave_obs, self.lnwave_j_eff, Sj_eff)
        Si_mod = poly(pixel - self.xcen, norm) * Si_eff
        return Si_mod

    def fit(self, pixel, spec_obs, par, sig=None):
        varykeys, varyvals = zip(*par.vary().items())

        def S_model(x, *params):
            return self(x, **(par + dict(zip(varykeys, params))))

        sigma = sig if sig is not None and len(sig) else None
        params, e_params = curve_fit(
            S_model,
            pixel,
            spec_obs,
            p0=varyvals,
            sigma=sigma,
            absolute_sigma=False,
            epsfcn=1e-12,
        )
        pnew = par + dict(zip(varykeys, params))
        for k, v in zip(varykeys, np.sqrt(np.diag(e_params))):
            pnew[k].unc = v
        return pnew, e_params


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


class _nameddict(dict):
    __getattr__ = dict.__getitem__


flag = _nameddict(ok=0, nan=1, out=2, clip=4)


def local_sigma(resid, halfwin=50):
    """Sliding-window MAD-based sigma estimate (scaled to Gaussian std)."""
    n = len(resid)
    sigma = np.empty(n)
    for i in range(n):
        lo = max(0, i - halfwin)
        hi = min(n, i + halfwin + 1)
        window = resid[lo:hi]
        valid = window[np.isfinite(window)]
        if len(valid) > 3:
            sigma[i] = np.nanmedian(np.abs(valid - np.nanmedian(valid))) * 1.4826
        else:
            sigma[i] = np.inf
    return sigma


# ---------------------------------------------------------------------------
# Atmospheric model loading
# ---------------------------------------------------------------------------


def load_atm_models(atm_data_dir, wl_min_nm, wl_max_nm, molecules="all"):
    """Load atmospheric transmission models for the given wavelength range.

    Parameters
    ----------
    atm_data_dir : str
        Path to directory containing stdAtmos_*.fits files.
    wl_min_nm, wl_max_nm : float
        Wavelength range in nm.
    molecules : str or list
        'all' for automatic selection, or list of molecule names.

    Returns
    -------
    specs_molec_all : dict of lists
    wave_atm_all : dict of lists
    molec : list of str
    """
    wl_min_A = wl_min_nm * 10
    wl_max_A = wl_max_nm * 10

    bands_all = ["vis", "J", "H", "K", "L", "M"]
    wave_band = np.array([0, 9000, 14000, 18500, 28000, 40000])

    w0 = wl_min_A - wave_band
    w1 = wl_max_A - wave_band
    bands = bands_all[np.argmin(w0[w0 >= 0]) : int(np.argmin(w1[w1 >= 0]) + 1)]

    specs_molec_all = defaultdict(list)
    wave_atm_all = defaultdict(list)
    molec_names = None
    if isinstance(molecules, list) and molecules != ["all"]:
        molec_names = molecules

    for band in bands:
        filepath = os.path.join(atm_data_dir, f"stdAtmos_{band}.fits")
        with fits.open(filepath) as hdu:
            cols = hdu[1].columns.names
            data = hdu[1].data
            if molec_names is None:
                molec_names = [col for col in cols if col != "lambda"]
            for mol in molec_names:
                if mol in cols:
                    specs_molec_all[mol].extend(data[mol])
                    wave_atm_all[mol].extend(data["lambda"] * (1 + (-0.249 / 3e5)))

    molec = list(specs_molec_all.keys())
    return specs_molec_all, wave_atm_all, molec


# ---------------------------------------------------------------------------
# Segment fitting
# ---------------------------------------------------------------------------


def fit_segment(
    wl_nm,
    spec,
    err,
    specs_molec_all,
    wave_atm_all,
    molec,
    *,
    deg_norm=3,
    deg_wave=3,
    ip_type="g",
    kapsig=6.0,
    vcut=100,
    deg_bkg=0,
    tell_bic=10,
    iphs=50,
    telluric="add",
):
    """Fit telluric model to a single spectral segment.

    Parameters
    ----------
    wl_nm : ndarray
        Wavelength in nm.
    spec, err : ndarray
        Spectrum and error arrays.
    specs_molec_all, wave_atm_all, molec
        Atmospheric models from load_atm_models().
    deg_norm : int
        Polynomial degree for continuum normalization.
    deg_wave : int
        Polynomial degree for wavelength calibration.
    ip_type : str
        Instrumental profile model ('g', 'sg', 'ag', 'bg', ...).
    kapsig : float
        Kappa-sigma clipping threshold (0 to disable).
    vcut : float
        Velocity cutoff for edge trimming [km/s].
    deg_bkg : int
        Background polynomial degree (0 for none).
    tell_bic : float
        BIC improvement required to prefer telluric model.
    iphs : int
        IP half-size in pixels.
    telluric : str
        Telluric mode ('add' or 'add2').

    Returns
    -------
    dict or None
        Keys: 'wl_new', 'tell_model', 'cont_model', 'rms', 'params'.
        None if fitting failed.
    """
    wl = wl_nm * 10  # nm -> Angstroms
    pixel = np.arange(len(spec), dtype=float)

    flag_obs = np.zeros(len(spec), dtype=int)
    flag_obs[np.isnan(spec)] |= flag.nan

    # clip extreme outliers before any model is computed
    ok_mask = flag_obs == 0
    if np.sum(ok_mask) < 10:
        return None
    p17, median_flux, p83 = np.percentile(spec[ok_mask], [17, 50, 83])
    sig_crude = (p83 - p17) / 2
    flag_obs[spec > median_flux + 6 * sig_crude] |= flag.clip

    wl_min, wl_max = float(np.nanmin(wl)), float(np.nanmax(wl))

    # log-wavelength grid for this segment
    lnwave_j = np.arange(
        np.log(wl_min) + vcut / c,
        np.log(wl_max) - vcut / c,
        200 / 3e8,
    )
    if len(lnwave_j) < 2 * iphs + 2:
        return None

    # trim edge pixels
    flag_obs[np.log(wl) < lnwave_j[0]] |= flag.out
    flag_obs[np.log(wl) > lnwave_j[-1]] |= flag.out

    i_ok = np.where(flag_obs == 0)[0]
    if len(i_ok) < deg_wave + deg_norm + 5:
        return None
    pixel_ok = pixel[i_ok]
    spec_ok = spec[i_ok]

    xcen = np.nanmean(pixel_ok) + 18

    # set up atmospheric models on the log-wavelength grid
    specs_molec = np.zeros((0, len(lnwave_j)))
    par_atm = []
    for mol in molec:
        wave_mol = np.array(wave_atm_all[mol])
        spec_mol_raw = np.array(specs_molec_all[mol])
        s_mol = slice(*np.searchsorted(wave_mol, [wl_min, wl_max]))
        if len(wave_mol[s_mol]) > 0:
            spec_mol = np.interp(lnwave_j, np.log(wave_mol[s_mol]), spec_mol_raw[s_mol])
            specs_molec = np.r_[specs_molec, [spec_mol]]
            if np.nanstd(spec_mol) > 0.0001:
                par_atm.append((1, np.inf))
            else:
                par_atm.append((np.nan, 0))
        else:
            specs_molec = np.r_[specs_molec, [lnwave_j * 0 + 1]]
            par_atm.append((np.nan, 0))

    if telluric == "add2" and len(molec) > 1:
        par_atm_arr = np.asarray(par_atm)
        is_H2O = np.array(molec) == "H2O"
        if any(is_H2O):
            specs_molec = np.array(
                [
                    specs_molec[is_H2O][0],
                    np.nanprod(
                        specs_molec[~is_H2O]
                        * par_atm_arr[~is_H2O][:, 0].reshape(-1, 1),
                        axis=0,
                    ),
                ]
            )
            par_atm = [(1, np.inf), (1, np.inf)]
        else:
            specs_molec = np.nanprod(
                specs_molec * par_atm_arr[:, 0].reshape(-1, 1), axis=0
            )[np.newaxis]
            par_atm = [(1, np.inf)]

    def S_star(x):
        return 0 * x + 1

    IP_func = IPs[ip_type]

    S_mod = model(S_star, lnwave_j, specs_molec, IP_func, IP_hs=iphs, xcen=xcen)

    par = Params()
    par.rv = (0, 0)  # fixed, no stellar template
    par.norm = [np.nanmean(spec_ok)] + [0] * deg_norm
    par.wave = np.polyfit(pixel_ok - xcen, wl[i_ok], deg_wave)[::-1]
    parguess = Params(par)
    par.ip = [1.5]
    par.atm = par_atm
    if deg_bkg:
        par.bkg = [0]

    if ip_type in ("sg", "ag", "agr", "bg"):
        par.ip += (
            [2.0]
            if ip_type in ("sg",)
            else [1.0 if ip_type in ("ag", "agr") else par.ip[-1]]
        )

    sig = np.ones_like(spec)
    par3 = par

    # pre-fit kappa-sigma clip at 15-sigma
    if kapsig > 0:
        try:
            smod = S_mod(pixel, **par3)
            resid = spec - smod
            resid[flag_obs != 0] = np.nan
            flag_obs[np.abs(resid) >= 15 * local_sigma(resid)] |= flag.clip
            i_ok = np.where(flag_obs == 0)[0]
            pixel_ok = pixel[i_ok]
            spec_ok = spec[i_ok]
        except Exception:
            pass

    # pre-fit with Gaussian IP for complex IP types
    if ip_type in ("sg", "ag", "agr", "bg"):
        try:
            S_modg = model(S_star, lnwave_j, specs_molec, IP_g, IP_hs=iphs, xcen=xcen)
            par1 = Params(par, ip=par.ip[0:1])
            par2, _ = S_modg.fit(pixel_ok, spec_ok, par1, sig=sig[i_ok])
            par = par + par2.flat()
        except Exception:
            pass
    par3 = par

    # main fit
    par.wave = parguess.wave
    try:
        par4, e_params = S_mod.fit(pixel_ok, spec_ok, par, sig=sig[i_ok])
        par = par4
    except Exception as exc:
        print(f"    fit failed: {exc}")
        return None

    # post-fit kappa-sigma clip and refit
    if kapsig > 0:
        smod = S_mod(pixel, **par)
        resid = spec - smod
        resid[flag_obs != 0] = np.nan
        nr_k1 = np.count_nonzero(flag_obs)
        flag_obs[np.abs(resid) >= kapsig * local_sigma(resid)] |= flag.clip
        nr_k2 = np.count_nonzero(flag_obs)

        if nr_k1 != nr_k2:
            i_ok = np.where(flag_obs == 0)[0]
            pixel_ok = pixel[i_ok]
            spec_ok = spec[i_ok]
            try:
                par5, e_params = S_mod.fit(pixel_ok, spec_ok, par3, sig=sig[i_ok])
                par = par5
            except Exception:
                pass

    # BIC: compare telluric vs no-telluric model
    if len(specs_molec):
        fmod_tell = S_mod(pixel_ok, **par)
        rss_tell = np.sum((spec_ok - fmod_tell) ** 2)
        k_tell = len(par.vary())
        n_data = len(pixel_ok)
        bic_tell = n_data * np.log(rss_tell / n_data) + k_tell * np.log(n_data)

        S_mod_notell = model(S_star, lnwave_j, [], IP_func, IP_hs=iphs, xcen=xcen)
        par_notell = Params(parguess)
        if deg_bkg:
            par_notell.bkg = [0]
        par_notell.ip = [1.5]

        try:
            par_nt, _ = S_mod_notell.fit(pixel_ok, spec_ok, par_notell, sig=sig[i_ok])
            fmod_notell = S_mod_notell(pixel_ok, **par_nt)
            rss_notell = np.sum((spec_ok - fmod_notell) ** 2)
            k_notell = len(par_nt.vary())
            bic_notell = n_data * np.log(rss_notell / n_data) + k_notell * np.log(
                n_data
            )

            if bic_notell <= bic_tell + tell_bic:
                print(
                    "    BIC: no-telluric preferred (%.1f vs %.1f)"
                    % (bic_notell, bic_tell)
                )
                par_nt.atm = [(np.nan, 0)] * len(par_atm)
                par = par_nt
                S_mod = S_mod_notell
        except Exception:
            pass

    # compute outputs
    fmod = S_mod(pixel_ok, **par)
    res = spec_ok - fmod
    rms = np.nanstd(res) / np.nanmean(fmod) * 100

    wl_fit_A = poly(pixel - xcen, [float(w) for w in par.wave])
    wl_new = wl_fit_A / 10  # Angstroms -> nm

    cont_model = poly(pixel - xcen, [float(n) for n in par.norm])

    full_model = S_mod(pixel, **par)
    with np.errstate(divide="ignore", invalid="ignore"):
        tell_model = np.where(cont_model != 0, full_model / cont_model, 1.0)

    return {
        "wl_new": wl_new,
        "tell_model": tell_model,
        "cont_model": cont_model,
        "rms": rms,
        "params": par,
    }


# ---------------------------------------------------------------------------
# File processing
# ---------------------------------------------------------------------------


def _fit_dv_polynomial(results, pipeline_wl, prms_max, deg):
    """Fit a polynomial dv(wavelength) in km/s from orders with good telluric fits.

    The chips are at fixed relative positions, so the residual velocity offset
    between the vipere telluric solution and the pipeline calibration varies
    smoothly across all 3 chips. This function fits that offset so it can be
    extrapolated to orders where the telluric fit failed or was rejected.

    Returns (coeffs, n_good) where coeffs is in numpy.polyfit order (highest
    degree first), or (None, n_good) if not enough samples.
    """
    wl_samples = []
    dv_samples = []
    n_good = 0
    for key, res in results.items():
        if res["rms"] > prms_max:
            continue
        if not any(np.isfinite(p.value) for p in res["params"].atm):
            continue
        n_good += 1
        wl_pipe = pipeline_wl[key]
        wl_vip = res["wl_new"]
        n = len(wl_pipe)
        idx = np.linspace(int(0.05 * n), int(0.95 * n), 20).astype(int)
        ok = np.isfinite(wl_pipe[idx]) & np.isfinite(wl_vip[idx]) & (wl_pipe[idx] > 0)
        if not np.any(ok):
            continue
        wp = wl_pipe[idx[ok]]
        wv = wl_vip[idx[ok]]
        wl_samples.extend(wp)
        dv_samples.extend((wv - wp) / wp * c)
    if n_good < 3 or len(wl_samples) < deg + 2:
        return None, n_good
    coeffs = np.polyfit(wl_samples, dv_samples, deg)
    return coeffs, n_good


def tellcorr(
    filename,
    atm_data_dir,
    *,
    output_dir=".",
    deg_norm=3,
    deg_wave=3,
    ip_type="g",
    kapsig=6.0,
    vcut=100,
    deg_bkg=0,
    tell_bic=10,
    iphs=50,
    telluric="add",
    molecules="all",
    wl_interp_prms_max=30.0,
    wl_interp_deg=1,
    wl_interp_max_dv=100.0,
):
    """Run telluric correction on a CRIRES+ FITS file.

    Parameters
    ----------
    filename : str
        Input FITS file path.
    atm_data_dir : str
        Path to directory with stdAtmos_*.fits files.
    output_dir : str
        Output directory for the result file.
    wl_interp_prms_max : float
        Per-segment rms% threshold below which a telluric fit is considered
        good enough to refine the wavelength solution. Orders above this are
        considered unfitted and get a wavelength solution interpolated from
        the good orders.
    wl_interp_deg : int
        Polynomial degree for the dv(wavelength) interpolation across orders.
    wl_interp_max_dv : float
        Maximum allowed |dv| [km/s] of the interpolated correction; if
        exceeded the order keeps its pipeline wavelength solution.
    Other parameters are passed to fit_segment().

    Returns
    -------
    output_file : str
        Path to the output FITS file.
    """
    fit_kwargs = dict(
        deg_norm=deg_norm,
        deg_wave=deg_wave,
        ip_type=ip_type,
        kapsig=kapsig,
        vcut=vcut,
        deg_bkg=deg_bkg,
        tell_bic=tell_bic,
        iphs=iphs,
        telluric=telluric,
    )

    with fits.open(filename) as hdul:
        # discover order names from first detector
        chip1_cols = hdul["CHIP1.INT1"].data.dtype.names
        order_names = sorted(
            set(col.rsplit("_", 1)[0] for col in chip1_cols if col.endswith("_SPEC"))
        )
        n_pixels = len(hdul["CHIP1.INT1"].data)

    # determine global wavelength range across all detectors and orders
    wl_all = []
    with fits.open(filename) as hdul:
        for chip in ("CHIP1.INT1", "CHIP2.INT1", "CHIP3.INT1"):
            for order in order_names:
                wl = hdul[chip].data[f"{order}_WL"]
                if np.any(np.isfinite(wl)):
                    wl_all.extend([np.nanmin(wl), np.nanmax(wl)])
    wl_min_nm, wl_max_nm = np.nanmin(wl_all), np.nanmax(wl_all)

    mol_arg = (
        molecules.split(",")
        if isinstance(molecules, str) and molecules != "all"
        else molecules
    )
    specs_molec_all, wave_atm_all, molec = load_atm_models(
        atm_data_dir, wl_min_nm, wl_max_nm, mol_arg
    )
    print(f"Loaded atmospheric models: {molec}")

    chips = ("CHIP1.INT1", "CHIP2.INT1", "CHIP3.INT1")
    results = {}
    pipeline_wl = {}

    for chip in chips:
        with fits.open(filename) as hdul:
            data = hdul[chip].data
        for order in order_names:
            wl = data[f"{order}_WL"].astype(np.float64)
            spec = data[f"{order}_SPEC"].astype(np.float64)
            err = data[f"{order}_ERR"].astype(np.float64)
            pipeline_wl[(chip, order)] = wl

            print(f"  {chip} {order}: ", end="", flush=True)
            result = fit_segment(
                wl,
                spec,
                err,
                specs_molec_all,
                wave_atm_all,
                molec,
                **fit_kwargs,
            )
            if result is not None:
                print(f"rms={result['rms']:.2f}%")
                results[(chip, order)] = result
            else:
                print("skipped")

    # Fit a smooth dv(wavelength) across orders with good telluric fits, so
    # that orders without a usable fit can still be wavelength-corrected.
    dv_coeffs, n_good = _fit_dv_polynomial(
        results, pipeline_wl, wl_interp_prms_max, wl_interp_deg
    )
    if dv_coeffs is not None:
        print(
            f"dv(wavelength) interpolation: {n_good} good orders, deg={wl_interp_deg}"
        )
    else:
        print(
            f"dv(wavelength) interpolation: only {n_good} good orders "
            "(need >=3); unfitted orders keep pipeline WL"
        )

    # write output
    base = os.path.basename(filename)
    name, ext = os.path.splitext(base)
    output_file = os.path.join(output_dir, f"{name}_tellcorr{ext}")

    def _is_good(res):
        return res["rms"] <= wl_interp_prms_max and any(
            np.isfinite(p.value) for p in res["params"].atm
        )

    with fits.open(filename) as hdul:
        for chip in chips:
            orig_cols = hdul[chip].columns

            new_col_list = []
            for order in order_names:
                key = (chip, order)
                if key in results and _is_good(results[key]):
                    tell = results[key]["tell_model"]
                    cont = results[key]["cont_model"]
                else:
                    tell = np.ones(n_pixels)
                    cont = np.ones(n_pixels)
                new_col_list.append(
                    fits.Column(name=f"{order}_TELL", format="D", array=tell)
                )
                new_col_list.append(
                    fits.Column(name=f"{order}_CONT", format="D", array=cont)
                )

            all_cols = fits.ColDefs(orig_cols) + fits.ColDefs(new_col_list)
            new_hdu = fits.BinTableHDU.from_columns(all_cols, header=hdul[chip].header)
            # update WL: vipere fit for good orders, dv-interpolated pipeline
            # WL for the rest (when a dv polynomial is available).
            for order in order_names:
                key = (chip, order)
                if key in results and _is_good(results[key]):
                    new_hdu.data[f"{order}_WL"] = results[key]["wl_new"]
                elif dv_coeffs is not None:
                    wl_pipe = pipeline_wl[key]
                    dv = np.polyval(dv_coeffs, wl_pipe)
                    if np.nanmax(np.abs(dv)) > wl_interp_max_dv:
                        print(
                            f"  {chip} {order}: dv={np.nanmedian(dv):.1f} km/s "
                            "exceeds limit, keeping pipeline WL"
                        )
                        continue
                    new_hdu.data[f"{order}_WL"] = wl_pipe * (1 + dv / c)
                    print(
                        f"  {chip} {order}: dv-interpolated "
                        f"({np.nanmedian(dv):+.2f} km/s)"
                    )

            hdul[chip] = new_hdu

        hdul[0].header["HISTORY"] = (
            "Telluric correction applied by cr2res_util_tellcorr"
        )
        hdul.writeto(output_file, overwrite=True)

    return output_file


def plot_tellcorr(filename, output_png=None):
    """Plot observed spectra with telluric model overlay.

    Parameters
    ----------
    filename : str
        Tellcorr output FITS file.
    output_png : str, optional
        Output PNG path. Defaults to filename with .png extension.
    """
    import matplotlib.pyplot as plt

    if output_png is None:
        output_png = os.path.splitext(filename)[0] + ".png"

    with fits.open(filename) as hdul:
        chips = ["CHIP1.INT1", "CHIP2.INT1", "CHIP3.INT1"]
        cols0 = hdul[chips[0]].columns.names
        orders = sorted(set(c.rsplit("_", 1)[0] for c in cols0 if c.endswith("_SPEC")))

        fig, (ax1, ax2) = plt.subplots(
            2,
            1,
            figsize=(16, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [3, 1]},
        )

        for chip in chips:
            data = hdul[chip].data
            for order in orders:
                wl = data[f"{order}_WL"]
                spec = data[f"{order}_SPEC"]
                tell = data[f"{order}_TELL"]
                cont = data[f"{order}_CONT"]
                if np.all(np.isnan(spec)):
                    continue
                ax1.plot(wl, spec, color="0.5", lw=0.4)
                ax1.plot(wl, cont * tell, color="C3", lw=0.6)
                ax2.plot(wl, tell, color="C0", lw=0.5)

        ax1.plot([], [], color="0.5", lw=0.8, label="observed")
        ax1.plot([], [], color="C3", lw=0.8, label="model")
        ax1.legend(fontsize=9)
        ax1.set_ylabel("Flux")
        ax2.set_ylabel("Telluric transmission")
        ax2.set_xlabel("Wavelength [nm]")
        ax2.set_ylim(-0.05, 1.1)
        fig.suptitle(os.path.basename(filename), fontsize=11)
        fig.tight_layout()
        fig.savefig(output_png, dpi=150)
        plt.close(fig)

    print(f"Saved {output_png}")
