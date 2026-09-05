"""Shared routines: the model library, the measured offsets, and h1/h2.

Definitions follow Maxted (2018, 2023) and Kostogryz et al. (2024):

    h1 = I(mu = 2/3) / I(mu = 1)
    h2 = I(mu = 2/3) / I(mu = 1)  -  I(mu = 1/3) / I(mu = 1)

Both are dimensionless and of order 0.85 and 0.19 for a solar type star in the
Kepler band. Offsets between models and observations are quoted throughout in
absolute units of h1 and h2, never as relative percentages. Kostogryz et al.
(2024) write those absolute offsets multiplied by 100 and label them with a per
cent sign, which is easy to misread. See `data/maxted2023_table3.csv` for the
row by row check that fixes the convention.
"""

from __future__ import annotations

import csv
import gzip
import pathlib
from typing import ClassVar

import numpy as np
from scipy.interpolate import RegularGridInterpolator

DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
PASSBANDS = ("Kepler", "TESS", "CHEOPS", "PLATO")

# Sample means of the transit hosts whose limb darkening was measured.
# Maxted (2023), section 4.3.1.
SAMPLE = {"Teff": 6355.0, "logg": 4.39, "FeH": 0.23}

__all__ = [
    "SAMPLE",
    "LimbDarkeningLibrary",
    "clv_intensity",
    "h_from_power2",
    "load_facular_contrast",
    "load_sample_uncertainties",
    "load_solar_clv",
    "load_table3",
]


def h_from_power2(c, alpha):
    """
    Steepness parameters implied by a power-2 law ``I = 1 - c (1 - mu**alpha)``.

    Examples
    --------
    >>> h1, h2 = h_from_power2(0.732, 0.645)
    >>> float(round(h1, 4))
    0.8315
    """
    c = np.asarray(c, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    i23 = 1.0 - c * (1.0 - (2.0 / 3.0) ** alpha)
    i13 = 1.0 - c * (1.0 - (1.0 / 3.0) ** alpha)
    return i23, i23 - i13


class LimbDarkeningLibrary:
    """
    The MPS-ATLAS limb-darkening library of Kostogryz et al. (2022).

    Set 1 is the reference the 2024 magnetic paper calls REFLD. Set 2 is the
    same code with the other chemical abundance scale. The library is on a fine
    grid in metallicity, effective temperature and surface gravity, which is
    the point of it: the authors built it that way so that interpolation error
    would not dominate.

    Parameters
    ----------
    which : {"set1", "set2"}
        Which abundance set to load.
    passband : str
        One of `PASSBANDS`.
    """

    FILES: ClassVar[dict[str, str]] = {
        "set1": "kostogryz2022/set1_ldcoeffs.dat.gz",
        "set2": "kostogryz2022/set2_ldcoeffs.dat.gz",
    }
    CLV_FILES: ClassVar[dict[str, str]] = {
        "set1": "kostogryz2022/set1_clv.dat.gz",
        "set2": "kostogryz2022/set2_clv.dat.gz",
    }
    # The 24 disc positions the library tabulates, from the CDS byte-by-byte
    # description. mu = 1 is the disc centre column and is implicitly 1.0.
    MU_GRID: ClassVar[tuple[float, ...]] = (
        1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.22, 0.2,
        0.17, 0.15, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.03, 0.02, 0.01,
    )

    def __init__(self, which="set1", passband="Kepler", method="clv"):
        if method not in ("clv", "power2"):
            raise ValueError("method must be 'clv' or 'power2'")
        self.method = method
        if which not in self.FILES:
            raise ValueError(f"which must be one of {sorted(self.FILES)}")
        if passband not in PASSBANDS:
            raise ValueError(f"passband must be one of {PASSBANDS}")

        self.which = which
        self.passband = passband

        rows = []
        with gzip.open(DATA / self.FILES[which], "rt") as handle:
            for line in handle:
                if line[15:21].strip() != passband:
                    continue
                rows.append((
                    float(line[0:5]), float(line[6:10]), float(line[11:14]),
                    float(line[22:30]), float(line[31:39]),
                ))
        if not rows:
            raise ValueError(f"no rows found for passband {passband!r}")

        feh, teff, logg, c, alpha = (np.array(x) for x in zip(*rows, strict=True))
        self._axes = (np.unique(feh), np.unique(teff), np.unique(logg))
        shape = tuple(len(a) for a in self._axes)
        index = tuple(np.searchsorted(a, v) for a, v in
                      zip(self._axes, (feh, teff, logg), strict=True))

        grid_c = np.full(shape, np.nan)
        grid_a = np.full(shape, np.nan)
        grid_c[index] = c
        grid_a[index] = alpha
        if np.isnan(grid_c).any():
            raise ValueError("the grid has holes, which the readers assume it does not")

        self._c = RegularGridInterpolator(self._axes, grid_c,
                                          bounds_error=False, fill_value=None)
        self._a = RegularGridInterpolator(self._axes, grid_a,
                                          bounds_error=False, fill_value=None)
        self._clv = None
        if method == "clv":
            self._load_clv(passband)

    def _load_clv(self, passband):
        """
        Load the tabulated intensities at 24 disc positions.

        This is what Maxted (2023) uses for the Kostogryz comparison. His
        section 4.3.1 says so explicitly: for these models he takes the
        tabulated I(mu) directly with linear interpolation, rather than going
        through the fitted power-2 coefficients. Matching that matters, because
        the two routes differ by about 0.0005 in h1, which is a tenth of the
        offset under discussion.
        """
        rows = []
        with gzip.open(DATA / self.CLV_FILES[self.which], "rt") as handle:
            for line in handle:
                if line[15:21].strip() != passband:
                    continue
                values = line[35:].split()
                if len(values) < 23:
                    raise ValueError("unexpected number of intensity columns")
                rows.append((float(line[0:5]), float(line[6:10]), float(line[11:14]),
                             [1.0, *[float(v) for v in values[:23]]]))

        feh, teff, logg, clv = (list(x) for x in zip(*rows, strict=True))
        feh = np.array(feh); teff = np.array(teff); logg = np.array(logg)
        clv = np.array(clv)
        axes = (np.unique(feh), np.unique(teff), np.unique(logg))
        shape = (*[len(a) for a in axes], clv.shape[1])
        index = tuple(np.searchsorted(a, v) for a, v in
                      zip(axes, (feh, teff, logg), strict=True))
        grid = np.full(shape, np.nan)
        grid[index] = clv
        if np.isnan(grid).any():
            raise ValueError("the intensity grid has holes")
        self._clv = RegularGridInterpolator(axes, grid,
                                            bounds_error=False, fill_value=None)

    @property
    def axes(self):
        """The native grid, as (metallicity, temperature, gravity) arrays."""
        return self._axes

    def coefficients(self, feh, teff, logg):
        """Power-2 coefficients, linearly interpolated on the native grid."""
        point = np.atleast_2d([feh, teff, logg])
        return float(self._c(point)[0]), float(self._a(point)[0])

    def h(self, feh, teff, logg):
        """
        The two steepness parameters at the given stellar parameters.

        With ``method="clv"``, the default, the tabulated intensities are
        interpolated in the stellar parameters and then in mu, which is what
        Maxted (2023) does. With ``method="power2"`` the fitted power-2
        coefficients are used instead.

        Returns
        -------
        (float, float)
            ``h1`` and ``h2``.
        """
        if self.method == "power2":
            c, alpha = self.coefficients(feh, teff, logg)
            h1, h2 = h_from_power2(c, alpha)
            return float(h1), float(h2)

        profile = self._clv(np.atleast_2d([feh, teff, logg]))[0]
        mu = np.array(self.MU_GRID)
        order = np.argsort(mu)
        i23 = float(np.interp(2.0 / 3.0, mu[order], profile[order]))
        i13 = float(np.interp(1.0 / 3.0, mu[order], profile[order]))
        return i23, i23 - i13

    def derivative(self, feh, teff, logg, parameter, step):
        """
        Central difference of h1 and h2 with respect to one parameter.

        Parameters
        ----------
        parameter : {"FeH", "Teff", "logg"}
        step : float
            Half-width of the central difference, in the units of the
            parameter. Keep it at or above the native grid spacing, which is
            0.05 dex in metallicity, 100 K in temperature and 0.1 dex in
            gravity near solar values.

        Returns
        -------
        (float, float)
            ``dh1/dparameter`` and ``dh2/dparameter``, in absolute h units.
        """
        base = {"FeH": feh, "Teff": teff, "logg": logg}
        if parameter not in base:
            raise ValueError(f"parameter must be one of {sorted(base)}")

        low, high = dict(base), dict(base)
        low[parameter] -= step
        high[parameter] += step
        h1a, h2a = self.h(low["FeH"], low["Teff"], low["logg"])
        h1b, h2b = self.h(high["FeH"], high["Teff"], high["logg"])
        return (h1b - h1a) / (2 * step), (h2b - h2a) / (2 * step)


def _read_csv(name):
    with (DATA / name).open() as handle:
        return list(csv.DictReader(line for line in handle if not line.startswith("#")))


def load_table3():
    """
    Maxted (2023) Table 3, as published: observed minus model.

    A positive ``dh1`` means the star is brighter at mu = 2/3 than the model
    says, which is the direction Kostogryz et al. attribute to magnetic fields.
    The values in the CSV are taken from the journal version. The arXiv
    preprint, which is v1 only and predates acceptance, carries the opposite
    sign, so do not mix the two.
    """
    out = []
    for row in _read_csv("maxted2023_table3.csv"):
        out.append({
            "passband": row["passband"],
            "library": row["model_library"],
            "dh1": float(row["dh1"]),
            "dh1_err": float(row["dh1_err"]),
            "sigma_ext1": float(row["sigma_ext1"]),
            "dh2": float(row["dh2"]),
            "dh2_err": float(row["dh2_err"]),
            "sigma_ext2": float(row["sigma_ext2"]),
            "n_stars": int(row["n_stars"]),
            "note": row["note"],
            "is_refld": "Set 1" in row["model_library"],
        })
    return out


def load_sample_uncertainties():
    """Sample mean parameters and their typical per-star uncertainties."""
    return {row["quantity"]: {"value": float(row["value"]),
                              "uncertainty": float(row["typical_uncertainty"]),
                              "units": row["units"], "note": row["note"]}
            for row in _read_csv("maxted2023_sample.csv")}


def load_facular_contrast():
    """Measured contrast anchors for bright magnetic features."""
    return [{"feature": row["feature"],
             "mu_peak": float(row["mu_peak"]),
             "contrast_peak": float(row["contrast_peak"]),
             "contrast_centre": float(row["contrast_disc_centre"]),
             "source": row["source"], "note": row["note"]}
            for row in _read_csv("facular_contrast.csv")]


def load_solar_clv():
    """The two measured solar centre-to-limb polynomials at 579.88 nm."""
    return {row["source"]: {
        "wavelength_nm": float(row["wavelength_nm"]),
        "coefficients": np.array([float(row[f"a{i}"]) for i in range(6)]),
    } for row in _read_csv("solar_clv_polynomials.csv")}


def clv_intensity(mu, coefficients):
    """Intensity of a measured solar profile, normalised to disc centre."""
    mu = np.asarray(mu, dtype=float)
    return sum(c * mu**k for k, c in enumerate(coefficients)) / sum(coefficients)
