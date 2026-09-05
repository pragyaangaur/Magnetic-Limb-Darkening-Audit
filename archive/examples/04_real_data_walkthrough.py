"""Run the whole pipeline on real SDO data.

This is the script the rest of the package exists for. It downloads a sequence
of full disk images, builds a disk integrated light curve, measures limb
darkening from a quiet frame, gets the geometry from JPL Horizons, and fits the
occultation with the size held fixed by the ephemeris.

It needs sunpy, astropy, astroquery and matplotlib, and it needs a network
connection. Install those with:

    pip install -e ".[solar,plot]" astroquery

Two things to know before running it.

The first is that this script has never been executed end to end, because sunpy
does not import in the environment it was written in. Treat it as a starting
point rather than a finished tool, and check each step's output.

The second is that step four decides everything. Run it with ``--quiet-day``
first, on a date with no occultation, and look at the scatter in the light
curve. That number is the photometric noise floor, and it says which events are
within reach. A lunar transit removes most of the light and will work through
almost anything. A Venus transit is about 1180 parts per million deep. A
Mercury transit is 41 parts per million and may well be out of reach.

Examples:

    python examples/04_real_data_walkthrough.py --quiet-day 2019-11-08
    python examples/04_real_data_walkthrough.py --event 2019-11-11 --body Mercury
"""

import argparse
import pathlib

import numpy as np

OUTPUT = pathlib.Path(__file__).resolve().parents[1] / "figures"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--quiet-day",
        help="date with no occultation, used to measure the noise floor",
    )
    group.add_argument("--event", help="date of a transit or a lunar crossing")
    parser.add_argument(
        "--body", default="Moon", help="occulting body, as JPL Horizons names it"
    )
    parser.add_argument(
        "--hours", type=float, default=4.0, help="length of the window in hours"
    )
    parser.add_argument(
        "--cadence", type=int, default=300, help="seconds between samples"
    )
    parser.add_argument(
        "--wavelength",
        type=int,
        default=1600,
        help="AIA channel in angstrom. Use 1600 or 1700 for photospheric work. "
        "The extreme ultraviolet channels see optically thin coronal emission "
        "and do not follow a photospheric limb darkening law.",
    )
    return parser.parse_args()


def fetch_maps(start, hours, cadence, wavelength):
    """Download a sequence of AIA full disk images."""
    import astropy.units as u
    import sunpy.map
    from astropy.time import Time
    from sunpy.net import Fido
    from sunpy.net import attrs as a

    begin = Time(start)
    end = begin + hours * u.hour

    query = Fido.search(
        a.Time(begin, end),
        a.Instrument.aia,
        a.Wavelength(wavelength * u.angstrom),
        a.Sample(cadence * u.second),
    )
    print(query)

    files = Fido.fetch(query, path=str(OUTPUT.parent / "data" / "{file}"))
    if not files:
        raise SystemExit("nothing was downloaded, so there is nothing to do")

    print(f"Downloaded {len(files)} files")
    return sunpy.map.Map(sorted(files))


def measure_noise_floor(curve):
    """Scatter left after removing a slow trend, which is the number that matters."""
    flux = np.asarray(curve["flux"], dtype=float)
    index = np.arange(flux.size, dtype=float)

    # A quadratic is enough to take out thermal drift over a few hours without
    # eating a transit that lasts a similar time. Check the fit by eye before
    # trusting the residual.
    trend = np.polynomial.Polynomial.fit(index, flux, deg=2)
    residual = flux - trend(index)

    rms = float(np.std(residual))
    print()
    print(f"Samples                     {flux.size}")
    print(f"Peak to peak, raw           {np.ptp(flux):.3e}")
    print(f"Scatter after a quadratic   {rms:.3e}  ({rms * 1e6:.0f} ppm)")
    print()
    print("What this means for each event:")
    for label, depth in [
        ("lunar transit", 0.8),
        ("Venus transit", 1.18e-3),
        ("Mercury transit", 4.1e-5),
    ]:
        ratio = depth / rms if rms > 0 else np.inf
        verdict = "easy" if ratio > 100 else "reachable" if ratio > 10 else "hard"
        print(f"  {label:<18} depth / scatter = {ratio:8.1f}   {verdict}")
    return rms


def run_quiet_day(args):
    from solarbatman.photometry import light_curve

    print(f"Quiet day noise floor test, {args.quiet_day}, AIA {args.wavelength}")
    maps = fetch_maps(args.quiet_day, args.hours, args.cadence, args.wavelength)

    # reuse_mask is safe here because the pointing does not move over a few
    # hours. Do not use it across months.
    curve = light_curve(maps, reuse_mask=True)
    rms = measure_noise_floor(curve)

    save_curve_plot(curve, f"noise-floor-{args.quiet_day}", rms=rms)


def run_event(args):
    from solarbatman.fitting import fit_occultation
    from solarbatman.geometry import (
        radius_ratio_from_ephemeris,
        separation_from_ephemeris,
    )
    from solarbatman.limbdark import fit_intensity_profile
    from solarbatman.photometry import light_curve, radial_profile

    print(f"Occultation fit, {args.event}, body {args.body}, AIA {args.wavelength}")
    maps = fetch_maps(args.event, args.hours, args.cadence, args.wavelength)
    maps = list(maps)

    print("\nStep 1. Disk integrated light curve.")
    curve = light_curve(maps, reuse_mask=True)
    print(f"  {len(curve)} samples, scatter {np.std(curve['flux']):.3e}")

    print("\nStep 2. Limb darkening from the first frame.")
    mu, profile, counts = radial_profile(maps[0], mu_min=0.05)
    measured = fit_intensity_profile(mu, profile, law="quadratic")
    print(f"  u1 = {measured['coeffs'][0]:.4f}, u2 = {measured['coeffs'][1]:.4f}")
    print(f"  residual rms {measured['residual_rms']:.3e} over {measured['n_points']} bins")
    print(
        "  If the residual is much larger than the pixel noise, the quadratic "
        "law is the wrong shape and another law will fit better."
    )

    print("\nStep 3. Geometry from JPL Horizons. Nothing is fitted here.")
    times = curve["time"]
    z = separation_from_ephemeris(times, args.body, observer="SDO")
    ratio = radius_ratio_from_ephemeris(times, args.body, observer="SDO")
    ratio_mid = float(np.median(ratio))
    print(f"  closest approach   z = {np.min(z):.4f} solar radii")
    print(f"  radius ratio       {ratio_mid:.5f}")
    if np.min(z) > 1 + ratio_mid:
        raise SystemExit(
            "the body never reaches the solar disk in this window, so there is "
            "no occultation to fit. Check the date and the observer."
        )

    print("\nStep 4. Fit, with the size fixed by the ephemeris.")
    result = fit_occultation(
        z,
        np.asarray(curve["flux"], dtype=float),
        ratio_guess=ratio_mid,
        fit_ratio=False,
        coeffs_guess=measured["coeffs"],
    )
    print(f"  u1 = {result['coeffs'][0]:.4f}  (image gave {measured['coeffs'][0]:.4f})")
    print(f"  u2 = {result['coeffs'][1]:.4f}  (image gave {measured['coeffs'][1]:.4f})")
    print(f"  baseline    {result['baseline']:.6f}")
    print(f"  residual    {result['residual_rms']:.3e}")
    print()
    print(
        "The two limb darkening measurements are independent. One comes from "
        "the shape of the disk in an image, the other from the shape of a light "
        "curve. If they agree, the whole chain is consistent."
    )

    save_curve_plot(curve, f"{args.body.lower()}-{args.event}", model=result["model"])


def save_curve_plot(curve, name, model=None, rms=None):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib is not installed, so no figure was written.")
        return

    OUTPUT.mkdir(parents=True, exist_ok=True)
    flux = np.asarray(curve["flux"], dtype=float)
    index = np.arange(flux.size)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(index, flux, ".", ms=3, color="0.4", label="disk integrated flux")
    if model is not None:
        ax.plot(index, model, lw=1.5, color="C3", label="fitted model")
    ax.set_xlabel("sample")
    ax.set_ylabel("relative flux")
    title = name
    if rms is not None:
        title = f"{name}, scatter {rms * 1e6:.0f} ppm"
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()

    path = OUTPUT / f"{name}.png"
    fig.savefig(path, dpi=130)
    print(f"\nFigure written to {path}")


def main():
    args = parse_args()
    if args.quiet_day:
        run_quiet_day(args)
    else:
        run_event(args)


if __name__ == "__main__":
    main()
