# Archive

The occultation modelling package from the first phase of this project. It is kept because the next step of the research needs it, not for its own sake.

## Why it is still here

The plan in the main README ends with characterising the HMI point spread function using lunar transits. Doing that means forward modelling the occultation: building the light curve, or the image, that the Moon crossing the solar disc should produce for a given limb-darkening profile, and comparing it against what HMI recorded. That forward model is `src/solarbatman/`, and it is tested.

The pieces that matter for the next step are `kernel.py`, which calls the batman occultation kernels with a separation supplied from an ephemeris rather than from a Keplerian orbit, and `geometry.py`, which computes that separation for a body seen from SDO.

## Reinstalling it

It was moved here, so an earlier editable install will no longer resolve.

```bash
pip install -e archive/
```

It needs `batman-package`, which the audit itself does not, so the reproduce workflow in `.github/workflows/` does not touch this directory. Run its tests yourself after installing it.

```bash
python -m pytest archive/tests -q
```

## What it established

The batman occultation kernel is accurate to better than seven parts in ten million against direct integration, across radius ratios from 0.0058 up to 1.03, which includes the total eclipse case. Exoplanet work only tests below about 0.2, so this was worth checking before relying on it for lunar geometry.

The exponential limb-darkening law has a pole at the limb, from a term in `1/(1 - exp(mu))`. For a large occulter this makes the light curve non-monotonic. Do not use it for solar work.

`examples/06`, `07` and `08` reproduce results already published by Howarth (2011) and Espinoza and Jordán (2015). They were written before that literature check. They are correct and they are not new, and they are kept as validation of the machinery rather than as findings.
