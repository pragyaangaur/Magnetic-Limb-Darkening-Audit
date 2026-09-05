# Data

Every number used anywhere in this project lives here, with its source written next to it. Nothing is invented and nothing is hard coded into an analysis script.

## `solar_clv_polynomials.csv`

The measured centre-to-limb variation of the Sun at 579.88 nm, from the two independent determinations that the field relies on. Pierce and Slaughter (1977) and Neckel and Labs (1994), transcribed from Table 1 of Hestroffer and Magnan (1998).

These two are the calibration anchor for the whole project. They agree very well over most of the disc and disagree near the limb, and that disagreement is the subject of the main result.

## `kostogryz2024_offsets.csv`

The limb-darkening offsets that Kostogryz et al. (2024) attribute to small-scale surface magnetic fields, taken from the Methods section of arXiv:2403.00118. These are the numbers the solar measurement has to be compared against.

## What is missing and should be added

The MPS-ATLAS and MURaM limb-darkening tables released with Kostogryz et al. (2026), arXiv:2606.21912. The paper says a public database of synthetic spectra at ten disc positions was released. Find it, download it, and put it here. Every result in `analysis/` currently uses only the two measured solar profiles, which is deliberate, but the comparison against their models needs their tables.

A modern solar centre-to-limb variation measurement. There is not one in this directory because as far as I can tell there is not a modern published one with a proper near-limb error budget. That absence is the point of the project.
