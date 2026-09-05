# A systematics audit of the magnetic limb-darkening result

Kostogryz et al. (2024, Nature Astronomy 8, 929) found that limb darkening from non-magnetic stellar atmosphere models is too steep compared with transit observations, and showed that small-scale surface magnetic fields account for the difference. They present the simultaneous fit to both steepness parameters, h1 and h2, as the strength of the result.

This repository audits that claim. The h1 half survives every test applied here. The h2 half is the size of an analysis systematic that Maxted (2023), the paper the observed values come from, measures and warns about explicitly, and that the 2024 paper does not mention.

This is independent work and it has not been peer reviewed. Everything in it is computed from numbers published by other people, and `analysis/verify.py` checks the code against fourteen of those numbers before any result is produced. All fourteen pass, including an exact reproduction of a metallicity correction Maxted quotes. If something here is wrong, the fastest way to show it is to run the scripts and change the inputs, which is why the inputs are all in `data/` with their sources written next to them.

![Systematics budget](results/figures/systematics_budget.png)

## Start here

Read [`notes/00-result.md`](notes/00-result.md) for the full argument with every number. Read [`notes/01-explain-it-simply.md`](notes/01-explain-it-simply.md) for the same argument with no equations. Read [`notes/02-literature-audit.md`](notes/02-literature-audit.md) to see what is already published, and [`notes/03-sources.md`](notes/03-sources.md) for the papers this rests on and how much of each was read.

## Results

**h2 does not carry the weight placed on it.** Maxted (2023) section 4.2 reanalysed the sixteen stars common to his 2018 and 2023 studies, changing only the data processing and the fitting code. h1 shifted by 0.000 ± 0.008 and h2 shifted by 0.010 ± 0.002, and he concludes that h2 "may be affected by systematic errors ~0.01 depending on the details of the analysis". The h2 offset being interpreted is 0.012. Separately, his section 4.3.4 quotes the MURaM 100 G prediction from Norris et al. (2017) as +0.007 in h1 and −0.005 in h2, against observed values of +0.006 and −0.012, so h1 agrees to fifteen per cent while h2 is out by a factor of 2.4. See `analysis/a06_h2_reliability.py`.

**Stellar parameter systematics cannot do it.** Producing the offset needs a metallicity scale error of −0.44 dex or a temperature scale error of +247 K, against SWEET-Cat uncertainties of 0.05 dex and 60 K. The values required by the Kepler and TESS samples agree with each other, which is what a genuine scale error looks like. Ten times the formal error in the same direction for two independently assembled samples is not credible. See `analysis/a02_stellar_parameter_systematics.py`.

**Bright magnetic features on the transit chord cannot do it.** This is the objection Maxted (2023) raised in reply to his referee, and it has not been quantified until now. The offset requires the star to be 0.71 per cent brighter than the model at μ = 2/3 and 2.74 per cent brighter at μ = 1/3. With the facular contrast Pietrow et al. (2026) measured using HMI, that needs 69 per cent of the chord covered in faculae, and with the largest contrast in the literature it is still 27 per cent. Neither is possible for the bright, quiet stars this sample is made of. See `analysis/a03_facular_chord_test.py`.

Both of those results support the magnetic interpretation, and neither appears to have been shown before.

**The reference model is the largest term.** The offset is a residual against one library. Maxted (2023) Table 3 gives the same measurement against nine non-magnetic libraries in the Kepler band, and Δh1 spans 0.057 across them, which is 9.5 times the signal. One library gives the opposite sign. Nothing in the transit data chooses between them. See `analysis/a01_signal_versus_library_spread.py`.

| | Δh1 |
|---|---|
| offset attributed to magnetic fields | 0.006 |
| abundance-scale choice within MPS-ATLAS alone | 0.004 |
| spread across nine published non-magnetic libraries | 0.057 |
| spread between the two measured solar profiles | 0.0011 |

**The Sun could rank the libraries, and has been used on only two of them.** Norris et al. (2017) Figure 4 already compares the measured solar profile against MURaM from field-free to 100 G, and Kostogryz et al. (2024) Figure 2 repeats it for their SSD case, so the solar check on MURaM exists and comes out in their favour. The other eight libraries in Maxted's Table 3 have not been checked against the Sun, and it is the choice among them that sets the size of the inferred field. The two published solar profiles agree to 0.0011 in h1, a fiftieth of the library spread, so the measurement is sharp enough to settle the ranking. See `analysis/a04_solar_anchor.py`.

**What stands between the Sun and that ranking is the passband conversion, and it is worth six times the signal.** The measured solar profiles are tabulated by wavelength and the libraries by passband. In MPS-ATLAS at solar parameters, model h1 moves by 0.0384 across the Kepler, TESS, CHEOPS and PLATO bands, which is 6.4 times the offset under discussion. Carrying the measurement into the band is therefore the whole measurement, and it is the object Kostogryz et al. (2022) section 5.2 built for their own library under the name NL_Kepler and applied to no other. No solar residual is quoted anywhere in this repository, because a 579.88 nm measurement against a broad-band model has a difference of unknown sign and size. See part one of `analysis/a07_solar_residual.py`.

**A solar correction is applied to a sample that is not solar.** Kostogryz et al. run MURaM for a solar atmosphere only, on the stated grounds that the sample has near-solar parameters. The sample mean is 583 K hotter and 0.23 dex more metal rich than the Sun. Inside the reference library, with the same passband on both sides so that no wavelength conversion enters, that parameter gap moves non-magnetic limb darkening by +0.0159 in h1 and −0.0153 in h2 in the Kepler band, which is 2.6 and 1.3 times the offsets being explained. The Kepler and TESS values agree, so this is a property of the atmospheres rather than an artefact of one passband. This does not explain the offset away, because Maxted computes the model prediction at each star's own parameters and the residual is already differential in that respect. What it says is that the magnetic correction is computed in a solar box and carried unchanged across a gap over which the non-magnetic structure moves by several times the effect, with no bound given on the error in that step. See part two of `analysis/a07_solar_residual.py`.

## Three things to know before reading anything

**The offsets are absolute, not relative.** Kostogryz et al. write "Δh1' = 0.6 % ± 0.2 %". Maxted's abstract states the same quantity as "a small but significant offset Δh1' ~ 0.006". They are absolute differences multiplied by 100, and all four Kepler and TESS values match Table 3 exactly with the sign flipped.

**The magnetic explanation is Maxted's, and it is not new in 2024.** His 2023 abstract ascribes the offset to "a mean vertical magnetic field strength ~100 G that is expected in the photospheres of these inactive solar-type stars", and he suggested it in 2018. Norris et al. (2017) had already found the measured solar gradient shallower than field-free MURaM and attributed it to residual quiet-Sun field, with Krivova and Solanki among the authors. The contribution of Kostogryz et al. (2024) is the full MURaM treatment across wavelength and magnetisation, which is a different and real contribution, and they cite Maxted properly.

**The sample is not solar.** Maxted's means are Teff = 6355 K, log g = 4.39 and [Fe/H] = +0.23, for 43 FGK stars dominated by F types. Kostogryz et al. restrict their MURaM simulations to a solar atmosphere on the stated grounds that the sample has near-solar parameters. It is 583 K hotter than the Sun and metal rich, and their own 2026 follow-up finds the magnetic effect grows towards hotter and more metal-rich stars.

## Reproducing it

```bash
pip install numpy scipy matplotlib
```

Run the verification first, because it checks the table readers, the interpolation and the unit convention, and nothing downstream is worth reading until it passes.

```bash
python analysis/verify.py
```

Each analysis script is independent of the others. Run them from anywhere.

```bash
python analysis/a01_signal_versus_library_spread.py
```

| script | what it establishes |
|---|---|
| `verify.py` | fourteen cross-checks against numbers printed in the papers |
| `a01_signal_versus_library_spread.py` | the signal against the disagreement between nine libraries |
| `a02_stellar_parameter_systematics.py` | rules out a parameter scale error |
| `a03_facular_chord_test.py` | rules out facular coverage of the chord |
| `a04_solar_anchor.py` | the Sun is precise enough to rank the libraries |
| `a05_figure.py` | the figure above |
| `a06_h2_reliability.py` | whether the h2 offset is larger than its own systematic |
| `a07_solar_residual.py` | what the passband conversion and the near-solar assumption are worth |

`analysis/ldlib.py` holds the shared readers, the interpolation over the model grid and the h1 and h2 definitions. Every script writes into `results/`, and `results/README.md` maps each output file back to the script that produces it. The workflow in `.github/workflows/` runs all of this on every push and fails if the committed tables no longer match the ones the code builds.

## Data

Everything in `data/` carries its source in the file header, and nothing is invented or hard coded into an analysis script.

| file | what it is |
|---|---|
| `kostogryz2022/` | both abundance sets of the MPS-ATLAS limb-darkening coefficients from CDS (J/A+A/666/A60), the library the 2024 comparison is measured against |
| `maxted2023_table3.csv` | the measured offsets against nine libraries, transcribed from the published table |
| `maxted2023_sample.csv` | the sample mean parameters and their uncertainties |
| `maxted2023_h2_systematic.csv` | the section 4.2 reanalysis of sixteen common stars, which is the h2 analysis systematic |
| `kostogryz2024_offsets.csv` | the offsets as reported in the 2024 Methods section |
| `facular_contrast.csv` | the measured contrast of bright magnetic features |
| `solar_clv_polynomials.csv` | the two measured solar centre-to-limb profiles |

The papers themselves are not in this repository, because they are published copies that are not mine to redistribute. [`notes/03-sources.md`](notes/03-sources.md) lists every one of them with an arXiv number or a journal reference, together with the reading notes taken from it.

## Limitations

The solar numbers are at 579.88 nm and the comparison still has to be redone in the Kepler and TESS passbands, weighted by the response functions and the disc-centre continuum. `a07` measures how much that conversion is worth and the answer is 6.4 times the signal, so it is the next piece of work rather than a caveat on the end of the existing one. Doing it needs the full wavelength tables of Neckel and Labs (1994) and Pierce and Slaughter (1977), which are not in `data/`.

Maxted's Table 3 has been verified against the journal version at academic.oup.com and all fourteen rows match. The arXiv preprint is v1 only and carries the opposite sign, so the two versions should not be mixed.

The derivatives in `a02` are central differences at the sample mean. The proper calculation is star by star with each object's own parameters, which needs Maxted's Table 1 extracted from the PDF.

The 10 per cent facular contrast attributed to Yeo et al. (2013) lives in their Figure 14 and is Pietrow et al.'s reading of it. `a03` does not depend on it, because it scans contrast continuously and the required filling factor stays above 0.09 even at 30 per cent.

The MURaM 100 G values of +0.007 and −0.005 used in `a06` are the width of the shaded band in Figure 4 of Norris et al. (2017), read off a plot by Maxted. They are not tabulated anywhere and they carry an unstated reading error.

How much of each paper was read is recorded at the end of [`notes/03-sources.md`](notes/03-sources.md), so that no claim here rests on more reading than actually happened.

## Layout

    data/        every input number, with its source in the file header
    analysis/    seven analyses, a verification suite and a shared module
    results/     figures and tables, all regenerated by the scripts
    notes/       the result, the plain version, the literature audit and the sources
    archive/     earlier occultation modelling, kept for the lunar transit forward
                 model, which the solar calibration step would need

## Licence

BSD 3-Clause, in [`LICENSE`](LICENSE). Citation metadata is in [`CITATION.cff`](CITATION.cff).
