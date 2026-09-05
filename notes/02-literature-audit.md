# Literature audit

What is already published, so that nobody in this project spends a week rediscovering it. Written 3 September 2026. I read the first three in full and only the abstracts of the rest, which is marked below.

## Already published, do not claim these

**Transit-fitted limb-darkening coefficients are not comparable with profile-fitted ones.** Howarth (2011), MNRAS 418, 1165. Read in part. He shows the two procedures minimise different things, that the difference grows with impact parameter, and introduces SPAM as the fix: generate a synthetic light curve from the model profile and fit it the same way you fit data, so you compare like with like. Espinoza and Jordán extended this to MC-SPAM.

What is left open is the size of the effect for the real Sun, measured rather than modelled, in the specific h1′ and h2′ parametrisation that Kostogryz et al. use. Nothing in this repository does that, because it needs the solar profile carried into a passband first, which is the gap `analysis/a07` measures.

**Radius bias from limb-darkening errors, and its dependence on impact parameter.** Espinoza and Jordán (2015), MNRAS 450, 1879. Read in full. Section 4.3 and Figures 10 to 14. Fixing the coefficients to wrong values biases the radius ratio by up to about 3 percent and floating them by about 1 percent. They also report biases up to 10 percent in the scaled semi-major axis and 5 percent in inclination, and they state that the impact parameter itself varies significantly with the limb-darkening assumption.

Anything computed about limb-darkening-induced transit parameter bias should be checked here first. `archive/examples/06`, `07` and `08` reproduce parts of this. They were written before I read the paper.

**Which limb-darkening law to prefer.** Espinoza and Jordán (2016), MNRAS 457, 3573. Abstract only. Three-parameter, square-root and logarithmic laws beat linear and quadratic.

**Model limb darkening is too steep, and magnetic fields explain it.** Kostogryz et al. (2024), Nature Astronomy 8, 929. Read in full. The offsets are Δh1′ = 0.6 ± 0.2 percent and Δh2′ = −1.2 ± 0.4 percent for Kepler, and 0.4 ± 0.3 and −0.9 ± 0.4 for TESS. Points at μ = 2/3 are observed 0.5 percent brighter than non-magnetic models predict, and points at μ = 1/3 are 1.5 percent brighter. MURaM simulations with 100, 200 and 300 G added fields reproduce this.

Note what it says at the end of the abstract: limb darkening could be used to measure the small-scale magnetic field on stars with transiting planets. That is the forward-looking claim, and nobody has checked it on the Sun.

**The magnetic effect across the main sequence.** Kostogryz et al. (2026), arXiv:2606.21912, submitted 20 June 2026. Abstract only. Covers 3200 to 6800 K and a range of metallicity. Strongest for hotter and more metal-rich stars, negligible for M dwarfs. Releases a public database of synthetic spectra at ten disc positions. Does not convert the effect into transit parameter biases and does not test against the Sun.

**Transiting synthetic planets across real HMI images.** Morris, Bobra, Agol, Lee, Hawley (2020), MNRAS 493, 5489. Abstract and summary only. Earth-sized planets across HMI continuum images, solar variability modelled with a Gaussian process, granulation and oscillation amplitude around 100 ppm, and a resulting 3.6 percent radius uncertainty driven by a degeneracy between radius ratio, limb darkening and impact parameter.

Do not propose putting a synthetic planet across an HMI image. It has been done.

**Solar cycle variation of a limb-darkening ratio.** Kim et al. (2017), JASS 34, 99. Summary only. Centre to limb intensity ratio at 0.9 solar diameter from SOHO/MDI (1996 to 2011) and SDO/HMI (2010 to 2016). MDI values run 0.65 to 0.70 and correlate with activity. HMI values sit near 0.65 and decline after 2014, which they attribute to CCD ageing.

This matters a lot. Somebody has already looked for a solar cycle signal in SDO limb darkening and found a trend they could not separate from instrumental degradation. That is precisely the problem the lunar transit calibration is meant to solve, and it is evidence that the problem is real rather than hypothetical.

**Lunar transits as an HMI stray light and point spread function diagnostic.** arXiv:2511.13348. Summary only. Builds an artificial lunar transit image from the ephemeris and a limb-darkening profile, convolves with a candidate point spread function, and compares against the observed transit.

Read this carefully before starting step two, and watch for the circularity: the method assumes a limb-darkening profile in order to get the point spread function, and we want the point spread function in order to measure the limb-darkening profile. Solving for both together is part of the work.

**Solar centre-to-limb variation measurements.** Pierce and Slaughter (1977), Solar Physics 52, 179, and Neckel and Labs (1994), Solar Physics 153, 91. Coefficients at 579.88 nm taken from Table 1 of Hestroffer and Magnan (1998), A&A 333, 338, which I read in full. Hestroffer and Magnan make the same point this project keeps hitting: the two agree excellently for r/R below 0.9 and the polynomial coefficients themselves diverge, because a polynomial cannot represent the profile at the limb where the radial derivative goes to infinity.

## Open, as far as I can tell

**How the magnetic signal compares with the disagreement between the libraries it is measured against.** Maxted (2023) published the numbers that make this computable and did not frame it this way. `analysis/a01`. This is the main result.

**Whether a stellar parameter scale error could produce the offset.** Ruled out, `analysis/a02`. It would need half a dex in metallicity or 240 K in temperature.

**Whether bright magnetic features on the transit chord could produce it.** Ruled out, `analysis/a03`. Maxted raises the objection qualitatively in his section 4.3.2 and nobody attached a number to it. It needs a filling factor of at least a quarter of the chord.

**Ranking the libraries against the Sun.** Not done. The solar profiles are tabulated by wavelength, the libraries by passband, and only Kostogryz et al. (2022) built the conversion, for their own library. `analysis/a04` shows the solar measurement is a fiftieth of the library spread, so it would settle the ranking easily.

**Whether the h2 offset is larger than its own analysis systematic.** Maxted (2023) section 4.2 measures that systematic and warns about it in the paper the offsets come from, and Kostogryz et al. (2024) do not mention it. `analysis/a06`. This is the finding most likely to matter to the authors.

**What the near-solar assumption is worth.** The magnetic correction is computed in a solar MURaM box and applied to a sample 583 K hotter and 0.23 dex more metal rich. `analysis/a07` measures the gap inside the reference library at 2.6 times the h1 offset, in one passband, so no wavelength conversion enters.

**Whether solar limb darkening varies over the cycle in the way the magnetic mechanism predicts.** Kim et al. (2017) is the nearest attempt and could not separate a real signal from CCD ageing. Still open, but it is a longer project than the ranking and it should come second.

## Read properly in the second pass

**Maxted (2023), MNRAS 519, 3723, arXiv:2212.09117.** The source of the offsets. Table 3, sections 4.3.1 and 4.3.2. His sign convention is model minus observed. His sample means are Teff 6355 K, log g 4.39, [Fe/H] +0.23, taken from SWEET-Cat, so the sample is F-star dominated and metal rich rather than solar. Section 4.3.2 contains the exchange with his referee about whether the models were tuned to solar measurements, which is where the facular objection comes from.

**Kostogryz et al. (2022), A&A 666, A60, arXiv:2206.06641.** The library REFLD comes from. Section 5.2 describes NL_Kepler, the Neckel and Labs profile carried into the Kepler band, which is the object needed to rank the libraries. Used only to validate their own models.

**Pietrow et al. (2026), arXiv:2606.28887.** The newest facular contrast from HMI. Peaks near 4 per cent at mu = 0.3 then declines, contradicting the older monotonic-rise reports. Stated by the authors to be a lower limit.

## Things I got wrong along the way, recorded so they are not repeated

I used physical radius ratios where the angular ratio is what matters. The Moon over the Sun is near 1.0 as seen from an Earth orbit, not 0.0025. The fix is in `archive/src/solarbatman/geometry.py`, which works in angular radii throughout.

I proposed transiting synthetic planets across HMI images as a novel project. It is Morris et al. (2020).

I computed the radius bias from wrong limb-darkening coefficients and presented it as new. It is Espinoza and Jordán (2015), and their number and mine agree at about 3 percent, which at least validates the machinery.

I wrote that a free baseline absorbs the effect of an unocculted starspot. It does not, because the light curve is already normalised to the dimmed star.

I read the offsets in Kostogryz et al. (2024) as relative percentages. They are absolute differences in h1 and h2 multiplied by 100. Every conclusion drawn before that was corrected was wrong by about a factor of five, including a claim that the solar reference was the limiting term. It is not. It is a fifth of the signal and a fiftieth of the library spread.

I proposed a project to test the magnetic mechanism against the solar cycle before checking what the actual bottleneck was. The bottleneck is the choice of reference library, not the precision of the solar measurement.

The common thread is that I computed before I read. The order should have been the other way round.
