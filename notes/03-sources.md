# Sources

Every paper this audit rests on, with the reading notes taken from it. The PDFs themselves are not in this repository, because they are published copies that are not mine to redistribute. Each entry carries an arXiv number or a journal reference so you can fetch your own.

The last section records how much of each paper was read, in full or in part or in abstract only, so that no claim here rests on more reading than actually happened.

## The result being audited

**Kostogryz, Shapiro, Witzke, Cameron, Gizon, Krivova, Ludwig, Maxted, Seager, Solanki, Valenti (2024), Nature Astronomy 8, 929.** *Magnetic origin of the discrepancy between stellar limb-darkening models and observations.* arXiv:2403.00118.

The paper everything here responds to. Non-magnetic model atmospheres predict limb darkening that is too steep, and MURaM simulations with small-scale surface fields reproduce the observations. The offsets used in `data/kostogryz2024_offsets.csv` come from the Methods section, under "Limb darkening without magnetic fields".

The definitions to know are h1' = I(2/3)/I(1) and h2' = I(2/3)/I(1) − I(1/3)/I(1), where the argument is mu.

**Maxted (2023), MNRAS 519, 3723.** *Limb darkening measurements from TESS and Kepler light curves of transiting exoplanets.* arXiv:2212.09117.

This is where the observed offsets actually come from. Table 3 gives the same measurement against nine libraries and is transcribed into `data/maxted2023_table3.csv`. Section 4.3.1 gives the sample means and the interpolation method. Section 4.2 reports the reanalysis of sixteen stars that supplies the h2 analysis systematic, which is the centre of this audit. Section 4.3.2 is the exchange with the referee about solar calibration, which is the origin of the facular objection quantified in `analysis/a03`.

His sign convention is model minus observed, the reverse of Kostogryz. Confirming that his Table 3 values are the same numbers Kostogryz quote as percentages is what fixes the units for everything in this repository.

**Kostogryz, Witzke, Shapiro, Solanki, Maxted, Kurucz, Gizon (2022), A&A 666, A60.** *Stellar limb darkening. A new MPS-ATLAS library for Kepler, TESS, CHEOPS, and PLATO passbands.* arXiv:2206.06641.

The library the 2024 comparison is measured against, and the reason `data/kostogryz2022/` exists. Section 5.1 compares it to Neckel and Labs. Section 5.2 is the important one, because it describes NL_Kepler, the Neckel and Labs solar profile carried into the Kepler band. That is exactly the object needed to rank the other libraries against the Sun, and it was built for their own models only.

## The measurements the audit leans on

**Norris, Beeck, Unruh, Solanki, Krivova, Yeo (2017), A&A 605, A45.** *Spectral variability of photospheric radiation due to faculae I. The Sun and Sun-like stars.* arXiv:1705.04455.

The source of the MURaM 100 G numbers Maxted quotes. Figure 4 plots the Neckel and Labs polynomials against a shaded band running from field-free to 100 G MURaM, at 391, 611, 1097.5 and 1597.5 nm. Maxted's +0.007 and −0.005 are the width of that band read off the plot, and they are not tabulated anywhere. Cite them as read from the figure.

The text beside that figure matters as much as the figure. They say the measured solar gradient is slightly shallower than field-free MURaM, and that the difference may in part be because we never observe the truly field-free quiet Sun. That is the Kostogryz et al. (2024) effect, stated for the Sun, seven years earlier, by an author list including Krivova and Solanki.

**Pietrow, Sumra, Petit dit de la Roche, Loessnitz, Denker, Assmus (2026).** *Center-to-limb variations of solar active regions.* arXiv:2606.28887.

The newest measured facular and network contrast, from HMI at 6173 Å with stray light corrected. Faculae peak near 4 per cent at mu around 0.3 and then decline towards the limb, which contradicts older reports of a monotonic rise. The authors state their values are a lower limit. This supplies the contrast anchors in `data/facular_contrast.csv`.

Their closing point is worth noting for the longer project. PHOENIX and ATLAS models with adjusted effective temperature fail to reproduce the observed centre-to-limb behaviour of active regions at all.

**Yeo, Solanki, Krivova (2013), A&A 550, A95.** *Intensity contrast of solar network and faculae.* arXiv:1302.1442.

SDO/HMI at 6173 Å. Confirms that continuum contrast is weakest near disc centre, peaks, and then declines towards the limb, and that the peak rises with magnetogram signal. The peak contrast is plotted in their Figure 14 against ⟨Bl⟩/mu and is not stated as a number in the text. The 10 per cent figure attributed to this paper is Pietrow et al.'s reading of that figure, and `analysis/a03` does not depend on it.

**Hestroffer and Magnan (1998), A&A 333, 338.**

Table 1 is the source of the two solar coefficient sets in `data/solar_clv_polynomials.csv`, which come originally from Pierce and Slaughter (1977) and Neckel and Labs (1994). The paper is worth reading for its argument that polynomial fits to the solar profile are unreliable near the limb, which is the same point this project keeps running into.

## Context on transit limb darkening

These are the papers that stop this project from rediscovering published results. `notes/02-literature-audit.md` records what each one already settles.

**Howarth (2011), MNRAS 418, 1165.** arXiv:1106.4659. *On stellar limb darkening and exoplanetary transits.* Shows that limb-darkening coefficients fitted to a transit light curve are not the same quantity as coefficients fitted to a model intensity profile, and introduces SPAM as the fix.

**Espinoza and Jordán (2015), MNRAS 450, 1879.** arXiv:1503.07020. *Limb darkening and exoplanets: testing stellar model atmospheres and identifying biases in transit parameters.* Section 4.3 and Figures 10 to 14 quantify the bias in the radius ratio, the scaled semi-major axis and the inclination from fixing or floating the coefficients, at impact parameters of 0.3 and 0.8. Up to 3 per cent in the radius ratio. Anything computed here about limb-darkening-induced transit parameter bias should be checked against this paper first, because it is probably already in it.

**Morris, Bobra, Agol, Lee, Hawley (2020), MNRAS 493, 5489.** arXiv:2002.08072. *The stellar variability noise floor for transiting exoplanet photometry with PLATO.* Transits synthetic planets across real HMI continuum images to get the granulation and oscillation noise floor, and finds a 3.6 per cent radius uncertainty driven by a degeneracy between the radius ratio, limb darkening and impact parameter. Read this before proposing anything that puts a synthetic planet across an HMI image.

**Kostogryz, Shapiro, Witzke, Bhatia, Solanki, Kuhlemann, Vasilyev, Unruh (2026).** arXiv:2606.21912. *Effect of surface magnetic fields on limb darkening in main-sequence stars.* The follow-up to the 2024 paper. It maps the magnetic effect from 3200 to 6800 K and across metallicity, and it reports that the effect grows towards hotter and more metal-rich stars, which is the direction the Maxted sample lies in. It does not convert the effect into transit parameter biases and it does not test against the Sun. The paper says a public database of synthetic spectra at ten disc positions was released, and finding that database is the first thing to do for `analysis/a07`.

**Maxted (2018), A&A 616, A39.** Introduces the h1, h2 reparametrisation of the power-2 law, and suggests the magnetic explanation for the offset.

## The route for the next step

**Stray light correction for the Helioseismic and Magnetic Imager.** arXiv:2511.13348.

The method for using lunar transits to characterise the HMI point spread function. An artificial lunar transit image is built from the ephemeris and an assumed limb-darkening profile, convolved with a candidate point spread function, and compared against the observed transit. This is the technical route for measuring a modern solar profile, and `archive/` holds the forward model it would need.

There is a circularity here to be careful about. The method assumes a limb-darkening profile in order to derive the point spread function, and the point spread function is what is wanted in order to measure the limb-darkening profile. Solving for both together is a real part of that work.

## Not obtained, and what each would settle

**Neckel and Labs (1994), Solar Physics 153, 91.** The full wavelength tables. These are what the passband conversion in `analysis/a07` needs, and without them the solar comparison cannot be carried into the Kepler and TESS bands.

**Pierce and Slaughter (1977), Solar Physics 52, 179.** The other measured solar profile, for the same reason.

**Kim et al. (2017), JASS 34, 99.** *Variation in solar limb darkening coefficient estimated from solar images taken by SOHO and SDO.* Measured a centre-to-limb intensity ratio from MDI and HMI over 1996 to 2016 and attributed the SDO-era trend to CCD ageing. This is the closest existing attempt at the solar cycle measurement, and it shows why the point spread function has to be handled first.

## How much of each was read

Read in full, or in every section that bears on this audit: Kostogryz et al. (2024), Kostogryz et al. (2022), Maxted (2023), Espinoza and Jordán (2015), Pietrow et al. (2026), Hestroffer and Magnan (1998).

Read in the relevant sections only: Norris et al. (2017), Yeo et al. (2013), Howarth (2011).

Read in abstract only: Kostogryz et al. (2026), Maxted (2018), Morris et al. (2020), Kim et al. (2017).

Not obtained: Neckel and Labs (1994), Pierce and Slaughter (1977), whose coefficients are used here as transcribed by Hestroffer and Magnan (1998).
