J/A+A/666/A60       Stellar Models and Limb Darkening         (Kostogryz+, 2022)
================================================================================
Stellar limb darkening. 
A new MPS-ATLAS library for Kepler, TESS, CHEOPS, and PLATO passbands.
    Kostogryz N.M., Witzke V., Shapiro A.I., Solanki S.K., Maxted P.F.L.,
    Kurucz R.L., Gizon L.
    <Astron. Astrophys. 666, A60 (2022)>
    =2022A&A...666A..60K      (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Models, atmosphere ; Photometry ; Optical
Keywords: radiative transfer - methods: numerical - Sun: atmosphere - 
          stars: atmospheres 

Abstract:
    The detection of the first exoplanet paved the way for the era of
    transit-photometry space missions with revolutionary photometric
    precision, whose aim is to discover new exoplanetary systems around
    different types of stars. With this high precision, it is possible to
    derive the radii of exoplanets very accurately, which is crucial for
    constraining their type and composition. However, it requires an
    accurate description of their host stars, especially their
    center-to-limb variation of intensities (so-called limb darkening) as
    it affects the planet-to-star radius ratio determination.

    We aim to improve the accuracy of limb-darkening calculations for
    stars with a wide range of fundamental parameters.

    We used the recently developed 1D Merged Parallelized Simplified ATLAS
    (MPS-ATLAS) code to compute model atmosphere structures and to
    synthesize stellar limb darkening on a very fine grid of stellar
    parameters. For the computations, we utilized the most accurate
    information on chemical element abundances and mixing-length
    parameters, including convective overshoot. The stellar limb darkening
    was fitted using the two most accurate limb darkening laws: the
    power-2 and 4-parameter nonlinear laws.

    We present a new extensive library of stellar model atmospheric
    structures, the synthesized stellar limb darkening curves, and the
    coefficients of parameterized limb-darkening laws on a very fine grid
    of stellar parameters in the Kepler, TESS, CHEOPS, and PLATO
    passbands. The fine grid allows the sizable errors, introduced by the
    need to interpolate, to be overcome. Our computations of solar limb
    darkening are in a good agreement with available solar measurements at
    different view angles and wavelengths. Our computations of stellar limb
    darkening agree well with available measurements of Kepler stars. A
    new grid of stellar model structures, limb darkening, and their fitted
    coefficients in different broad passbands are provided.

Description:
    We present stellar model atmospheric structures, limb darkening, and
    limb darkening coefficients according to the description in Table 3 of
    our paper.

    files/tablee1.txt and files/tablee2.txt contain the entire grid of our
    atmospheric structures.
     First line in each structure gives the corresponding stellar parameters.

     Second line gives the number of depth points for the structure.

     Then we provide the structure which is described in Description of file

     In the end of the each structure we give the input parameters used to
     compute the structure (see description in Witzke et al.,
     2021A&A...653A..65W).

    A phython program is provided to read
     files/tablee3.txt and files/tablee4.txt.

    def read_limb_darkening_curves(filename: str, mh: float, teff: float, \
         logg: float, passband: str):
    mh_grid, teff_grid, logg_grid, i0 = np.loadtxt(filename, \
                  usecols=(0, 1, 2, 4), unpack = True)
    clv = np.loadtxt(filename, usecols=(range(5,28)))
    passbands_grid = np.array(["Kepler", "TESS", "CHEOPS", "PLATO"])
    ind_pb = np.squeeze(np.argwhere(passbands_grid==passband))
    Nmh = len(np.unique(mh_grid))
    Nteff = len(np.unique(teff_grid))
    Nlogg = len(np.unique(logg_grid))
    Npb = len(passbands_grid)
    mh_grid = mh_grid.reshape((Nmh, Nteff,Nlogg, Npb))[:, 0, 0, ind_pb]
    teff_grid = teff_grid.reshape((Nmh, Nteff, Nlogg, Npb))[0, :, 0, ind_pb]
    logg_grid = logg_grid.reshape((Nmh, Nteff, Nlogg, Npb))[0, 0, :, ind_pb]
    i0 = i0.reshape((Nmh, Nteff, Nlogg, Npb))[:, :, :, ind_pb]
    clv = clv.reshape((Nmh, Nteff, Nlogg, Npb, clv.shape[1]))[:,:,:,ind_pb,:]
    indMH = np.argmin((abs(mh_grid - mh)))
    indTeff = np.argmin((abs(teff_grid - teff)))
    indLogg = np.argmin((abs(logg_grid - logg)))
    return mh_grid[indMH], teff_grid[indTeff], logg_grid[indLogg], \
      i0[indMH, indTeff, indLogg], clv[indMH, indTeff, indLogg, :]

    This function can be called for the given stellar parameters and
    passband, and it returns the closest stellar parameters in our grid,
    absolute value of intensity at the disk center, and the corresponding
    stellar limb darkening curves at 24 disk positions.


    A phython program is provided to read files/tablee5.txt and
     files/tablee6.txt.

   def read_limb_darkening_coefficients(filename: str, mh: float, teff: float, \
   logg: float, passband: str, ld_law: str):
   mh_grid, teff_grid, logg_grid, c, alpha, a1, a2, a3, a4 =\
    np.loadtxt(filename, usecols=(0, 1, 2, 4, 5, 6, 7, 8, 9), unpack = True)
   passbands_grid = np.array(["Kepler", "TESS", "CHEOPS", "PLATO"])
   ind_pb = np.squeeze(np.argwhere(passbands_grid==passband))
   Nmh = len(np.unique(mh_grid))
   Nteff = len(np.unique(teff_grid))
   Nlogg = len(np.unique(logg_grid))
   Npb = len(passbands_grid)
   mh_grid = mh_grid.reshape((Nmh, Nteff, Nlogg, Npb))[:, 0, 0, ind_pb]
   teff_grid = teff_grid.reshape((Nmh, Nteff, Nlogg, Npb))[0, :, 0, ind_pb]
   logg_grid = logg_grid.reshape((Nmh, Nteff, Nlogg, Npb))[0, 0, :, ind_pb]
   indMH = np.argmin((np.min(mh_grid - mh)))
   indTeff = np.argmin((abs(teff_grid - teff)))
   indLogg = np.argmin((abs(logg_grid - logg)))
   if ld_law == "4_coefficients":
       a1 = a1.reshape((Nmh, Nteff, Nlogg, Npb))[:, :, :, ind_pb]
       a2 = a2.reshape((Nmh, Nteff, Nlogg, Npb))[:, :, :, ind_pb]
       a3 = a3.reshape((Nmh, Nteff, Nlogg, Npb))[:, :, :, ind_pb]
       a4 = a4.reshape((Nmh, Nteff, Nlogg, Npb))[:, :, :, ind_pb]
       return mh_grid[indMH], teff_grid[indTeff], logg_grid[indLogg], \
       a1[indMH, indTeff, indLogg], a2[indMH, indTeff, indLogg], \
       a3[indMH, indTeff, indLogg], a4[indMH, indTeff, indLogg]
   c = c.reshape((Nmh, Nteff, Nlogg, Npb))[:, :, :, ind_pb]
   alpha = alpha.reshape((Nmh, Nteff, Nlogg, Npb))[:, :, :, ind_pb]
   return mh_grid[indMH], teff_grid[indTeff], logg_grid[indLogg], \
   c[indMH, indTeff, indLogg], alpha[indMH, indTeff, indLogg]

    This function can be called for the given stellar parameters,
    passband(["Kepler","TESS","CHEOPS","PLATO"]), and limb darkening laws
    ("power-2" or "4-coefficients"), and it returns the closest stellar
    parameters in our grid and the corresponding limb-darkening
    coefficients.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
list.dat          97        6   List of the tables (table3 of the paper)
files/*            .        6   Individual tables in original format
tablee3.dat      241   136640   Set1 of stellar CLVs in different passbands
tablee4.dat      241   136640   Set2 of stellar CLVs in different passbands
tablee5.dat       79   136640   Set1 of stellar limb darkening coefficients
tablee6.dat       79   136640   Set1 of stellar limb darkening coefficients
--------------------------------------------------------------------------------

Description of file: files/table1.txt, files/table2.txt

 For each blok of data:
 ------------------------------------------------------------------------------
   Bytes  Format  Units    Label     Explanations
 ------------------------------------------------------------------------------
   2- 14  E13.12  g/cm2    rhox      Column mass
  18- 24  F6.1    K        Tgas      Gas temperature
  26- 35  E10.4   erg/cm3  Pgas      Gas pressure
  37- 46  E10.4   cm-1     ne        Electron number density
  48- 57  E10.4   cm2/g    kRoss     Mean Rosseland opacity
  59- 68  E10.4   erg/cm3  Prad      Radiation pressure
  70- 79  E10.4   cm/s     vturb     Turbulent velocity
 ------------------------------------------------------------------------------

Byte-by-byte Description of file: list.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 11  A11   ---     FileName  Name of the file, in subdirectory files (1)
  13- 91  A79   ---     Title     Description of the table
  93- 97  A5    ---     Set       Grid set
--------------------------------------------------------------------------------
Note (1): For files tablee3.txt, tablee4.txt, tablee5.txt and tablee6.txt,
 versions for the python programs.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: tablee3.dat tablee4.dat
--------------------------------------------------------------------------------
   Bytes Format Units            Label          Explanations
--------------------------------------------------------------------------------
   1-  5  F5.2  ---              M/H            [-5.0/1.5] Metallicity
   7- 10  I4    K                Teff           [3500/9000] Effective
                                                 temperature
  12- 14  F3.1  [cm/2]           logg           [3.0/5.0] gravity
  16- 21  A6    ---              passband       Instrumental passband
                                                 (CHEOPS, Kepler, PLATO, TESS)
  23- 34  E12.6 10-7W/(cm3*s*sr) I(1.0)         Disk center intensity
  36- 43  F8.6  ---              I(0.9)/I(1.0)  Limb darkening at mu=0.9
  45- 52  F8.6  ---              I(0.8)/I(1.0)  Limb darkening at mu=0.8
  54- 61  F8.6  ---              I(0.7)/I(1.0)  Limb darkening at mu=0.7
  63- 70  F8.6  ---              I(0.6)/I(1.0)  Limb darkening at mu=0.6
  72- 79  F8.6  ---              I(0.5)/I(1.0)  Limb darkening at mu=0.5
  81- 88  F8.6  ---              I(0.4)/I(1.0)  Limb darkening at mu=0.4
  90- 97  F8.6  ---              I(0.35)/I(1.0) Limb darkening at mu=0.35
  99-106  F8.6  ---              I(0.3)/I(1.0)  Limb darkening at mu=0.3
 108-115  F8.6  ---              I(0.25)/I(1.0) Limb darkening at mu=0.25
 117-124  F8.6  ---              I(0.22)/I(1.0) Limb darkening at mu=0.22
 126-133  F8.6  ---              I(0.2)/I(1.0)  Limb darkening at mu=0.2
 135-142  F8.6  ---              I(0.17)/I(1.0) Limb darkening at mu=0.17
 144-151  F8.6  ---              I(0.15)/I(1.0) Limb darkening at mu=0.15
 153-160  F8.6  ---              I(0.12)/I(1.0) Limb darkening at mu=0.12
 162-169  F8.6  ---              I(0.1)/I(1.0)  Limb darkening at mu=0.1
 171-178  F8.6  ---              I(0.09)/I(1.0) Limb darkening at mu=0.09
 180-187  F8.6  ---              I(0.08)/I(1.0) Limb darkening at mu=0.08
 189-196  F8.6  ---              I(0.07)/I(1.0) Limb darkening at mu=0.07
 198-205  F8.6  ---              I(0.06)/I(1.0) Limb darkening at mu=0.06
 207-214  F8.6  ---              I(0.05)/I(1.0) Limb darkening at mu=0.05
 216-223  F8.6  ---              I(0.03)/I(1.0) Limb darkening at mu=0.03
 225-232  F8.6  ---              I(0.02)/I(1.0) Limb darkening at mu=0.02
 234-241  F8.6  ---              I(0.01)/I(1.0) Limb darkening at mu=0.01
--------------------------------------------------------------------------------

Byte-by-byte Description of file: tablee5.dat tablee6.dat
--------------------------------------------------------------------------------
  Bytes Format Units   Label    Explanations
--------------------------------------------------------------------------------
  1-  5  F5.2  [-]     M/H      [-5.0/1.5] Metallicity
  7- 10  I4    K       Teff     [3500/9000] Effective temperature
 12- 14  F3.1  [cm/s]  logg     [3.0/5.0] gravity
 16- 21  A6    ---     passband Instrumental passband
                                 (CHEOPS, Kepler, PLATO, TESS)
 23- 30  F8.6  ---     c        coefficient of power-2 limb darkening law
 32- 39  F8.6  ---     alpha    coefficient of power-2 limb darkening law
 41- 49  F9.6  ---     a1       coefficient of 4-coefficients limb darkening law
 51- 59  F9.6  ---     a2       coefficient of 4-coefficients limb darkening law
 61- 69  F9.6  ---     a3       coefficient of 4-coefficients limb darkening law
 71- 79  F9.6  ---     a4       coefficient of 4-coefficients limb darkening law
--------------------------------------------------------------------------------

Acknowledgements:
    Nadiia Kostogryz, kostogryz(at)mps.mpg.de

================================================================================
(End)                                                         [CDS]  09-Sep-2022
