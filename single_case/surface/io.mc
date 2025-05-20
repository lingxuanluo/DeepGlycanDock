##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue May 20 00:53:46 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path single_case_temp
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 805 atoms
Valist_getStatistics:  Max atom coordinate:  (138.487, 109.614, 120.147)
Valist_getStatistics:  Min atom coordinate:  (105.226, 82.323, 95.117)
Valist_getStatistics:  Molecule center:  (121.856, 95.9685, 107.632)
NOsh_setupCalcMGAUTO(/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1855):  coarse grid center = 121.856 95.9685 107.632
NOsh_setupCalcMGAUTO(/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1860):  fine grid center = 121.856 95.9685 107.632
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.633002, 0.523086, 0.494098
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.580688, 0.516031, 0.494098
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.917355, 0.986512, 1 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1972):  coarse mesh center = 121.856 95.9685 107.632
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 152.241 121.077 131.349
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = 91.4724 70.8604 83.9153
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 149.73 120.738 131.349
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = 93.9835 71.199 83.9153
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 149.73 120.738 131.349
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = 93.9835 71.199 83.9153
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 21.2206
Vpbe_ctor2:  solute dimensions = 35.746 x 29.539 x 27.902
Vpbe_ctor2:  solute charge = 2
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 71 x 59 x 55 table
Vclist_ctor2:  Using 71 x 59 x 55 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (44.337, 38.367, 36.106)
Vclist_setupGrid:  Grid lower corner = (99.688, 76.785, 89.579)
Vclist_assignAtoms:  Have 1176764 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.044609
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.513130e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 2.549000e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.749900e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 3.730720e-01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.427325e-01
Vprtstp: contraction number = 1.427325e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 2.111670e-02
Vprtstp: contraction number = 1.479460e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 3.617577e-03
Vprtstp: contraction number = 1.713136e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 6.876070e-04
Vprtstp: contraction number = 1.900739e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 1.396392e-04
Vprtstp: contraction number = 2.030800e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 2.976058e-05
Vprtstp: contraction number = 2.131248e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 6.612148e-06
Vprtstp: contraction number = 2.221780e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 1.522015e-06
Vprtstp: contraction number = 2.301846e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 3.613235e-07
Vprtstp: contraction number = 2.373982e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 5.078180e-01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 6.146970e-01
Vpmg_setPart:  lower corner = (91.4724, 70.8604, 83.9153)
Vpmg_setPart:  upper corner = (152.241, 121.077, 131.349)
Vpmg_setPart:  actual minima = (91.4724, 70.8604, 83.9153)
Vpmg_setPart:  actual maxima = (152.241, 121.077, 131.349)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 2.408907513214E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.815000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 21.2206
Vpbe_ctor2:  solute dimensions = 35.746 x 29.539 x 27.902
Vpbe_ctor2:  solute charge = 2
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 71 x 59 x 55 table
Vclist_ctor2:  Using 71 x 59 x 55 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (44.337, 38.367, 36.106)
Vclist_setupGrid:  Grid lower corner = (99.688, 76.785, 89.579)
Vclist_assignAtoms:  Have 1176764 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = 93.9835, 71.199, 83.9153
VPMG::focusFillBound -- New mesh maxs = 149.73, 120.738, 131.349
VPMG::focusFillBound -- Old mesh mins = 91.4724, 70.8604, 83.9153
VPMG::focusFillBound -- Old mesh maxs = 152.241, 121.077, 131.349
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (93.9835, 71.199, 83.9153)
Vpmg_setPart:  upper corner = (149.73, 120.738, 131.349)
Vpmg_setPart:  actual minima = (91.4724, 70.8604, 83.9153)
Vpmg_setPart:  actual maxima = (152.241, 121.077, 131.349)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (93.9835, 71.199, 83.9153)
VPMG::extEnergy    Disj part upper corner = (149.73, 120.738, 131.349)
VPMG::extEnergy    Old lower corner = (91.4724, 70.8604, 83.9153)
VPMG::extEnergy    Old upper corner = (152.241, 121.077, 131.349)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 0.06663 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.046120
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.825940e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 2.548700e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.738100e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.276473e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.437504e-01
Vprtstp: contraction number = 1.437504e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.990559e-02
Vprtstp: contraction number = 1.384732e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 3.040729e-03
Vprtstp: contraction number = 1.527576e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 5.015294e-04
Vprtstp: contraction number = 1.649372e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 8.832023e-05
Vprtstp: contraction number = 1.761018e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 1.645356e-05
Vprtstp: contraction number = 1.862944e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 3.215539e-06
Vprtstp: contraction number = 1.954311e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 6.543835e-07
Vprtstp: contraction number = 2.035067e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 4.586080e-01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 5.653030e-01
Vpmg_setPart:  lower corner = (93.9835, 71.199, 83.9153)
Vpmg_setPart:  upper corner = (149.73, 120.738, 131.349)
Vpmg_setPart:  actual minima = (93.9835, 71.199, 83.9153)
Vpmg_setPart:  actual maxima = (149.73, 120.738, 131.349)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 2.553552786378E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.826000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 2.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
printEnergy:  Performing global reduction (sum)
Vcom_reduce:  Not compiled with MPI, doing simple copy.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 2.031592e+00
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue May 20 00:53:48 2025

##############################################################################
Vgrid_readDX:  Grid dimensions 97 x 97 x 97 grid
Vgrid_readDX:  Grid origin = (93.9835, 71.199, 83.9153)
Vgrid_readDX:  Grid spacings = (0.580688, 0.516031, 0.494098)
Vgrid_readDX:  allocating 97 x 97 x 97 doubles for storage
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue May 20 00:55:53 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path single_case_temp
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 805 atoms
Valist_getStatistics:  Max atom coordinate:  (138.487, 109.614, 120.147)
Valist_getStatistics:  Min atom coordinate:  (105.226, 82.323, 95.117)
Valist_getStatistics:  Molecule center:  (121.856, 95.9685, 107.632)
NOsh_setupCalcMGAUTO(/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1855):  coarse grid center = 121.856 95.9685 107.632
NOsh_setupCalcMGAUTO(/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1860):  fine grid center = 121.856 95.9685 107.632
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.633002, 0.523086, 0.494098
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.580688, 0.516031, 0.494098
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.917355, 0.986512, 1 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1972):  coarse mesh center = 121.856 95.9685 107.632
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 152.241 121.077 131.349
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = 91.4724 70.8604 83.9153
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 149.73 120.738 131.349
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = 93.9835 71.199 83.9153
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 149.73 120.738 131.349
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = 93.9835 71.199 83.9153
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 21.2206
Vpbe_ctor2:  solute dimensions = 35.746 x 29.539 x 27.902
Vpbe_ctor2:  solute charge = 2
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 71 x 59 x 55 table
Vclist_ctor2:  Using 71 x 59 x 55 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (44.337, 38.367, 36.106)
Vclist_setupGrid:  Grid lower corner = (99.688, 76.785, 89.579)
Vclist_assignAtoms:  Have 1176764 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.045394
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.502260e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 2.488100e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.650600e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 3.692230e-01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.427325e-01
Vprtstp: contraction number = 1.427325e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 2.111670e-02
Vprtstp: contraction number = 1.479460e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 3.617577e-03
Vprtstp: contraction number = 1.713136e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 6.876070e-04
Vprtstp: contraction number = 1.900739e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 1.396392e-04
Vprtstp: contraction number = 2.030800e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 2.976058e-05
Vprtstp: contraction number = 2.131248e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 6.612148e-06
Vprtstp: contraction number = 2.221780e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 1.522015e-06
Vprtstp: contraction number = 2.301846e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 3.613235e-07
Vprtstp: contraction number = 2.373982e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 5.200710e-01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 6.251980e-01
Vpmg_setPart:  lower corner = (91.4724, 70.8604, 83.9153)
Vpmg_setPart:  upper corner = (152.241, 121.077, 131.349)
Vpmg_setPart:  actual minima = (91.4724, 70.8604, 83.9153)
Vpmg_setPart:  actual maxima = (152.241, 121.077, 131.349)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 2.408907513214E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.797000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 21.2206
Vpbe_ctor2:  solute dimensions = 35.746 x 29.539 x 27.902
Vpbe_ctor2:  solute charge = 2
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 71 x 59 x 55 table
Vclist_ctor2:  Using 71 x 59 x 55 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (44.337, 38.367, 36.106)
Vclist_setupGrid:  Grid lower corner = (99.688, 76.785, 89.579)
Vclist_assignAtoms:  Have 1176764 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = 93.9835, 71.199, 83.9153
VPMG::focusFillBound -- New mesh maxs = 149.73, 120.738, 131.349
VPMG::focusFillBound -- Old mesh mins = 91.4724, 70.8604, 83.9153
VPMG::focusFillBound -- Old mesh maxs = 152.241, 121.077, 131.349
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (93.9835, 71.199, 83.9153)
Vpmg_setPart:  upper corner = (149.73, 120.738, 131.349)
Vpmg_setPart:  actual minima = (91.4724, 70.8604, 83.9153)
Vpmg_setPart:  actual maxima = (152.241, 121.077, 131.349)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (93.9835, 71.199, 83.9153)
VPMG::extEnergy    Disj part upper corner = (149.73, 120.738, 131.349)
VPMG::extEnergy    Old lower corner = (91.4724, 70.8604, 83.9153)
VPMG::extEnergy    Old upper corner = (152.241, 121.077, 131.349)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 0.06663 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.046431
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.803730e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 2.483100e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.650000e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.282743e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.437504e-01
Vprtstp: contraction number = 1.437504e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.990559e-02
Vprtstp: contraction number = 1.384732e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 3.040729e-03
Vprtstp: contraction number = 1.527576e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 5.015294e-04
Vprtstp: contraction number = 1.649372e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 8.832023e-05
Vprtstp: contraction number = 1.761018e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 1.645356e-05
Vprtstp: contraction number = 1.862944e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 3.215539e-06
Vprtstp: contraction number = 1.954311e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 6.543835e-07
Vprtstp: contraction number = 2.035067e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 4.622890e-01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 5.673070e-01
Vpmg_setPart:  lower corner = (93.9835, 71.199, 83.9153)
Vpmg_setPart:  upper corner = (149.73, 120.738, 131.349)
Vpmg_setPart:  actual minima = (93.9835, 71.199, 83.9153)
Vpmg_setPart:  actual maxima = (149.73, 120.738, 131.349)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 2.553552786378E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.810000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 2.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
printEnergy:  Performing global reduction (sum)
Vcom_reduce:  Not compiled with MPI, doing simple copy.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 2.033781e+00
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue May 20 00:55:55 2025

##############################################################################
Vgrid_readDX:  Grid dimensions 97 x 97 x 97 grid
Vgrid_readDX:  Grid origin = (93.9835, 71.199, 83.9153)
Vgrid_readDX:  Grid spacings = (0.580688, 0.516031, 0.494098)
Vgrid_readDX:  allocating 97 x 97 x 97 doubles for storage
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue May 20 00:58:05 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path single_case_temp
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing file (got QUIT)
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue May 20 09:30:19 2025

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path 8ufh_y75_prepared_temp
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 4, dime = (161, 193, 225)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 3142 atoms
Valist_getStatistics:  Max atom coordinate:  (143.545, 166.042, 170.467)
Valist_getStatistics:  Min atom coordinate:  (87.616, 89.705, 84.535)
Valist_getStatistics:  Molecule center:  (115.581, 127.874, 127.501)
NOsh_setupCalcMGAUTO(/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1855):  coarse grid center = 115.581 127.874 127.501
NOsh_setupCalcMGAUTO(/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1860):  fine grid center = 115.581 127.874 127.501
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.616664, 0.697407, 0.666711
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.487744, 0.514406, 0.481469
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.790939, 0.737598, 0.722155 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1972):  coarse mesh center = 115.581 127.874 127.501
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 164.914 194.825 202.173
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = 66.2473 60.9224 52.8294
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 154.6 177.257 181.425
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = 76.561 78.4905 73.5765
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 154.6 177.257 181.425
NOsh_setupCalcMGAUTO (/install/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = 76.561 78.4905 73.5765
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 56.5893
Vpbe_ctor2:  solute dimensions = 58.039 x 78.766 x 87.849
Vpbe_ctor2:  solute charge = 18
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (67.005, 87.413, 97.008)
Vclist_setupGrid:  Grid lower corner = (82.078, 84.167, 78.997)
Vclist_assignAtoms:  Have 1282217 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.178621
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.031595e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (161, 193, 225)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 2.056900e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (081, 097, 113)
Vbuildops: Galer: (041, 049, 057)
Vbuildops: Galer: (021, 025, 029)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 1.077641e+00
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.458694e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.116769e-01
Vprtstp: contraction number = 1.116769e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.276598e-02
Vprtstp: contraction number = 1.143117e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.547555e-03
Vprtstp: contraction number = 1.212249e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.982197e-04
Vprtstp: contraction number = 1.280857e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.738260e-05
Vprtstp: contraction number = 1.381426e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 4.130212e-06
Vprtstp: contraction number = 1.508335e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 6.809082e-07
Vprtstp: contraction number = 1.648604e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 4.551390e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 5.951462e+00
Vpmg_setPart:  lower corner = (66.2473, 60.9224, 52.8294)
Vpmg_setPart:  upper corner = (164.914, 194.825, 202.173)
Vpmg_setPart:  actual minima = (66.2473, 60.9224, 52.8294)
Vpmg_setPart:  actual maxima = (164.914, 194.825, 202.173)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 6.148808315566E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.561100e-02
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 56.5893
Vpbe_ctor2:  solute dimensions = 58.039 x 78.766 x 87.849
Vpbe_ctor2:  solute charge = 18
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 75 table
Vclist_ctor2:  Using 75 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (67.005, 87.413, 97.008)
Vclist_setupGrid:  Grid lower corner = (82.078, 84.167, 78.997)
Vclist_assignAtoms:  Have 1282217 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = 76.561, 78.4905, 73.5765
VPMG::focusFillBound -- New mesh maxs = 154.6, 177.257, 181.425
VPMG::focusFillBound -- Old mesh mins = 66.2473, 60.9224, 52.8294
VPMG::focusFillBound -- Old mesh maxs = 164.914, 194.825, 202.173
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (76.561, 78.4905, 73.5765)
Vpmg_setPart:  upper corner = (154.6, 177.257, 181.425)
Vpmg_setPart:  actual minima = (66.2473, 60.9224, 52.8294)
Vpmg_setPart:  actual maxima = (164.914, 194.825, 202.173)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (76.561, 78.4905, 73.5765)
VPMG::extEnergy    Disj part upper corner = (154.6, 177.257, 181.425)
VPMG::extEnergy    Old lower corner = (66.2473, 60.9224, 52.8294)
VPMG::extEnergy    Old upper corner = (164.914, 194.825, 202.173)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 5.42537 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.153781
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.317973e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (161, 193, 225)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.759030e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (081, 097, 113)
Vbuildops: Galer: (041, 049, 057)
Vbuildops: Galer: (021, 025, 029)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 9.430970e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 9.619195e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.235109e-01
Vprtstp: contraction number = 1.235109e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.443475e-02
Vprtstp: contraction number = 1.168702e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.769091e-03
Vprtstp: contraction number = 1.225578e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 2.314125e-04
Vprtstp: contraction number = 1.308087e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 3.324570e-05
Vprtstp: contraction number = 1.436642e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 5.321610e-06
Vprtstp: contraction number = 1.600691e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 9.424800e-07
Vprtstp: contraction number = 1.771043e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 4.536115e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 5.762683e+00
Vpmg_setPart:  lower corner = (76.561, 78.4905, 73.5765)
Vpmg_setPart:  upper corner = (154.6, 177.257, 181.425)
Vpmg_setPart:  actual minima = (76.561, 78.4905, 73.5765)
Vpmg_setPart:  actual maxima = (154.6, 177.257, 181.425)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 1.069613914748E+05 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.461200e-02
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 0.000000e+00
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
printEnergy:  Performing global reduction (sum)
Vcom_reduce:  Not compiled with MPI, doing simple copy.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 1.605484e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Tue May 20 09:30:35 2025

##############################################################################
Vgrid_readDX:  Grid dimensions 161 x 193 x 225 grid
Vgrid_readDX:  Grid origin = (76.561, 78.4905, 73.5765)
Vgrid_readDX:  Grid spacings = (0.487744, 0.514406, 0.481469)
Vgrid_readDX:  allocating 161 x 193 x 225 doubles for storage
