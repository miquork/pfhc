# pfhc
Particle Flow Hadron Correction (PFHC)

Quick how-to:
=============

`root -l -b -q mk_compile.C`

`root -l -b -q mk_piongun.C` [edit file list if needed]

`root -l -b -q drawPionGun.C+g` [if needed]

This will produce `piongun.root` used by `drawPionGun.root`. The latter will store results in `drawPionGun.root` and `piongun.txt`, in addition to many plots in pdf folder.

The text file can be added with `#include "piongun.txt"` in `PFEnergyCalibrationFromMikko.C` to test closure of the new corrections by re-running the code above after switching on `bool applyPFEC_Neutral = true;` in `piongun.C`.