# pfhc
Particle Flow Hadron Correction (PFHC)

Quick how-to:
=============

`root -l -b -q mk_compile.C`

`root -l -b -q mk_piongun.C` [edit file list if needed]

After filling histograms and storing them in `piongun.root`,  it runs these plotting steps that can also be repeated manually:

`root -l -b -q drawPiongun.C+g`

`root -l -b -q drawPiongunResolution.C+g`

`root -l -b -q drawPiongunEfficiency.C+g`

The `drawPiongun.root` takes `piongun.root` as input and creates `drawPiongun.root` and `piongut.txt`, in addition to many plots in `pdf` folder.

The text file can be added with `#include "piongun.txt"` in `PFEnergyCalibrationFromMikko.C` to test closure of the new corrections by re-running the code above after switching on `bool applyPFEC_Neutral = true;` in `piongun.C`.