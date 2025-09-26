## CalculateDiamondEfficiency.C

This macro generates efficiencies per diamond (and plots with "squares") based on the timing efficiency root file produced by automation and JSON with diamonds' coordinates.

The macro takes following arguments:

- `infile`: path to the efficiency root file (this has to be a path relative to `/eos/cms/store/group/dpg_ctpps/comm_ctpps/pps-automation`)
- `coordsJson`: path to the JSON file with diamonds' coordinates
- `runNumber`: run number
- `fourthDiamondInCyl`: 1 if there are 4 diamonds in the cylindrical RPs, 0 otherwise
- `outputPath`: path to the folder where the output files will be generated (if you want the data to be generated where you are just leave blank)


