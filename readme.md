Alternative Stable States of the Forest Mycobiome
================

This R project contains all code to replicate analyses and figures
presented in Averill et al. 2022. Alternative Stable States of the
Forest Mycobiome. *Nature Ecology & Evolution*.

The file `paths.r` is sourced at the top of all analysis and figures
scripts. This file tracks all file paths within the project. To run this
code, the user needs to change the path object `storage.dir` within the
paths.r file to match the location of the data directory on their
computer. If you run this analysis on multiple computers you will need a
few conditional statements to specify where the data directory is, as
well as a way to sync files across the computers (more on this below).
Furthermore, you will need a copy of the US Forest Inventory and
Analysis Database (FIA database), version 7. This project uses a pSQL
version of this database. The path to this directory of files needs to
be specified in `paths.r`. I have kept the FIA7 directory out of the
overall storage directory because the FIA database is huge.

This project also depends on several custom functions housed in
`project_functions/` as well as a few small reference data files, housed
in `required_products_utilities/`. Finally, this analysis depends on the
CASTNET nitrogen deposition product, and the worldclim2 climate product.
Paths to the rasters from these products need to be specified when using
the functions `worldclim2_grab()` and `extract_ndep()`.

The directory `1._data_construction/` contains all code to work up data
within the FIA7 database into analysis ready products. These scripts are
numbered and must be run in order.

The directory `2._data_analysis/` contains all code to replicate
analyses and simulations described in the paper. These scripts are
numbered and must be run in order. The hystersis simulation is
computationally intensive, and takes several hours to run on a 28-core
computing cluster. Some scripts have a submission script, indicated by
`3q._script_name.r`, which submits the script `3._script_name.r`. These
are design for my particular computing cluster, and you will need to
generate your own batch submission scripts for your computing resources
if your computer cannot run these jobs interactively.

All analysis output is stored in the data directory specified in the
paths.r script. All figures are stored within the project repository sub
directory `figures/`. This subdirectory is ignored for the purposes of
git tracking.

I run this project on multiple computers, and use in the files within
`transfer_scripts` to sync files across computers. This is a pretty
simple way to do this, and file transfers need to be triggered
“manually”. There are more sophisticated ways to do this (check out
Rsync), but this way works. This won’t be necessary if you perform all
computation on a single machine.

This project was built under R version 3.6.1 and depends on the
following packages: `betareg` `boot` `data.table` `doParallel` `ggalt`
`ggplot2` `ggpubr` `magick` `mgcv` `purrr` `raster`
