
# SkyGlyMerger

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of SkyGlyMerger is to merge SKyline data with GlyCounter data...


## Installation

1.  Install [R version 4.5.3](https://cran.r-project.org/bin/windows/base/old/) and [Rstudio](https://posit.co/download/rstudio-desktop/) on your computer.
    (R 4.5.3 can be installed alongside other R versions).
2.  Install [Rtools 4.5](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html) using the official _Rtools45_ installer and keep the default settings.
    This is required to build some packages from source.
3.  Configure RStudio to use R 4.5.3 (`Tools → Global Options → General → R version → Change…`). Apply the changes and close RStudio.
4.  Download the source code of the `main` branch as a zip file, then unzip and store the `SkyGlyMerger-main` folder somewhere.
5.  Double click `SkyGlyMerger.Rproj`, this will open the project in RStudio. You will see a message indicating that the `renv` package was installed.
6.  In the R console, run `renv::restore()`. This will download and install all required R packages exactly as specified for this project.\
   ⏳ *This may take several minutes.*
8.  In RStudio, open the file `dev/run_dev.R`.
9.  With `dev/run_dev.R` open, press `Ctrl+Shift+Enter` to run the dashboard.\
    If prompted to install extra packages, confirm the installation.
10. After these steps, the application should start automatically in your RStudio session.
