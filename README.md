This is an add-on package to the R package [grpSLOPE](https://cran.r-project.org/package=grpSLOPE). It contains Monte Carlo based methods for the estimation of the regularizing sequence.

**This package contains exactly one function, `lambdaMC`. This package is designed to be used in conjunction with the R package [grpSLOPE](https://cran.r-project.org/package=grpSLOPE).**

## Installation

Your R configuration must allow for a working Rcpp. This is generally not a problem on Unix/Linux, but setting it up on Windows may require some work.

The easiest way to install the latest development version of `grpSLOPEMC` is by using the R package `devtools`. Just open up an R session and run:

```R
# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("agisga/grpSLOPEMC")
```

If you don't want to use `devtools`, you can install `grpSLOPEMC` by downloading the source code and then following these steps:

0. Install the R packages `Rcpp` and `RcppEigen` if you don't have them installed already.
1. Go to the directory that contains the `grpSLOPEMC` directory (which contains the `grpSLOPEMC` source code).
2. Open an R session and run `Rcpp::compileAttributes("./grpSLOPEMC")`. Then quit R.
3. Run `R CMD build grpSLOPEMC`. You should then have a file like `grpSLOPEMC_0.1.0.tar.gz`.
4. Run `R CMD INSTALL grpSLOPEMC_0.1.0.tar.gz` to install the package.

## Contributing

### Code style

Variable names are all lower case with words separated by dots.
Function names begin with a lower case letter and are written in camel case.
Constants names are all caps.
Otherwise, I try to follow [Google's R style guide](https://google.github.io/styleguide/Rguide.xml).

### Workflow

0. Modify the code.
1. Open `grpSLOPEMC.Rproj` with RStudio.
2. Run `devtools::document()`.
3. Do "Build and Reload" from the menu (or CTRL-Shift-B).
4. Do `devtools::test()` to run the unit tests.
5. Install with `devtools::install()`
