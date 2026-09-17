# Installation requirements and instructions

## Requirements

Most recent versions of the following software are always preferred, however the minimal requirements are also specified. The common software required by three main-stream operative systems are provided below, and any specific requirement is treated separately:

* [R software](http://www.r-project.org/): the minimal requirement is R version 4.0.0.

* [R packages](http://cran.r-project.org/) (it is recommended if the most recent versions of the following packages are installed):

  * `dplyr`, for extending data manipolation in R;

  * `purrr`, for extending functional programming in R;

  * `tidyr`, for extending functionalities for several data-storage formats in R;

To install the required packages, the following R code should be exectued before the installation:

```{R}
pkgnames <- c("dplyr", "purrr", "tidyr")
pkgnames <- pkgnames[!pkgnames %in% .packages(TRUE)]
if (length(pkgnames)) install.packages(pkgnames)
```

To update all the installed packages to the last version available, the following command line should be typed into an R console:

```{R}
update.packages(ask = FALSE)
```

## Package management

### Installing the stable release of the HRTnomaly package

The use of the following R command is highly suggested to install the **HRTnomaly** package:

```{R}
install.packages("HRTnomaly")
```

The other alternative to install an R package is from its source-code compressed as a tarball archive. This can be done by entering the following command into a terminal session on **Linux** and **(Mac) OS X**  :

```{bash}
R CMD INSTALL HRTnomaly_26.9.13.tar.gz
```

On **Windows**, by opening the command prompt (`cmd.exe`), it is possible to point to the proper directory with `cd`, and then install the package via `Rcmd.exe` with the following command:

```{bash}
Rcmd.exe INSTALL HRTnomaly_26.9.13.tar.gz
```

More details can be found on the "Installing packages" section of the [R-admin](https://cran.r-project.org/doc/manuals/R-admin.html) manual.

### Updating the HRTnomaly package

To update the **HRTnomaly** package, it is necessary to type the following code from the R console:

```{R}
update.packages("HRTnomaly")
```

### Removing the HRTnomaly package

To remove **HRTnomaly** from the list of R packages, it is necessary to type the following code from the R console:

```{R}
remove.packages("HRTnomaly")
```

## Third-party notice

The library [mimalloc](https://www.github.com/microsoft/mimalloc) is included in the source code folder [`src`](https://github.com/drwolf85/HRTnomaly/tree/master/src) and redistributed under [MTI license](https://github.com/drwolf85/HRTnomaly/blob/master/src/mimalloc/LICENSE).

Snippets of code from the [pcg-c-basic](https://github.com/imneme/pcg-c-basic) is included in the file [`src/post_thresh.c`](https://github.com/drwolf85/HRTnomaly/tree/master/src/post_thresh.c) and redistributed under [APACHE-2.0 license](https://github.com/drwolf85/HRTnomaly/blob/master/inst/LICENSE-PCG32.txt).
