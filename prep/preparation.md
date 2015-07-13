------------------------------

Dear SISMID participant,

Ed and I are looking forward to having you in this year's [**Summer Institute in Statistics and Modeling in Infectious Diseases**](http://sismid.uw.edu) **Module 10: Simulation-based Inference for Epidemiological Dynamics**.  We'll be teaching some powerful and flexible approaches that expand the range of models that can be used in inference.

You'll be learning about these methods through computational exercises.  Accordingly, please bring a laptop.

This message is intended to make sure that you have the software you'll need for the module.  Please read the appropriate sections below, install the software as needed, and run the test scripts.  **We will not have time in Seattle to spend on installation of software, so please make sure this is all done before you arrive.**  If you run into problems, send me a note with a detailed description including the OS you’re running, the versions of `R`, `Rstudio`, and `pomp` you’re attempting to install, what you’ve done, and what problems you’ve encountered.

#### All users:

`R` and `Rstudio` are free and open-source.  You’ll need at least version 3.2.0 of `R`.  The latest version is 3.2.1, so if you need to update, go ahead and install version 3.2.1.  Source code and binaries are available on CRAN (http://cran.r-project.org).  Install the latest version of `Rstudio` (http://www.rstudio.com/products/rstudio/download/).  

Once you’ve installed these, open a session in `R` or `Rstudio` and run the following:

```
> source("http://kinglab.eeb.lsa.umich.edu/SBIED/prep/packages.R",echo=TRUE)

> source("http://kinglab.eeb.lsa.umich.edu/SBIED/prep/pompTest.R",echo=TRUE)
```

These commands run two scripts located on my website.  [The `>` is the command prompt; it’s not part of the command.  Also, depending on your email client program, you may need to replace the quotation marks with plain keyboard double quotes.]  The first script installs the packages we’ll need.  As this script runs, you may be prompted to choose a CRAN mirror: choose one near you.  The second checks whether you can work with pomp, which is the principal R package we’ll be using.  The installation script above will install the development version from my own R package repository (http://kinglab.eeb.lsa.umich.edu/R/) in preference to the stable version from CRAN.

#### Linux and unix users:

If you have trouble with either script above, make sure you have the GNU compiler collection (GCC) installed on your computer.  Linux distributions typically include this by default but it is not impossible that you have somehow avoided this.

#### MacOSX users:

So that you can compile C code and dynamically link it into an `R` session, you will need to make sure you have the Xcode app installed before running the second script above.  This is free and can be installed via the App Store or downloaded from https://developer.apple.com/xcode/downloads/.

If you have trouble with the first command trying to install `pomp` from source, receiving the error,

```
make: gfortran-4.8: No such file or directory
```

then it is likely that you do not have the necessary version of gfortran installed.  Have a look at [these instructions](./scripts/mac-fortran.html) and contact me if these don’t work for you.

#### Windows users:

You will need the ability to compile C code and dynamically link it into an `R` session.  To do this, you’ll need to install the `Rtools` suite.  Download the latest frozen version (http://cran.r-project.org/bin/windows/Rtools) and install it.  I had some problems with the unfrozen version (`Rtools33`) but none with the last frozen version (`Rtools32`).  I also had some difficulties initially with the latest version of `Rstudio` but these went away when I installed version 0.98.1103 (downloaded from https://support.rstudio.com/hc/en-us/articles/206569407-Older-Versions-of-RStudio-Desktop) and now I am having no problems with version `Rstudio` version 0.99.447.

If you use some other operating system, let me know!

In addition, if you can arrange to be able to run remote computations at your home institution, this will be helpful.  For example, if you have remote access to a multi-core computer on which to run computations, make sure it too is set up according to the instructions above and that you can log in and run computations.

Looking forward to seeing you in Seattle,

Aaron

------------------------------
