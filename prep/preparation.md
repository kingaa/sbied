------------------------------

# Software installation instructions

Please read the appropriate sections below, install the software as needed, and run the test scripts.  If you run into problems, send a note to *kingaa@umich.edu* with a detailed description of the problem you've encountered, including the OS you’re running, the versions of **R**, **Rstudio**, and **pomp** you’re attempting to install, what you’ve done, and what problems you’ve encountered.

#### All users:

**R** and **Rstudio** are free and open-source.  You’ll need at least version 3.2.0 of **R**.  The latest version is 3.2.1, so if you need to update, go ahead and install version 3.2.1.  Source code and binaries are available on CRAN (http://cran.r-project.org).  Install the latest version of **Rstudio** (http://www.rstudio.com/products/rstudio/download/).  

Once you’ve installed these, open a session in **R** or **Rstudio** and run the following:

```
> update.packages()
> source("http://kingaa.github.io/sbied/prep/packages.R",echo=TRUE)
> source("http://kinglab.eeb.lsa.umich.edu/SBIED/prep/pompTest.R",echo=TRUE)
```

[The `>` is the command prompt; it’s not part of the command.  Also, depending on your email client program, you may need to replace the quotation marks with plain keyboard double quotes.]  The first command updates your installed packages.  You may be prompted to specify a CRAN mirror: choose one near you.  The second command runs a script on my website.  It will install some needed packages if these are not already installed on your system.  The third command will attempt to install **pomp**, the principal R package we’ll be using, from source code and will check whether you can work with it.  The `pompTest.R` script will attempt to install a version that is more recent than that which is on CRAN.  The latter will not be enough for our purposes.

If the final command fails, try the following:
```
> source("http://kingaa.github.io/sbied/prep/hello.R",echo=TRUE)
```
If this fails to give the "Hello world!" message, you will need to follow the instructions below that correspond to your OS.

#### Linux and unix users:

If you have trouble with either script above, make sure you have the GNU compiler collection (GCC) installed on your computer.  Linux distributions typically include this by default but it is not impossible that you have somehow avoided this.

#### MacOSX users:

So that you can compile C code and dynamically link it into an **R** session, you will need to make sure you have the Xcode app installed before running the second script above.  This is free and can be installed via the App Store or downloaded from [https://developer.apple.com/xcode/downloads/].

If you have trouble with the first command trying to install `pomp` from source, receiving the error,

```
make: gfortran-4.8: No such file or directory
```

then it is likely that you do not have the necessary version of gfortran installed.  Have a look at [these instructions](./mac-fortran.html) and contact me if these don’t work for you.

#### Windows users:

You will need the ability to compile C code and dynamically link it into an **R** session.  To do this, you’ll need to install the `Rtools` suite.  Download the latest frozen version ([http://cran.r-project.org/bin/windows/Rtools]) and install it.  I had some problems with the unfrozen version (**Rtools33**) but none with the last frozen version (**Rtools32**).  I also had some difficulties initially with the latest version of **Rstudio** but these went away when I installed version 0.98.1103 (downloaded from [here](https://support.rstudio.com/hc/en-us/articles/206569407-Older-Versions-of-RStudio-Desktop)) and now I am having no problems with version **Rstudio** version 0.99.447.

If you use some other operating system, let me know!

------------------------------
