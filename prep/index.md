---
title: Software installation instructions
output:
  html_document:
    includes:
      after_body:
        - ../_includes/lesson_bottom.html
        - ../_includes/license.html        
---

<style type="text/css">
div .nb {
	background-color: #ffeca3;
	border-style: solid;
	border-width: 2;
	border-color: #00274c;
	padding: 1em;
}
hr {
	border-width: 3;
	border-color: #00274c;
}
</style>

------------------------------

Please read the appropriate sections below, which give specific instructions for installing and testing the software we will be using.
First follow the instructions for "all users", then those for your specific operating system (OS).

<div class="nb"> 

**Important note:** If you run into problems, send a direct message to `@Aaron A. King` on the course Slack channel with a detailed description of the problem you've encountered.
In this message, **be certain to include all of the following information**:

- the operating system you’re running,
- the version numbers of **R**, **Rstudio**, and **pomp** you’re attempting to install,
- what command(s) you've executed, and
- what error messages you've received.

In particular, it is often easiest to send a screenshot or transcript showing the commands you've entered and the error messages you've received.
In an **R** session, you can run 
```
> source("https://kingaa.github.io/scripts/diagnostics.R")
```
to get information on the operating system and version numbers of **R** and loaded packages.

</div>

----------------------

## All users

### Install **R** and **RStudio**

**R** and **RStudio** are free and open-source.
You’ll need at least version 4.2.0 of **R**.
Source code and binaries are available on CRAN (https://cran.r-project.org/).
Install *the latest version* of **RStudio** from [rstudio.com](https://www.rstudio.com/products/rstudio/download/).

**For Windows users**, there is a <a href="https://youtu.be/n6mnN3lGj4s" target="_blank">video tutorial on the installation of **R** and **Rstudio**</a>.

### Install needed packages

Open a session in **RStudio** and run the following:

```
> update.packages()
> source("https://kingaa.github.io/sbied/prep/packages.R")
```

*[The `>` is the command prompt; it is not part of the command.
Also, the quotation marks are plain keyboard double quotes.]*

The first command updates your installed packages.
You may be prompted to specify a CRAN mirror:
choose one geographically near you.
In **RStudio**, you can also push the "Update" button on the "Packages" tab to accomplish this.

The second command runs a script on my website.
It will install some needed packages if these are not already installed on your system.

-------------------------------

## Windows users

If your machine runs Windows, you must install **Rtools**.
This will give you the ability to compile C code and dynamically link it into an **R** session.

[Download **Rtools** from CRAN](https://cran.r-project.org/bin/windows/Rtools/) and install it.
A <a href="https://youtu.be/lmIhiT_QsPE" target="_blank">video tutorial demonstrating how to install **Rtools** is available</a>, as are [detailed installation instructions](https://cran.r-project.org/bin/windows/Rtools/).
*Note that, after installation, there is one more step to be completed.*

***It is essential that you install these tools before the course starts!***


-------------------------------

## MacOS users

So that you can compile C code and dynamically link it into an **R** session, you will need to have the **Xcode** app installed.
This is gratis and can be installed via the App Store or downloaded from [developer.apple.com](https://developer.apple.com/download/).
Video tutorials demonstrating how to [check if you need to install Xcode](https://youtu.be/j0kqPpMecp4),
[how to install Xcode](https://youtu.be/aryEseip6Mk), and
[how to install pomp once you have installed Xcode](https://youtu.be/ikvJcN3Zi2w) are available.

Note that you must go beyond merely installing the **Xcode** app.
After you've installed the app, open a unix terminal (listed as the **Terminal** app under "Utilities" in the Finder) and run the following line
```
xcode-select --install
```
This will install the "Command Line Tools" that are needed to compile native C codes. Running this command is tantamount to the part of the Xcode installation tutorial video where we complete the installation using the graphical user interface.

------------------------------------

## All users

### Test **pomp**

Open a session in **RStudio** and run the following:
```
> source("https://kingaa.github.io/scripts/pompTest.R")
```
This will check whether you can work with **pomp**.

If it fails, try the following:
```
> source("https://kingaa.github.io/scripts/helloC.R",echo=TRUE)
```
If this fails to give the "Hello!" message, you will need to follow the instructions below corresponding to your OS before re-trying the `pompTest.R` script.

---------------------------

## Troubleshooting

### Linux and unix users

If you have trouble with any of the scripts above, make sure you have the GNU compiler collection (GCC), including **gfortran**, installed on your computer.
Linux distributions typically include this by default but it is not impossible that you have somehow avoided this.

### MacOS users

If the `pompTest.R` script fails because you cannot load **pomp**, try installing it from source.
The easiest way to do this is to execute
```
install.packages("pomp",repos="https://kingaa.github.io/")
```
in an **R** session.
You can also use the **devtools** package.
Do
```
install.packages("devtools")
library(devtools)
install_github("kingaa/pomp")
```
If, while trying to install from source, you receive the error,
```
make: gfortran-8.2: No such file or directory
```
or one that otherwise refers to `gfortran`, then it is likely that you do not have the necessary version of **gfortran** installed.
Have a look at [these instructions](https://mac.r-project.org/tools/) and contact me at the address above if these don’t work for you.

Some users have reported receiving an error complaining that 
```
'stdlib.h' file not found
```
This indicates that the command-line tools from **Xcode** have not been properly installed.
In a unix terminal, run 
```
xcode-select --install
```
to correct this problem.

### Windows users

You have probably failed to install the **Rtools** correctly.
Revisit the [instructions above](#windows-users).

---------------------------

## Once you've finished...

...please fill out [this online form](https://forms.gle/xhCQ2mGWZoVcm7n18) to help us prepare.

------------------------------

[**Pre-course instructions for Windows users (Video)**](https://www.youtube.com/playlist?list=PLluGwj6FGt2R3iM5CAEfIIof5dQfHVRgz)  
[**Pre-course instructions for MacOS users (Video)**](https://www.youtube.com/playlist?list=PLluGwj6FGt2S8GmrOF3s68LlsWoTRmuMQ)  
