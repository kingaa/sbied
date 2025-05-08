## Simulation-based Inference for Epidemiological Dynamics<br>Course Materials

[![DOI](https://zenodo.org/badge/39361090.svg)](https://doi.org/10.5281/zenodo.15364411)

These are notes from a short course given by [Ed
Ionides](https://ionides.github.io/) and [Aaron
King](https://kinglab.eeb.lsa.umich.edu/king/) at the [Summer Institute
in Statistics and Modeling in Infectious Diseases
(SISMID)](https://sismid.uw.edu). From the [main
page](https://kingaa.github.io/sbied/), links lead to pages on a number
of specific topics, culminating in four case studies that exemplify the
methods and raise key issues. For each such page, there is a
corresponding **R** script, which contains the codes needed to
recapitulate the calculations. This script is a starting point for
students to follow, explore, and modify the analysis according to their
own curiosity and interest.

------------------------------------------------------------------------

### Required software

The codes require, at a minimum, [**R**](https://cran.r-project.org/)
version 4.2 and [**pomp**](https://kingaa.github.io/pomp/) version 5.2.
Windows users must also have the appropriate version of
[**Rtools**](https://cran.r-project.org/bin/windows/Rtools/) installed.
The `prep` directory contains scripts that will install other needed
packages and test the user’s installation.

The following **R** packages are needed to build the materials, but are
not used *in* the materials:

-   **revealjs**
-   **diagram**
-   **mapproj**

The `contacts` lesson requires
[**panelPomp**](https://github.com/cbreto/panelPomp). An installation
script is provided there.

A [Github Action](https://github.com/kingaa/sbied/actions) checks that
these installations and actions succeed on a variety of current and
legacy platforms:

[![install-test](https://github.com/kingaa/sbied/actions/workflows/install-test.yml/badge.svg)](https://github.com/kingaa/sbied/actions/workflows/install-test.yml)

------------------------------------------------------------------------

### Manifest

-   index.md, index.html: Landing page
-   acknowledge.md, acknowledge.html: Acknowledgements
-   subdirectories:
    -   prep: Lesson 0, on preparing for the course
    -   intro: Lesson 1
    -   stochsim: Lesson 2
    -   pfilter: Lesson 3
    -   mif: Lesson 4
    -   measles: Lesson 5
    -   polio: Lesson 6
    -   ebola: Lesson 7
    -   contacts: Lesson 8
    -   od: Lesson from an earlier version of the course
    -   misc: Miscellaneous lessons
    -   \_includes: latex, HTML, **R** files included in other files
    -   graphics: figures in various formats
-   LICENSE.md: Full text of license
-   TODO.md: To-do list
-   Makefile, rules.mk: GNU Make files
-   sbied.bib: Bibliography database
-   jss.bst, jss.csl: Bibliography style files
-   CLUSTER.R: example `CLUSTER.R` file for site-specific cluster setup

------------------------------------------------------------------------

### Compilation

Full compilation of all the materials can be accomplished by running
`make` in the root directory. This requires substantial resources. For
this reason, the most expensive computations are archived using the
facilities provided for the purpose in **pomp**. Compilation with these
archives in place requires much less time than does compilation from
scratch. The following gives an indication of the size of the archives
and the time required for their computation on a linux cluster with
250 cpus.

<table>
<thead>
<tr class="header">
<th style="text-align: left;">directory</th>
<th style="text-align: right;">size (kB)</th>
<th style="text-align: right;">time (min)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">measles/results</td>
<td style="text-align: right;">340</td>
<td style="text-align: right;">224.0</td>
</tr>
<tr class="even">
<td style="text-align: left;">mif/results</td>
<td style="text-align: right;">2244</td>
<td style="text-align: right;">67.6</td>
</tr>
<tr class="odd">
<td style="text-align: left;">polio/results</td>
<td style="text-align: right;">41687</td>
<td style="text-align: right;">65.6</td>
</tr>
<tr class="even">
<td style="text-align: left;">polio/initial-values-exercise/results</td>
<td style="text-align: right;">41635</td>
<td style="text-align: right;">64.3</td>
</tr>
<tr class="odd">
<td
style="text-align: left;">polio/starting-values-exercise/results</td>
<td style="text-align: right;">41687</td>
<td style="text-align: right;">64.2</td>
</tr>
<tr class="even">
<td style="text-align: left;">contacts/results</td>
<td style="text-align: right;">206</td>
<td style="text-align: right;">13.7</td>
</tr>
<tr class="odd">
<td style="text-align: left;">pfilter/results</td>
<td style="text-align: right;">142</td>
<td style="text-align: right;">13.2</td>
</tr>
<tr class="even">
<td style="text-align: left;">polio/runlevel2/results</td>
<td style="text-align: right;">3266</td>
<td style="text-align: right;">4.8</td>
</tr>
<tr class="odd">
<td style="text-align: left;">polio/runlevel1/results</td>
<td style="text-align: right;">612</td>
<td style="text-align: right;">0.8</td>
</tr>
<tr class="even">
<td style="text-align: left;">ebola/results</td>
<td style="text-align: right;">2259</td>
<td style="text-align: right;">0.4</td>
</tr>
</tbody>
</table>

The archives amount to about 131 MB. About 126 MB of this is associated
with the polio lesson, in which large amounts of redundant information
are stored.

Full compilation, i.e., rebuilding the complete set of materials
following deletion of all archives, requires about 534 min on a
250-processor cluster. Full compilation regenerates the complete set of
archives. A finishing compilation, i.e., rebuilding with all archives in
place, but with the re-Making of all documents, requires about 12 min on
a 64-cpu workstation.

------------------------------------------------------------------------

### License

<figure>
<img src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png"
alt="CC-BY-NC" />
<figcaption aria-hidden="true">CC-BY-NC</figcaption>
</figure>

This work is licensed under the [Creative Commons
Attribution-NonCommercial 4.0 International License (CC BY-NC
4.0)](https://creativecommons.org/licenses/by-nc/4.0/). Under its terms,
you are free to:

-   Share — copy and redistribute the material in any medium or format
-   Adapt — remix, transform, and build upon the material

under the following terms:

-   Attribution — You must give appropriate credit, provide a link to
    the license, and indicate if changes were made. You may do so in any
    reasonable manner, but not in any way that suggests the licensor
    endorses you or your use.
-   NonCommercial — You may not use the material for commercial
    purposes.
-   No additional restrictions — You may not apply legal terms or
    technological measures that legally restrict others from doing
    anything the license permits.

The licensor cannot revoke these freedoms as long as you follow the
license terms.

------------------------------------------------------------------------
