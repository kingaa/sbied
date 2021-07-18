### Simulation-based Inference for Epidemiological Dynamics
#### [Course Notes](http://kingaa.github.io/sbied/)

These are notes from a short course given by Ed Ionides and Aaron King at the [Summer Institute in Statistics and Modeling in Infectious Diseases (SISMID)](http://sismid.uw.edu).
From the [main page](http://kingaa.github.io/sbied/), links lead to pages on a number of specific topics, culminating in four case studies that exemplify the methods and raise key issues.
For each such page, there is a corresponding **R** script, which contains the codes needed to recapitulate the calculations.
This script is a starting point for students to follow, explore, and modify the analysis according to their own curiosity and interest.

----------------------------

#### Manifest

- index.md, index.html: Landing page
- welcome.md, welcome.html: Course introduction slides
- acknowledge.md, acknowledge.html: Acknowledgements
- syllabus.md, syllabus.html: Course syllabus
- EDTsyllabus.Rmd, EDTsyllabus.html: Course syllabus with timezone offset
- LICENSE.md: Full text of license
- TODO.md: To-do list
- Makefile, rules.mk: GNU Make files
- sbied.bib: Bibliography database
- jss.bst, jss.csl: Bibliography style files
- header.tex: General latex preamble commands
- beamer.tex: Latex beamer preamble commands
- logo.jpg: SBIED logo
- packages.R: code for loading needed packages
- CLUSTER-kinglab.R: example `CLUSTER.R` file for site-specific cluster setup
- setup.Rmd, setup.Rnw: **knitr** child documents for R-markdown and R-noweb documents, respectively
- prep: Lesson 0, on preparing for the course
- intro: Lesson 1
- stochsim: Lesson 2
- pfilter: Lesson 3
- mif: Lesson 4
- measles: Lesson 5
- polio: Lesson 6
- ebola: Lesson 7
- contacts: Lesson 8
- od: Lesson from an earlier version of the course
- quiz: Quiz materials

----------------------------

#### Compilation

Full compilation of all the materials can be accomplished by running `make` in the root directory&ast;.
This requires substantial resources.
For this reason, the most expensive computations are archived using the facilities provided for the purpose in **pomp**.
Compilation with these archives in place requires much less time than does compilation from scratch.

The following is a log of compilations.

  - *Full* compilation refers to compilation following the deletion of all archives.
  Full compilation generates the complete set of archives.
  - *Finishing* compilation refers to compilation with all archives in place, but with the re-Making of all documents.
  
| Type    | Completion date | Time required | CPUs available | Archive size |
|:--------|:----------------|--------------:|---------------:|-------------:|
| Full    | 2021-07-19      |        587min |            250 |         59MB |
| Partial | 2021-07-19      |         11min |             36 |         59MB |

&ast;Some of the folders within the `polio` lesson cannot be reproduced using `make`.
These include `initial-values-exercise` and `starting-values-exercise`.

----------------------------

#### License

![CC-BY-NC](https://i.creativecommons.org/l/by-nc/4.0/88x31.png)

This work is licensed under the [Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)](http://creativecommons.org/licenses/by-nc/4.0/).
Under its terms, you are free to:

- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material

The licensor cannot revoke these freedoms as long as you follow the license terms.
Under the following terms:

- Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
- NonCommercial — You may not use the material for commercial purposes.
- No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

----------------------------
