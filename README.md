## Simulation-based Inference for Epidemiological Dynamics
### Course Materials

These are notes from a short course given by Ed Ionides and Aaron King at the [Summer Institute in Statistics and Modeling in Infectious Diseases (SISMID)](http://sismid.uw.edu).
From the [main page](http://kingaa.github.io/sbied/), links lead to pages on a number of specific topics, culminating in four case studies that exemplify the methods and raise key issues.
For each such page, there is a corresponding **R** script, which contains the codes needed to recapitulate the calculations.
This script is a starting point for students to follow, explore, and modify the analysis according to their own curiosity and interest.

----------------------------

### Manifest

- index.md, index.html: Landing page
- acknowledge.md, acknowledge.html: Acknowledgements
- syllabus.md, syllabus.html: Course syllabus
- welcome.md, welcome.html: Course introduction slides
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
- subdirectories:
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

### Compilation

Full compilation of all the materials can be accomplished by running `make` in the root directory.
This requires substantial resources.
For this reason, the most expensive computations are archived using the facilities provided for the purpose in **pomp**.
Compilation with these archives in place requires much less time than does compilation from scratch.

<!--
Get archive file sizes in kB:

```
library(tidyverse)

gettime <- function (path) {
  if (grepl(".rds",path)) {
    readRDS(path) |> attr("system.time") |> getElement(3)
  } else if (grepl(".rda",path)) {
    e <- new.env()
    load(path,envir=e)
    e$.system.time[3]
  } else {
    NA_real_
  }
}

list.files(pattern=r"{.+\.rd[as]$}",recursive=TRUE) |> 
  {\(x) data.frame(
          path=x,
          dir=dirname(x),
          file=basename(x),
          size=file.size(x)/2^10,
          time=sapply(x,gettime)/60
        )}() |>
  remove_rownames() -> dat
dat		

dat |>
  dplyr::group_by(dir) |>
  dplyr::summarize(size=sum(size),time=sum(time))

dat |> 
  dplyr::summarize(size=sum(size),time=sum(time))
```
-->

At last count, the archives amount to about 129MB.
About 123MB of this is six files in the polio lesson, in which large amounts of redundant information are stored.


The following is a log of compilations.

  - *Full* compilation refers to compilation following the deletion of all archives.
  Full compilation generates the complete set of archives.
  - *Finishing* compilation refers to compilation with all archives in place, but with the re-Making of all documents.
  
| Type          | Completion date | Time required | CPUs available |
|:--------------|:----------------|--------------:|---------------:|
| finishing     | 2021-07-20      |         11min |             36 |
| full          | 2021-07-21      |        488min |            250 |
| polio level 3 | 2021-07-21      |           min |            250 |


----------------------------

### License

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
