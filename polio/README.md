polio.Rmd : source file, run with 
```
Rscript -e 'rmarkdown::render("polio.Rmd")'
```
You can also extract the **R** code at the **R** prompt with
```
require(knitr)
purl("polio.Rmd")
```

polio.html : the resulting document

polio.R: the resulting **R** script

polio.bib : bibliography

polio_params.csv : results of various previous searches

polio_wisconsin.csv : data

polio_fig1A.png : Fig 1A from Martinez-Bakker et al (2015)
