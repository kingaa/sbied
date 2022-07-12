---
title: " Strategies for debugging **pomp** code"
author: "Aaron A. King and Edward L. Ionides"
output:
  html_document:
    toc: yes
    toc_depth: 2
bibliography: ../sbied.bib
csl: ../jss.csl
---


We review universal debugging considerations, and make suggestions specific to the structure of pomp and the workflow adopted in these lecture notes.

1. **Run the code from the top, in a clean R session, and identify the first error**. Rstudio facilitates testing code chunks in the middle of a source file, but if problems are not immediately solved it is better to start with a clean slate.

2. **Is the error (a) a standard R error; (b) a problem finding compilers needed by pomp; (c) a compilation error?** If you suspect (b), or want to rule it out, go back and re-run the pomp installation tests. If you have edited C snippets then (c) is likely. You will receive the C compiler error message in R.

3. **Does the code take annoyingly long to run, especially if running from the start as recommended above?** For effective debugging, your code should run in seconds. This may necessitate using very few particles (say, 10) and replications (say, 2). The polio case study introduces a flag system which can adjust the computational intensity simultaneously throughout the code by setting a single variable called `run_level`. Earlier lessons avoid this for simplicity of code. That means you have to edit each and every algorithmic parameter separately for convenient debugging.

4. **Cache errors**. bake and stew save results to avoid unnecessary recomputation while editing. This can cause problems for debugging. Remove the resulting rda or rds files to ensure proper recomputation.

