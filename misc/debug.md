---
title: Strategies for debugging **pomp** code
author: Aaron A. King and Edward L. Ionides
output:
  html_document:
    toc: yes
    toc_depth: 2
    includes:
      after_body:
      - ../_includes/supp_bottom.html
      - ../_includes/license.html
bibliography: ../sbied.bib
csl: ../jss.csl
---


We review universal debugging considerations, and make suggestions specific to the structure of pomp and the workflow adopted in these lecture notes.

1. **Run the code from the top, in a clean R session, and identify the first error**.
Rstudio facilitates testing code chunks in the middle of a source file, but if problems are not immediately solved it is better to start with a clean slate.

1. **Read the error message.**
Error messages can be difficult to parse, but important information can be hidden in the messages.
Beginning at the top, look for statements of the form `Error: <error message text>`.
Also, do you see any references to familiar-looking code?
This can give a hint as to where in your code the problem arises.


2. **Is the error (a) a standard R error; (b) a problem finding compilers needed by pomp; (c) a compilation error?**
If you suspect (b), or want to rule it out, go back and re-run the pomp installation tests.
If you have edited C snippets then (c) is possible: 
the C compiler error message will be included **R**'s error message.

3. **Does the code take annoyingly long to run, especially if running from the start as recommended above?**
For effective debugging, your code should run in seconds.
This may necessitate using very few particles (say, 10) and replications (say, 2).
The [polio case study](../polio/) introduces a flag system which can adjust the computational intensity simultaneously throughout the code by setting a single variable called `run_level`.
Earlier lessons avoid this to keep the exposition simple.
That means you have to edit each and every algorithmic parameter separately for convenient debugging.

4. **Cache errors**. bake and stew save results to avoid unnecessary recomputation while editing.
This can cause problems for debugging.
Remove the `.rda` or `.rds` archive files to ensure proper recomputation.
