---
title: Pro Tips
author: Aaron A. King and Edward L. Ionides
output:
  html_document:
    toc: yes
    toc_depth: 2
    toc_float:
      collapsed: TRUE
      smooth_scroll: TRUE
    highlight: haddock
    number_sections: FALSE
    includes:
      after_body:
      - ../_includes/supp_bottom.html
      - ../_includes/license.html
bibliography: ../sbied.bib
csl: ../jss.csl
---


## When using **Rstudio**, turn off automatic save and restore of global workspace.

By default, when you quit a session, **R** asks whether to save the global session to a hidden file, `.RData`, in the current working directory.
Presumably because this behavior is annoying, **Rstudio**, by default, answers the "Do you want to save to `.RData`?" question for you, in the affirmative.
Also by default, and with only an easily overlooked message in the startup banner, **R** loads such a file, re-establishing the workspace at the start of a session.
Because the file is hidden, and the behavior is easy to forget about, this can lead to errors that are difficult to track down.
[For example, situations where different results are obtained on different machines during a large parallel computation, despite all the code being precisely the same!]
For these reasons, it's best to put a stop to all of this skulduggery.

To do so, go to the "Tools" menu in **Rstudio** and select "Global Options".
Make sure the "Restore .RData into workspace at startup" box is unticked.
For good measure, set the "Save workspace to .RData on exit" to "Never".

If you ever do want to save your workspace, it's as easy as `save.image(file="<filename>.rda")`;
restoring the file is a matter of `load("<filename>.rda")`.
When you do this, the file you create will be visible, as of course it should be since you gain nothing by hiding things from yourself!


## Never use whitespace in file or directory names

A perpetual source of annoying and easily avoidable errors!
See, e.g., [**pomp** issue 115](https://github.com/kingaa/pomp/issues/115).

## Keeping a database of parameter-space explorations

Likelihood surfaces for dynamic models can be very complex and the computations needed to explore them can be expensive.
By keeping a record of all parameter points visited, along with the computed likelihood at each point, is a good way to ensure that you continually improve your picture of the likelihood surface.

Doing this can be as simple as maintaining a CSV file with one column for each parameter, plus the likelihood (and s.e.).
It can be useful to supplement this with an indication of the name of the model and any other qualifying information.

The [Lesson on Iterated Filtering](../mif/index.html) shows one way of setting up, maintaining, and using such a database.
