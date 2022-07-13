---
title: "Simulation-based Inference for Epidemiological Dynamics"
author: "Aaron A. King, Edward L. Ionides, Qianying Lin"
output:
  html_document:
    includes:
      after_body:
      - _includes/main_bottom.html
      - _includes/license.html
bibliography: sbied.bib
csl: jss.csl

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

----------------------

## Module description

This module introduces statistical inference techniques and computational methods for dynamic models of epidemiological systems.
The course will explore deterministic and stochastic formulations of epidemiological dynamics and develop inference methods appropriate for a range of models.
Special emphasis will be on exact and approximate likelihood as the key elements in parameter estimation, hypothesis testing, and model selection. Specifically, the course will cover sequential Monte Carlo and synthetic likelihood techniques.
Students will learn to implement these in R to carry out maximum likelihood and Bayesian inference. Knowledge of the material in Module 1 is assumed.
Students new to R should complete a [tutorial](https://kingaa.github.io/R_Tutorial/) before the module.

----------------------

## Course objectives

1. To introduce partially observed Markov process (POMP) models as tools for scientific investigation and public health policy.
1. To give students the ability to formulate POMP models of their own.
1. To teach efficient approaches for performing scientific inference using POMP models.
1. To familiarize students with the **pomp** package.
1. To give students opportunities to work with such inference methods.
1. To provide documented examples for student re-use.

----------------------

## Lessons

[**0. Instructions for preparing your laptop for the course exercises.**](./prep/)

[**1. Introduction: What is "Simulation-based Inference for Epidemiological Dynamics"?  POMPs and pomp.**](./intro/)

[**2. Simulation of stochastic dynamic models.**](./stochsim/)

[**3. Likelihood for POMPs: theory and practice.**](./pfilter/)

[**4. Iterated filtering: theory and practice.**](./mif/)

[**5. Case study: measles.  Recurrent epidemics, long time series, covariates, extra-demographic stochasticity, interpretation of parameter estimates.**](./measles/)

[**6. Case study: polio. Workflow for a real research problem.**](./polio/)

[**7. Case study: Ebola. Model diagnostics and forecasting.**](./ebola/)

[**8. Case study: HIV and fluctuating sexual contact rates. Panel data.**](./contacts/)

----------------------
