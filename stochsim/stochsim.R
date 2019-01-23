#' ---
#' title: "Simulation of stochastic dynamic models"
#' author: "Edward L. Ionides and Aaron A. King"
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 4
#' df_print: paged
#' bibliography: ../sbied.bib
#' csl: ../ecology.csl
#' nocite: >
#'   @Keeling2007
#' 
#' ---
#' 
#' \newcommand\prob[1]{\mathbb{P}\left[{#1}\right]}
#' \newcommand\expect[1]{\mathbb{E}\left[{#1}\right]}
#' \newcommand\var[1]{\mathrm{Var}\left[{#1}\right]}
#' \newcommand\dist[2]{\mathrm{#1}\left(#2\right)}
#' \newcommand\dlta[1]{{\Delta}{#1}}
#' \newcommand\lik{\mathcal{L}}
#' \newcommand\loglik{\ell}
#' \newcommand\Rzero{\mathfrak{R}_0}
#' 
#' --------------
#' 
#' 
#' [Licensed under the Creative Commons Attribution-NonCommercial license](http://creativecommons.org/licenses/by-nc/4.0/).
#' Please share and remix noncommercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
#' Produced with **R** version `r getRversion()` and **pomp2** version `r packageVersion("pomp2")`.
#' 
#' -------------
#' 
#' --------------
#' 
#' 
## ----prelims,echo=F,cache=F----------------------------------------------
library(plyr)
library(tidyverse)
library(pomp2)
theme_set(theme_bw())
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp2")>"2.0.9")
set.seed(594709947L)

#' 
#' ## Objectives
#' 
#' This tutorial develops some classes of dynamic models of relevance in biological systems, especially epidemiology.
#' We have the following goals:
#' 
#' 1. Dynamic systems can often be represented in terms of _flows_ between _compartments_.
#' We will develop the concept of a _compartmental model_ for which we specify _rates_ for the flows between compartments.
#' 1. We show how deterministic and stochastic versions of a compartmental model are derived and related.
#' 1. We introduce Euler's method to simulate from dynamic models, and we apply it to both deterministic and stochastic compartmental models.
#' 
#' <br>
#' 
#' ------
#' 
#' -----
#' 
#' ## Introduction
#' 
#' * Compartmental models are of great utility in many disciplines and very much so in epidemiology.
#' 
#' * We start by deriving deterministic and stochastic versions of the susceptible-infected-recovered (SIR) model of disease transmission dynamics in a closed population.
#' 
#' * In so doing, we will use notation that generalizes to more complex systems [[@breto09]](http://dx.doi.org/10.1214/08-AOAS201).
#' 
#' 
#' 
#' - Let $S(t)$, $I(t)$, and $R(t)$ represent, respectively, the number of susceptible hosts, the number of infected (and, by assumption, infectious) hosts, and the number of recovered or removed hosts. This corresponds to placing individuals into one of the compartments S, I and R in the diagram.
#' 
#' - We suppose that each arrow has an associated *per capita* rate, so here there is a rate $\mu_{SI}$ at which each individual in S transitions to I, and $\mu_{IR}$ at which each individual in I transitions to R. 
#' 
#' - To account for demography (birth/death/migration) we allow the possibility of a source and sink compartment, which is not represented on the flow diagram above.
#' - We write $\mu_{{\small\bullet} S}$ for a rate of births into S.
#' - Mortality rates are denoted by $\mu_{S{\small\bullet}}$, $\mu_{I{\small\bullet}}$, $\mu_{R{\small\bullet}}$.
#' 
#' - The rates may be either constant or varying. In particular, for a simple SIR model, the recovery rate $\mu_{IR}$ is a constant but the infection rate has the time-varying form $$\mu_{SI}(t)=\beta\,\frac{I(t)}{N(t)},$$ with $\beta$ being the _contact rate_ and $N$ the total size of the host population.
#' In the present case, since the population is closed, we set 
#' $$\mu_{{\small\bullet} S}=\mu_{S{\small\bullet}}=\mu_{I{\small\bullet}}=\mu_{R{\small\bullet}}=0.$$
#' 
#' - In general, it turns out to be convenient to keep track of the flows between compartments as well as the number of individuals in each compartment.
#' 
#' * Let $N_{SI}(t)$ count the number of individuals who have transitioned from $S$ to $I$ by time $t$. 
#' 
#' * We say that $N_{SI}(t)$ is a _counting process_.
#' 
#' * A similarly constructed process $N_{IR}(t)$ counts individuals transitioning from I to R.
#' 
#' * To include demography, we could keep track of birth and death events by the counting processes $N_{{\small\bullet} S}(t)$, $N_{S{\small\bullet}}(t)$, $N_{I{\small\bullet}}(t)$, $N_{R{\small\bullet}}(t)$.
#' - For discrete population compartment models, the flow counting processes are non-decreasing and integer valued.
#' - For continuous population compartment models, the flow counting processes are non-decreasing and real valued.
#' 
#' * The number of hosts in each compartment can be computed via these counting processes.
#' 
#' Ignoring demography, we have: 
#' $$\begin{aligned} 
#' S(t) &= S(0) - N_{SI}(t)\\
#' I(t) &= I(0) + N_{SI}(t) - N_{IR}(t)\\
#' R(t) &= R(0) + N_{IR}(t)\\
#' \end{aligned}$$
#' 
#' * These equations represent a kind of conservation law. Over any finite time interval $[t,t+\delta)$, we have
#' $$\begin{aligned} 
#' \dlta{S} &= -\dlta{N}_{SI}\\
#' \dlta{I} &= \phantom{-}\dlta{N}_{SI}-\dlta{N}_{IR}\\
#' \dlta{R} &= \phantom{-}\dlta{N}_{IR},
#' \end{aligned}$$
#' where the $\Delta$ notation indicates the increment in the corresponding process.
#' 
#' - Thus, for example $\dlta{N}_{SI}(t) = N_{SI}(t+\delta)-N_{SI}(t)$.
#' 
#' <br>
#' 
#' -------
#' 
#' ------
#' 
#' ## Compartmental models in theory
#' 
#' ### The deterministic version of the SIR model
#' 
#' Together with initial conditions specifying $S(0)$, $I(0)$ and $R(0)$, we just need to write down ordinary differential equations (ODE) for the flow counting processes.
#' These are,
#' $$\begin{gathered}
#' \frac{dN_{SI}}{dt} = \mu_{SI}(t)\,S(t), \qquad
#' \frac{dN_{IR}}{dt} = \mu_{IR}\,I(t).
#' \end{gathered}$$
#' 
#' <br>
#' 
#' -------
#' 
#' ------
#' 
#' ### The simple continuous-time Markov chain version of the SIR model
#' 
#' - Continuous-time Markov chains are the basic tool for building discrete population epidemic models.
#' 
#' - Recall that a _Markov chain_ is a discrete-valued stochastic process with the _Markov property_:
#' the future evolution of the process depends only on the current state.
#' 
#' - Surprisingly many models have this Markov property. 
#' If all important variables are included in the state of the system, then the Markov property appears automatically.
#' 
#' - The Markov property lets us specify a model by giving the transition probabilities on small intervals together with initial conditions.
#' 
#' - For the SIR model in a closed population, we have
#' $$\begin{aligned}
#' &\prob{N_{SI}(t+\delta)=N_{SI}(t)+1} &=& &\mu_{SI}(t)\,S(t)\,\delta + o(\delta)\\
#' &\prob{N_{SI}(t+\delta)=N_{SI}(t)} &=& &1-\mu_{SI}(t)\,S(t)\,\delta + o(\delta)\\
#' &\prob{N_{IR}(t+\delta)=N_{IR}(t)+1} &=& &\mu_{IR}(t)\,I(t)\,\delta + o(\delta)\\
#' &\prob{N_{IR}(t+\delta)=N_{IR}(t)} &=& &1-\mu_{IR}(t)\,I(t)\,\delta + o(\delta)\\
#' \end{aligned}$$
#' 
#' - A *simple* counting process is one for which no more than one event can occur at a time ([Wikipedia: point process](https://en.wikipedia.org/wiki/Point_process)).
#' 
#' - Thus, in a technical sense, the SIR Markov chain model we have written is simple.
#' 
#' - One may want to model the extra randomness resulting from multiple simultaneous events:
#' someone sneezing in a crowded bus, large gatherings at football matches, etc. 
#' This extra randomness may even be critical to match the variability in data. 
#' 
#' - We will see later, in the [measles case study](../measles/measles.html), a situation where this extra randomness plays an important role. 
#' The representation of the model in terms of counting processes turns out to be useful for this.
#' 
#' --------------------------
#' 
#' #### Optional Exercise: From Markov chain to ODE
#' 
#' Find the expected value of $N_{SI}(t+\delta)-N_{SI}(t)$ and $N_{IR}(t+\delta)-N_{IR}(t)$ given the current state, $S(t)$, $I(t)$ and $R(t)$.
#' Take the limit as $\delta\to 0$ and show that this gives the ODE model.
#' 
#' [Worked solution to the Optional Exercise](./exer_ctmc_ode.html#optional-exercise-from-markov-chain-to-ode)
#' 
#' --------------------------
#' 
#' ### Euler's method for an ODE
#' 
#' [Euler](https://en.wikipedia.org/wiki/Leonhard_Euler) took the following approach to numeric solution of an ODE:
#' 
#' - He wanted to investigate an ODE $\frac{dx}{dt}=h(x,t)$ with an initial condition $x(t_0)=x_0$.
#' 
#' - Since, in many cases of interest, the true solution $x(t)$ of this ODE cannot be worked out analytically, Euler wished to find a numerical approximation, $\tilde x(t)$, to $x(t)$.
#' 
#' - He divided time into small intervals of length $\delta$, with $\tilde t_k=t_0+k\,\delta$.
#' 
#' - He initialized the numerical approximation at the known starting value, $$\tilde x(t_0)=x(t_0)=x_0.$$
#' 
#' - He then reasoned that the gradient, $dx/dt$, is roughly constant over each small time interval $\tilde t_k\le t\le \tilde t_{k+1}$, for $k=1,2,\dots$. 
#' 
#' - Accordingly, he defined $$\tilde{x}(\tilde t_{k+1}) = \tilde{x}(\tilde t_k) + \delta\,h\big(\tilde{x}(\tilde t_{k}),\tilde t_{k}).$$
#' 
#' - This defines $\tilde x(t)$ only for $t=\tilde t_{k}$, but let's suppose $\tilde x(t)$ is constant between these discrete times (see the figure below).
#' 
#' - We now have a numerical scheme stepping forward in time increments of size $\delta$ which can be readily evaluated by computer.
#' 
#' - [Mathematical analysis of Euler's method](https://en.wikipedia.org/wiki/Euler_method) says that, as long as the function $h(x)$ is not too exotic, then $x(t)$ is well approximated by $\tilde x(t)$  when the discretization time-step, $\delta$, is sufficiently small.
#' 
#' - Euler's method is not the only numerical scheme to solve ODEs. More advanced schemes have better convergence properties, meaning that the numerical approximation is closer to $x(t)$.
#' 
#' - There are 3 reasons we choose to lean heavily on Euler's method:
#' 1. Euler's method is the simplest (the KISS principle).
#' 2. Euler's method extends naturally to stochastic models, both continuous-time Markov chains models and stochastic differential equation (SDE) models.
#' 3. In the context of data analysis, close approximation of the numerical solutions to a continuous-time model is less important than may be supposed, a topic discussed later.
#' 
#' --------------------------
#' 
#' 
#' <br>
#' 
#' --------
#' 
#' --------
#' 
#' 
#' ### Some comments on using continuous-time models and discretized approximations
#' 
#' - In some physical situations, a system follows an ODE model closely. 
#' For example, Newton's laws provide a very good approximation to the motions of celestial bodies.
#' 
#' - In many biological situations, ODE models become good approximations to reality only at relatively large scales. 
#' On small temporal scales, models cannot usually capture the full scope of biological variation and biological complexity.
#' 
#' - If we are going to expect substantial error in using $x(t)$ to model a biological system, maybe the numerical solution $\tilde x(t)$ represents the system being modeled as well as $x(t)$  does.
#' 
#' - If our model fitting, model investigation, and final conclusions are all based on our numerical solution  $\tilde x(t)$ (e.g., we are sticking entirely to simulation-based methods) then we are most immediately concerned with how well $\tilde x(t)$ describes the system of interest.  
#' $\tilde x(t)$ becomes more important than the original model, $x(t)$.
#' 
#' - When following this perspective, it is important that one fully describe the numerical model $\tilde x(t)$.
#' 
#' -  From this point of view, then, the main advantage of the continuous-time model $x(t)$ is then that it gives a succinct way to describe how $\tilde x(t)$ was constructed.
#' 
#' - All numerical methods are, ultimately, discretizations.
#' + For continuous-time modeling, we usually aim to set $\delta$ small compared to the timescale of the process being modeled, so that the choice of $\delta$ does not play an explicit role in the interpretation of the model.
#' + Epidemiologically, setting $\delta$ to be a day, or an hour, can be quite different from setting $\delta$ to be two weeks or a month. 
#' + Putting more emphasis on the scientific role of the numerical solution itself reminds you that the numerical solution has to do more than approximate a target model in some asymptotic sense: 
#' the numerical solution should be a sensible model in its own right. 
#' 
#' <br>
#' 
#' -----
#' 
#' -----
#' 
#' ### Euler's method for a discrete SIR model
#' 
#' - Recall the simple continuous-time Markov chain interpretation of the SIR model without demography:
#' $$\begin{aligned}
#' \prob{N_{SI}(t+\delta)=N_{SI}(t)+1} &= \mu_{SI}(t)\,S(t)\,\delta + o(\delta),\\
#' \prob{N_{IR}(t+\delta)=N_{IR}(t)+1} &= \mu_{IR}\,I(t)\,\delta + o(\delta).
#' \end{aligned}$$
#' 
#' - We look for a numerical solution with state variables $\tilde S(t_k)$, $\tilde I(t_k)$, $\tilde R(t_k)$. 
#' 
#' - The counting processes for the flows between compartments are $\tilde N_{SI}(t)$ and $\tilde N_{IR}(t)$. The counting processes are related to the numbers of individuals in the compartments by the same flow equations we had before:
#' $$\begin{aligned} 
#' \dlta{\tilde S} &= -\dlta{\tilde N}_{SI}\\
#' \dlta{\tilde I} &= \dlta{\tilde N}_{SI}-\dlta{\tilde N}_{IR}\\
#' \dlta{\tilde R} &= \dlta{\tilde N}_{IR},
#' \end{aligned}$$
#' 
#' - Let's focus $N_{SI}(t)$;
#' the same methods can also be applied to $N_{IR}(t)$.
#' 
#' - Here are three stochastic Euler schemes for $N_{SI}$:
#' 1. Poisson increments:
#' $$\dlta{\tilde N}_{SI}\;\sim\;\dist{Poisson}{\tilde \mu_{SI}(t)\,\tilde S(t)\,\delta},$$ where $\dist{Poisson}{\mu}$ is the Poisson distribution with mean $\mu$ and $$\tilde\mu_{SI}(t)=\beta\,\frac{\tilde I(t)}{N}.$$
#' 1. Binomial increments with linear probability:
#' $$\dlta{\tilde N}_{SI}\;\sim\;\dist{Binomial}{\tilde{S}(t),\tilde\mu_{SI}(t)\,\delta},$$ where $\dist{Binomial}{n,p}$ is the binomial distribution with mean $n\,p$ and variance $n\,p\,(1-p)$.
#' 1. Binomial increments with exponential decaying probability:
#' $$\dlta{\tilde{N}}_{SI}\;\sim\;\dist{Binomial}{\tilde{S}(t),1-e^{-\tilde{\mu}_{SI}(t)\,\delta}}.$$
#' 
#' - Note that these schemes all become equivalent as $\delta\to 0$.
#' 
#' - What are the advantages and disadvantages of these different schemes?
#' Conceptually, it is simplest to think of (1) or (2). 
#' Numerically, it is usually preferable to implement (3). 
#' 
#' <br>
#' 
#' ------
#' 
#' ------
#' 
#' ### Compartmental models via stochastic differential equations (SDE)
#' 
#' The Euler method extends naturally to stochastic differential equations. A natural way to add stochastic variation to an ODE $dx/dt=h(x)$ is
#' $$\frac{dX}{dt}=h(X)+\sigma\,\frac{dB}{dt}$$
#' where $B(t)$ is Brownian motion and so $dB/dt$ is Gaussian white noise.
#' The so-called Euler-Maruyama approximation $\tilde X$ is generated by 
#' $$\tilde X(\tilde t_{k+1}) = \tilde X(\tilde t_k) + \delta\,h\big(\tilde X(\tilde t_k)\big) + \sigma \sqrt{\delta} \, Z_k$$
#' where $Z_1,Z_2,\dots$ is a sequence of independent standard normal random variables, i.e., $Z_k\sim\dist{Normal}{0,1}$. 
#' Although SDEs are often considered an advanced topic in probability, the Euler approximation doesn't demand much more than familiarity with the normal distribution.
#' 
#' <br>
#' 
#' --------
#' 
#' --------
#' 
#' #### Optional Exercise: SDE version of the SIR model
#' 
#' Write down the Euler-Maruyama method for an SDE representation of the closed-population SIR model. 
#' Consider some difficulties that might arise with non-negativity constraints, and propose some practical way one might deal with that issue.
#' 
#' 
#' A useful method to deal with positivity constraints is to use Gamma noise rather than Brownian noise [@bhadra11;@He2010;@laneri10].
#' SDEs driven by Gamma noise can be investigated by Euler solutions simply by replacing the Gaussian noise by an appropriate Gamma distribution.
#' 
#' [Worked solution to the Optional Exercise](./exer_ctmc_ode.html#optional-exercise-sde-version-of-the-sir-model)
#' 
#' <br>
#' 
#' ------
#' 
#' -----
#' 
#' ### Euler's method vs.&nbsp;Gillespie's algorithm
#' 
#' - A widely used, exact simulation method for continuous time Markov chains is [Gillespie's algorithm](https://en.wikipedia.org/wiki/Gillespie_algorithm) [@Gillespie1977a].
#' 
#' - We do not put much emphasis on Gillespie's algorithm here. The Euler method is applicable to a wider class of models.
#' 
#' - When would you prefer an implementation of Gillespie's algorithm to an Euler solution?
#' 
#' - Numerically, Gillespie's algorithm is often approximated using so-called [$\tau$-leaping](https://en.wikipedia.org/wiki/Tau-leaping) methods [@Gillespie2001], which are closely related to Euler's approach with $\delta$ being $\tau$.
#' 
#' - Indeed, an Euler solution for a continuous time Markov chain is sometimes called a Gillespie tau-leaping method.
#' 
#' <br>
#' 
#' --------
#' 
#' -------
#' 
#' ## Compartmental models in **pomp2**.
#' 
#' ### The boarding-school flu outbreak
#' 
#' * As an example that we can probe in some depth, let's look at an isolated outbreak of influenza that occurred in a boarding school for boys in England [@Anonymous1978].
#' 
#' <!--- 763 boys were at risk, and ultimately 512 spent time away from class (either confined to bed or in convalescence. --->
#' 
#' * The data are provided in the package.
#' Examine them:
## ----flu-data1-----------------------------------------------------------
library(plyr)
library(tidyverse)
library(pomp2)

bsflu

#' 
#' * The variable `B` refers to boys confined to bed and `C` to boys in convalescence.
#' 
#' * Let's restrict our attention for the moment to the `B` variable.
#' 
## ----flu-data2,echo=F----------------------------------------------------
bsflu %>%
  ggplot(aes(x=day,y=B))+
  geom_line()+
  geom_point()

#' 
#' <br>
#' 
#' -------
#' 
#' ------
#' 
#' ### A first POMP model
#' 
#' * We start by treating $B$ as incidence data, counting the number of boys diagnosed with flu and confined to bed the preceding day.
#' 
#' * This lets us model the disease using the simple SIR model.
#' 
#' * Our tasks will be, first, to estimate the parameters of the SIR and, second, to decide whether or not the SIR model is an adequate description of these data.
#' 
#' * Below is a diagram of the SIR model.
#' The host population is divided into three classes according to their infection status: 
#' S, susceptible hosts; 
#' I, infected (and infectious) hosts; 
#' R, recovered and immune hosts. 
#' 
#' * The rate at which individuals move from S to I is the __force of infection__, $\lambda=\mu_{SI}=\beta\,I/N$, while that at which individuals move into the R class is $\mu_{IR}=\gamma$.
#' 
#' 
#' * Let's look at how we can view the SIR as a POMP model.
#' 
#' * The unobserved state variables, in this case, are the numbers of individuals, $S(t)$, $I(t)$, $R(t)$ in the S, I, and R compartments, respectively.
#' 
#' * It's reasonable in this case to view the population size $N=S(t)+I(t)+R(t)$, as fixed.
#' 
#' * The numbers that actually move from one compartment to another over any particular time interval are modeled as stochastic processes.
#' 
#' * In this case, we'll assume that the stochasticity is purely demographic, i.e., that each individual in a compartment at any given time faces the same risk of exiting the compartment.
#' 
#' * __Demographic stochasticity__ is the unavoidable randomness when dynamics arise from chance events occuring in a discrete and finite population.
#' 
#' <br>
#' 
#' -----
#' 
#' -----
#' 
#' ### Implementing the model
#' 
#' * To implement the model in **pomp2**, the first thing we need is a stochastic simulator for the unobserved state process.
#' 
#' * We've seen that there are several ways of approximating the process just described for numerical purposes.
#' 
#' * An attractive option here is to model the number moving from one compartment to the next over a very short time interval as a binomial random variable.
#' 
#' * In particular, we model the number, $\dlta{N_{SI}}$, moving from S to I over interval $\dlta{t}$ as $$\dlta{N_{SI}} \sim \dist{Binomial}{S,1-e^{-\beta\,I/N\dlta{t}}},$$ and the number moving from I to R as $$\dlta{N_{IR}} \sim \dist{Binomial}{I,1-e^{-\gamma\dlta{t}}}.$$
#' 
## ----rproc1R-------------------------------------------------------------
sir_step <- function (S, I, R, N, delta.t, Beta, gamma, ...) {
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-gamma*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  c(S = S, I = I, R = R)
}

#' 
#' * At day zero, we'll assume that $I=1$ and $R=0$, but we don't know how big the school is, so we treat $N$ as a parameter to be estimated and let $S(0)=N-1$.
#' 
## ----init1R--------------------------------------------------------------
sir_init <- function(N, ...) {
  c(S = N-1, I = 1, R = 0)
}

#' 
#' * We fold these basic model components, with the data, into a `pomp` object thus:
#' 
## ----pomp1R--------------------------------------------------------------
bsflu %>%
  pomp(times="day",t0=0,
    rprocess=euler(sir_step,delta.t=1/6),
    rinit=sir_init
  ) -> sir

#' 
#' * Now let's assume that the case reports, $B$, result from a process by which new infections result in confinement with probability $\rho$, which we can think of as the probability that an infection is severe enough to be noticed by the school authorities.
#' 
#' * Since confined cases have, presumably, a much lower transmission rate, let's treat $B$ as being a count of the number of boys who have moved from I to R over the course of the past day.
#' 
#' * We need a variable to track these daily counts.
#' Let's modify our rprocess function above, adding a variable $H$ to tally the incidence.
#' We'll then replace the rprocess component in `sir` with the new one.
#' In **pomp2** terminology, $H$ is an *accumulator variable*.
#' 
## ----rproc2R-------------------------------------------------------------
sir_step <- function (S, I, R, H, N, delta.t, Beta, gamma, ...) {
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-gamma*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  H <- H + dN_IR;
  c(S = S, I = I, R = R, H = H)
}

sir_init <- function (N, ...) {
  c(S = N-1, I = 1, R = 0, H = 0)
}

sir %>%
  pomp(rprocess=euler(sir_step,delta.t=1/6),rinit=sir_init) -> sir

#' 
#' * Now, we'll model the data, $B$, as a binomial process,
#' $$B_t \sim \dist{Binomial}{H(t)-H(t-1),\rho}.$$
#' 
#' * We have a problem: at time $t$, the variable `H` we've defined will contain $H(t)$, not $H(t)-H(t-1)$.
#' 
#' * We can overcome this by telling `pomp` that we want `H` to be set to zero immediately following each observation.
#' 
#' * We accomplish this using the `accumvars` argument to `pomp`:
#' 
## ----zero1R--------------------------------------------------------------
sir %>% pomp(accumvars="H") -> sir

#' 
#' * Now, to include the observations in the model, we must write both a `dmeasure` and an `rmeasure` component:
#' 
## ----meas-modelR---------------------------------------------------------
dmeas <- function (B, H, rho, log, ...) {
 dbinom(x=B, size=H, prob=rho, log=log)
}

rmeas <- function (H, rho, ...) {
  c(B=rbinom(n=1, size=H, prob=rho))
}

#' 
#' * We then put these into our `pomp` object:
#' 
## ----add-meas-modelR-----------------------------------------------------
sir %>% pomp(rmeasure=rmeas,dmeasure=dmeas) -> sir

#' 
## ----test_R_pomp,include=FALSE-------------------------------------------
sir %>% simulate(params=c(Beta=1.5,gamma=1,rho=0.9,N=2600))

#' 
#' #### Using C snippets
#' 
#' * Although we can always specify basic model components using **R** functions, as above, we'll typically want the computational speed up that we can obtain only by using compiled native code.
#' 
#' * **pomp2** provides a facility for doing so with ease, using *C snippets*.
#' 
#' * A C snippet is a small piece of C code used to specify a basic model component.
#' 
#' * For example, C snippets that encode the basic model components in `sir` are as follows.
#' 
## ----csnips--------------------------------------------------------------
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  H = 0;
")

dmeas <- Csnippet("
  lik = dbinom(B,H,rho,give_log);
")

rmeas <- Csnippet("
  B = rbinom(H,rho);
")

#' 
#' * A call to `pomp` replaces the basic model components with these, much faster, implementations:
#' 
## sir %>%

#' 
#' * Note that, when using C snippets, one has to tell `pomp` which of the variables referenced in the C snippets are state variables and which are parameters.
#' This is accomplished using the `statenames` and `paramnames` arguments.
#' 
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Testing the model: simulations
#' 
#' * Let's perform some simulations, just to verify that our codes are working as intended.
#' 
#' * To do so, we'll need some parameters. A little thought will get us some ballpark estimates.
#' 
#' * In the data, it looks like there were a total of `r sum(bsflu$B)` infections, so the population size, $N$, must be somewhat in excess of this number.
#' 
#' * We can use the final-size equation from the basic theory of SIR epidemics,
#' $$\Rzero = -\frac{\log{(1-f)}}{f},$$
#' where $f=R(\infty)/N$ is the final size of the epidemic and $\Rzero$ is the basic reproduction number giving the expected number of secondary infections from one primary infection in a fully susceptible population.
#' 
#' * Typically, $\Rzero$ for influenza is thought to be around 1.5, from which we can estimate that $f\approx 0.6$, whence $N\approx 2600$.
#' 
#' * If the infectious period is roughly 1&nbsp;da, then $1/\gamma \approx 1~\text{da}$ and $\beta = \gamma\,\Rzero \approx 1.5~\text{da}^{-1}$.
#' 
#' * Let's simulate the model at these parameters.
## ----sir_sim1------------------------------------------------------------
sir %>%
  simulate(params=c(Beta=1.5,gamma=1,rho=0.9,N=2600),
    nsim=20,format="data.frame",include.data=TRUE) -> sims

sims %>%
  ggplot(aes(x=day,y=B,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' <br>
#' 
#' ------------
#' 
#' -----------
#' 
#' #### Basic Exercise: Explore the SIR model
#' 
#' Fiddle with the parameters to see if you can't find parameters for which the data are a more plausible realization.
#' 
#' [Worked solution to the Exercise](./exercises.html#basic-exercise-explore-the-sir-model)
#' 
#' <br>
#' 
#' ----------
#' 
#' ----------
#' 
#' #### Basic Exercise: The SEIR model
#' 
#' Below is a diagram of the so-called SEIR model.
#' This differs from the SIR model in that infected individuals must pass a period of latency before becoming infectious.
#' 
#' 
#' Modify the codes above to construct a `pomp` object containing the flu data and an SEIR model.
#' Perform simulations as above and adjust parameters to get a sense of whether improvement is possible by including a latent period.
#' 
#' [Worked solution to the Exercise](./exercises.html#basic-exercise-the-seir-model)
#' 
#' <br>
#' 
#' ------------
#' 
#' ------------
#' 
#' #### Basic Exercise: Rethinking the boarding-school flu data
#' 
#' In the preceding, we've been assuming that $B_t$ represents the number of boys *sent* to bed on day $t$.
#' Actually, this isn't correct at all.
#' As described in the report [@Anonymous1978], $B_t$ represents the total number of boys *in bed* on day $t$, not the number of new infections on day $t$, i.e., the incidence, as we've been assuming.
#' Since boys were potentially confined for more than one day, the data count each infection multiple times.
#' On the other hand, we have information about the total number of boys at risk and the total number who were infected.
#' In fact, we know that $N=763$ boys were at risk and 512 boys in total spent between 3 and 7 days away from class (either in bed or convalescent).
#' Moreover, we have data on the number of boys, $C_t$, convalescent at day $t$.
#' Since $1540~\text{boy-da}/512~\text{boy} \approx 3~\text{da}$, we know that the average duration spent in bed was 3&nbsp;da and, since $\sum_t\!C_t=`r sum(bsflu$C)`$, we can infer that the average time spent convalescing was $`r sum(bsflu$C)`~\text{boy-da}/512~\text{boy} \approx `r signif(sum(bsflu$C)/512,2)`~\text{da}$.
#' 
#' 
#' Formulate a model with a latent class and both confinement and convalescent stages.
#' Implement it in **pomp2** using a compartmental model like that diagrammed below.
#' 
#' 
#' In general, there is interest in using all of the data to inform the model but to economize time, let's focus on the confinement data, $B$.
#' Give some thought as to how to model the relationship between $B$ and the state variables.
#' 
#' Although we might use the data to estimate all the model's parameters, to keep the exercise simple, we'll fix some of the parameters at reasonable values.
#' Which parameters can be estimated from general knowlege of influenza?
#' Fix these and focus on estimating the remainder.
#' 
#' [Partially worked solution to the Exercise](./seir3_exercise.html)
#' 
#' --------------------------
#' 
#' ## [Back to course homepage](../index.html)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/master/stochsim/stochsim.R)
#' 
#' --------------------------
#' 
#' ## References
#' 
