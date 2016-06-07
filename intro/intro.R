#' ---
#' title: "Introduction to Simulation-based Inference"
#' author: "Aaron A. King and Edward L. Ionides"
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 4
#' bibliography: ../sbied.bib
#' csl: ../ecology.csl
#' nocite: |
#'   @King2008, @Romero-Severson2015, @He2010, 
#'   @Laneri2010, @King2015
#' ---
#' 
#' \newcommand\prob[1]{\mathbb{P}\left[{#1}\right]}
#' \newcommand\expect[1]{\mathbb{E}\left[{#1}\right]}
#' \newcommand\var[1]{\mathrm{Var}\left[{#1}\right]}
#' \newcommand\dist[2]{\mathrm{#1}\left(#2\right)}
#' \newcommand\dlta[1]{{\Delta}{#1}}
#' \newcommand\lik{\mathcal{L}}
#' \newcommand\loglik{\ell}
#' 
#' [Licensed under the Creative Commons Attribution-NonCommercial license](http://creativecommons.org/licenses/by-nc/4.0/).
#' Please share and remix noncommercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
#' Produced in **R** version `r getRversion()` using **pomp** version `r packageVersion("pomp")`.
#' 
## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(stringsAsFactors=FALSE)
library(ggplot2)
theme_set(theme_bw())
set.seed(2028866059L)

#' 
#' ## Introduction: ecological and epidemiological dynamics
#' 
#' - Ecological systems are complex, open, nonlinear, and nonstationary.
#' - It is useful to model them as stochastic systems.
#' - "Laws of Nature" are unavailable except in the most general form.
#' - For any observable phenomenon, multiple competing explanations are possible.
#' - Central scientific goals:
#'     - Which explanations are most favored by the data?
#'     - Which kinds of data are most informative?
#' - Central applied goals:
#'     - How to design ecological or epidemiological intervention?
#'     - How to make accurate forecasts?
#' - Time series are particularly useful sources of data.
#' 
#' 
#' ### Noisy clockwork: Time series analysis of population fluctuations in animals
#' 
#' ##### Six problems of @Bjornstad2001
#' 
#' Obstacles for **ecological** modeling and inference via nonlinear mechanistic models:
#' 
#' 1. Combining measurement noise and process noise.
#' 2. Including covariates in mechanistically plausible ways.
#' 3. Using continuous-time models.
#' 4. Modeling and estimating interactions in coupled systems. 
#' 5. Dealing with unobserved variables.
#' 6. Modeling spatial-temporal dynamics.
#' 
#' The same issues arise for **epidemiological** modeling and inference via nonlinear mechanistic models.
#' 
#' ## Objectives
#' 
#' 1. To show how stochastic dynamical systems models can be used as scientific instruments.
#' 1. To give students the ability to formulate models of their own.
#' 1. To teach efficient approaches for performing scientific inference using POMP models.
#' 1. To familiarize students with the **pomp** package.
#' 1. To give students opportunities to work with such inference methods.
#' 1. To provide documented examples for student re-use.
#' 
#' 
#' ### Questions and answers
#' 
#' 1. [What roles are played by asymptomatic infection and waning immunity in cholera epidemics?](http://dx.doi.org/10.1038/nature07084)
#' 7. [Do subclinical infections of pertussis play an important epidemiological role?](http://dx.doi.org/10.1371/journal.pone.0072086)
#' 3. [What explains the seasonality of measles?](http://dx.doi.org/10.1098/rsif.2009.0151)
#' 2. [What is the contribution to the HIV epidemic of dynamic variation in sexual behavior of an individual over time? How does this compare to the role of heterogeneity between individuals?](http://dx.doi.org/10.1093/aje/kwv044)
#' 5. [What explains the interannual variability of malaria?](http://dx.doi.org/10.1371/journal.pcbi.1000898)
#' 6. [What will happen next in an Ebola outbreak?](http://dx.doi.org/10.1098/rspb.2015.0347)
#' 
#' -------------------------------
#' 
#' ## Partially observed Markov process (POMP) models
#' 
#' * Data $y^*_1,\dots,y^*_N$ collected at times $t_1<\dots<t_N$ are modeled as noisy and incomplete observations of a Markov process $\{X(t), t\ge t_0\}$.
#' * This is a __partially observed Markov process (POMP)__ model, also known as a hidden Markov model or a state space model.
#' * $\{X(t)\}$ is Markov if the history of the process, $\{X(s), s\le t\}$, is uninformative about the future of the process, $\{X(s), s\ge t\}$, given the current value of the process, $X(t)$. 
#' * If all quantities important for the dynamics of the system are placed in the __state__, $X(t)$, then the Markov property holds by construction.
#' * POMP models can include all the features desired by @Bjornstad2001.
#' 
#' ------------------------------
#' 
#' ### Schematic of the structure of a POMP
#' 
#' showing causal relations.
#'    
#' ![](pomp_schematic1.png)
#' 
#' 
#' **The key perspective to keep in mind is that the model is to be viewed as the process that generated the data.**
#' 
#' ------------------------------
#' 
#' #### The Markov assumption
#' 
#' - $\prob{X_n|X_0,\dots,X_{n-1}}=\prob{X_n|X_{n-1}}$.
#' - Interpretation: knowledge of the system's state at any point in time is sufficient to determine the distribution of possible futures.
#' - Alternative interpretation: the system's state is sufficiently rich so as to encompass all important features of the system's history
#' - Systems with delays can usually be rewritten as Markovian systems, at least approximately.
#' - An important special case: any system of differential equations is Markovian.
#' 
#' ------------------------------
#' 
#' #### Notation for partially observed Markov process models
#' 
#' * Write $X_n=X(t_n)$ and $X_{0:N}=(X_0,\dots,X_N)$. Let $Y_n$ be a random variable modeling the observation at time $t_n$.
#' 
#' * The one-step transition density, $f_{X_n|X_{n-1}}(x_n|x_{n-1};\theta)$, together with the measurement density, $f_{Y_n|X_n}(y_n|x_n;\theta)$ and the initial density, $f_{X_0}(x_0;\theta)$, specify the entire joint density via
#' 
#' $$f_{X_{0:N},Y_{1:N}}(x_{0:N},y_{1:N};\theta) = f_{X_0}(x_0;\theta)\,\prod_{n=1}^N\!f_{X_n | X_{n-1}}(x_n|x_{n-1};\theta)\,f_{Y_n|X_n}(y_n|x_n;\theta).$$
#' 
#' * The marginal density for sequence of measurements, $Y_{1:N}$, evaluated at the data, $y_{1:N}^*$, is
#' 
#' $$ f_{Y_{1:N}}(y^*_{1:N};\theta)=\int f_{X_{0:N},Y_{1:N}}(x_{0:N},y^*_{1:N};\theta)\, dx_{0:N}.$$
#' 
#' ------------------------------
#' 
#' ### Another POMP model schematic
#' 
#' showing dependence among model variables:
#' 
#' 
#' The state process, $X_n$, is Markovian, i.e.,
#' $$\prob{X_n|X_0,\dots,X_{n-1},Y_1,\dots,Y_{n-1}}=\prob{X_n|X_{n-1}}.$$
#' Moreover, the measurements, $Y_n$, depend only on the state at that time:
#' $$\prob{Y_n|X_0,\dots,X_{n},Y_1,\dots,Y_{n-1}}=\prob{Y_n|X_{n}},$$
#' for all $n=1,\dots,N$.
#' 
#' ------------------------------
#' 
#' To think algorithmically, we define some function calls:
#' 
#' * `rprocess( )`: a draw from $f_{X_n|X_{n-1}}(x_n| x_{n-1};\theta)$
#' * `dprocess( )`: evaluation of $f_{X_n|X_{n-1}}(x_n| x_{n-1};\theta)$
#' * `rmeasure( )`: a draw from $f_{Y_n|X_n}(y_n| x_n;\theta)$
#' * `dmeasure( )`: evaluation of $f_{Y_n|X_n}(y_n| x_n;\theta)$
#' * `initializer( )`: a draw from $f_{X_0}(x_0;\theta)$
#' 
#' ### What does it mean for methodology to be __simulation-based__?
#' 
#' * Simulating random processes is often much easier than evaluating their transition probabilities.
#' * In other words, we may be able to write `rprocess()` but not `dprocess()`.
#' *  __Simulation-based__ methods require the user to specify `rprocess()` but not `dprocess()`.
#' * __Plug-and-play__, __likelihood-free__ and __equation-free__ are alternative terms for "simulation-based" methods.
#' * Much development of simulation-based statistical methodology has occurred in the past decade.
#' 
#' 
#' ## The **pomp** package for POMP models
#' 
#' * **pomp** is an  **R**  package for data analysis using partially observed Markov process (POMP) models.
#' 
#' * Note the distinction: lower case **pomp** is a software package; 
#' upper case POMP is a class of models.
#' 
#' * **pomp** builds methodology for POMP models in terms of arbitrary user-specified `rprocess()`, `dprocess()`, `rmeasure()`, and `dmeasure()` functions.
#'  
#' * Following modern practice, most methodology in **pomp** is simulation-based, so does not require specification of `dprocess()`.
#' 
#' * **pomp** has facilities to help construct `rprocess()`, `rmeasure()`, and `dmeasure()` functions for model classes of epidemiological interest.
#' 
#' * **pomp** provides a forum for development, modification and sharing of models, methodology and data analysis workflows.
#' 
#' The following diagrams show the structure of a POMP model schematically.
#' 
#' 
#' ### Example: the Ricker model
#' 
#' #### The deterministic Ricker map
#' 
#' The Ricker map describes the deterministic dynamics of a simple population:
#' $$N_{t+1} = r\,N_{t}\,\exp(-N_{t})$$
#' Here, $N_t$ is the population density at time $t$ and $r$ is a fixed value (a parameter), related to the population's intrinsic capacity to increase.
#' $N$ is a *state variable*, $r$ is a *parameter*.
#' If we know $r$ and the *initial condition* $N_0$, the deterministic Ricker equation predicts the future population density at all times $t=1,2,\dots$.
#' We can view the initial condition, $N_0$ as a special kind of parameter, an *initial-value parameter*.
#' 
#' #### Process noise
#' 
#' We can model process noise in this system by making the growth rate $r$ into a random variable.
#' For example, if we assume that the intrinsic growth rate is log-normally distributed, $N$ becomes a stochastic process governed by
#' $$N_{t+1} = r\,N_{t}\,\exp(-N_{t}+\varepsilon_{t}), \qquad \varepsilon_{t}\;\sim\;\dist{Normal}{0,\sigma},$$
#' where the new parameter $\sigma$ is the standard deviation of the noise process $\varepsilon$.
#' 
#' #### Measurement error
#' 
#' Let's suppose that the Ricker model is our model for the dynamics of a real population.
#' However, we cannot know the exact population density at any time, but only estimate it through sampling.
#' 
#' Let's model measurement error by assuming the measurements, $y_t$, are Poisson with mean $\phi\,N_t$:
#' $$y_{t}\;\sim\;\dist{Poisson}{\phi\,N_{t}}$$
#' 
#' In this equation,
#' 
#' 1. $N_t$ is the true population density at time $t$,
#' 2. $y_t$ is the number of individuals sampled at time $t$,
#' 3. the choice of units for $N$ is peculiar and depends on the parameters (e.g., $N=\log(r)$ is the equilibrium of the deterministic model),
#' 4. the parameter $\phi$ is proportional to our sampling effort, and also has peculiar units.
#' 
#' ### Working with the Ricker model in **pomp**.
#' 
#' The  **R**  package **pomp** provides facilities for modeling POMPs, a toolbox of statistical inference methods for analyzing data using POMPs, and a development platform for implmenting new POMP inference methods.
#' The basic data-structure provided by **pomp** is the *object of class* `pomp`, alternatively known as a `pomp` object.
#' It is a container that holds real or simulated data and a POMP model, possibly together with other information such as model parameters, that may be needed to do things with the model and data.
#' 
#' Let's see what can be done with a `pomp` object.
#' First, we'll load some packages, including **pomp**.
## ----prelims,cache=F-----------------------------------------------------
library(ggplot2)
library(plyr)
library(reshape2)
library(pomp)
stopifnot(packageVersion("pomp")>="1.4.9")

#' 
#' A pre-built `pomp` object encoding the Ricker model comes included with the package.
#' Load it by
## ----load-ricker,cache=FALSE,results="hide"------------------------------
pompExample(ricker)

#' This has the effect of creating a `pomp` object named `ricker` in your workspace.
#' We can plot the data by doing
## ----plot-ricker---------------------------------------------------------
plot(ricker)

#' We can simulate by doing
## ----sim-ricker1---------------------------------------------------------
x <- simulate(ricker)

#' What kind of object have we created?
## ------------------------------------------------------------------------
class(x)
plot(x)

#' Why do we see more time series in the simulated `pomp` object?
#' 
#' We can turn a `pomp` object into a data frame:
## ------------------------------------------------------------------------
y <- as.data.frame(ricker)
head(y)
head(simulate(ricker,as.data.frame=TRUE))

#' 
#' We can also run multiple simulations simultaneously:
## ------------------------------------------------------------------------
x <- simulate(ricker,nsim=10)
class(x)
sapply(x,class)
x <- simulate(ricker,nsim=10,as.data.frame=TRUE)
head(x)
str(x)

#' Also,
## ----fig.height=8--------------------------------------------------------
x <- simulate(ricker,nsim=9,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=x,aes(x=time,y=y,group=sim,color=(sim=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(~sim,ncol=2)

#' 
#' We refer to the deterministic map as the "skeleton" of the stochastic map.
#' We can compute a trajectory of the the deterministic skeleton using `trajectory`:
## ----traj-ricker---------------------------------------------------------
y <- trajectory(ricker)
dim(y)
dimnames(y)
plot(time(ricker),y["N",1,],type="l")

#' 
#' Notice that `ricker` has parameters associated with it:
## ----coef-ricker---------------------------------------------------------
coef(ricker)

#' These are the parameters at which the simulations and deterministic trajectory computations above were done.
#' We can run these at different parameters:
## ------------------------------------------------------------------------
theta <- coef(ricker)
theta[c("r","N.0")] <- c(5,3)
y <- trajectory(ricker,params=theta)
plot(time(ricker),y["N",1,],type="l")
x <- simulate(ricker,params=theta)
plot(x,var="y")

#' 
#' We can also change the parameters stored inside of `ricker`:
## ------------------------------------------------------------------------
coef(ricker,c("r","N.0","sigma")) <- c(39,0.5,1)
coef(ricker)
plot(simulate(ricker),var="y")

#' 
#' In all of the above, it's possible to work with more than one set of parameters at a time.
#' For example:
## ----bifdiag-------------------------------------------------------------
p <- parmat(coef(ricker),500)
dim(p); dimnames(p)
p["r",] <- seq(from=2,to=40,length=500)
y <- trajectory(ricker,params=p,times=200:1000)
matplot(p["r",],y["N",,],pch=".",col='black',xlab='r',ylab='N',log='x')

#' 
#' How do you interpret the above plot?
#' This is called a *one-parameter bifurcation diagram".
#' 
#' More information on manipulating and extracting information from `pomp` objects can be viewed in the help pages (`methods?pomp`).
#' 
#' There are a number of other examples included with the package.
#' Do `pompExample()` to see a list of these.
#' More examples can be found in the **pompExamples** package:
#' 
#' ### Inference algorithms in **pomp**
#' 
#' **pomp** provides a wide range of inference algorithms.
#' We'll learn about these in detail soon, but for now, let's just look at some of their general features.
#' 
#' The `pfilter` function runs a simple particle filter.
#' It can be used to evaluate the likelihood at a particular set of parameters.
#' One uses the `Np` argument to specify the number of particles to use:
## ----pfilter1------------------------------------------------------------
pf <- pfilter(ricker,Np=1000)
class(pf)
plot(pf)
logLik(pf)

#' 
#' Note that `pfilter` returns an object of class `pfilterd.pomp`.
#' This is the general rule: inference algorithms return objects that are `pomp` objects with additional information.
#' The package provides tools to extract this information.
#' We can run the particle filter again by doing
## ----pfilter2------------------------------------------------------------
pf <- pfilter(pf)
logLik(pf)

#' which has the result of running the same computation again.
#' Note that, because the particle filter is a Monte Carlo algorithm, we get a slightly different estimate of the log likelihood.
#' 
#' Note that, by default, running `pfilter` on a `pfilterd.pomp` object causes the computation to be re-run with the same parameters as before.
#' Any additional arguments we add override these defaults.
#' This is the general rule in **pomp**.
#' For example,
## ----pfilter3------------------------------------------------------------
pf <- pfilter(pf,Np=100)
logLik(pf)

#' Here, the particle filtering has been performed with only `r unique(pf@Np)` particles.
#' 
#' ### Building a custom `pomp` object
#' 
#' A real **pomp** data analysis begins with constructing one or more `pomp` objects to hold the data and the model or models under consideration.
#' We'll illustrate this process a dataset of *Parus major* abundance in Wytham Wood, near Oxford [@McCleery1991].
#' 
#' Download and plot the data:
## ----parus-data----------------------------------------------------------
loc <- url("http://kingaa.github.io/sbied/intro/parus.csv")
dat <- read.csv(loc)
head(dat)
plot(pop~year,data=dat,type='o')

#' 
#' Let's suppose that we want to fit the stochastic Ricker model discussed above to these data.
#' 
#' The call to construct a `pomp` object is, naturally enough, `pomp`.
#' Documentation on this function can be had by doing `?pomp`. 
#' Learn about the various things you can do once you have a `pomp` object by doing `methods?pomp` and following the links therein.
#' Read an overview of the package as a whole with links to its main features by doing `package?pomp`.
#' A complete index of the functions in **pomp** is returned by the command `library(help=pomp)`.
#' Finally, the home page for the `pomp` project is http://kingaa.github.io/pomp;
#' there you have access to the complete source code, tutorials, manuals, issues page, news blog, etc.
#' 
#' Now, to construct our `pomp` object:
## ----parus-pomp1---------------------------------------------------------
library(pomp)
parus <- pomp(dat,times="year",t0=1959)

#' The `times` argument specifies that the column of `dat` labelled "year" gives the measurement times;
#' `t0` is the "zero-time", the time at which the state process will be initialized.
#' We've set it to one year prior to the beginning of the data.
#' Plot the new `pomp` object:
## ----parus-plot1---------------------------------------------------------
plot(parus)

#' 
#' #### Adding in the process model simulator
#' 
#' We can add the stochastic Ricker model to `parus` by writing a Csnippet that simulates one realization of the stochastic process, from an arbitary time $t$ to $t+1$, given arbitrary states and parameters.
#' We provide this to `pomp` in the form of a `Csnippet`, a little snippet of C code that performs the computation.
#' The following does this.
## ----parus-sim-defn------------------------------------------------------
stochStep <- Csnippet("
  e = rnorm(0,sigma);
  N = r*N*exp(-N+e);
")
pomp(parus,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
     paramnames=c("r","sigma"),statenames=c("N","e")) -> parus

#' Note that in the above, we use the `exp` and `rnorm` functions from the [**R** API](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#The-R-API).
#' In general any C function provided by **R** is available to you.
#' **pomp** also provides a number of C functions that are documented in the header file, `pomp.h`, that is installed with the package.
#' See the `Csnippet` documentation (`?Csnippet`) to read more about how to write them.
#' Note too that we use `discrete.time.sim` here because the model is a stochastic map.
#' We specify that the time step of the discrete-time process is `delta.t`, here, 1&nbsp;yr.
#' 
#' At this point, we have what we need to simulate the state process:
## ----ricker-first-sim----------------------------------------------------
sim <- simulate(parus,params=c(N.0=1,e.0=0,r=12,sigma=0.5),
                as.data.frame=TRUE,states=TRUE)
plot(N~time,data=sim,type='o')

#' 
#' #### Adding in the measurement model and parameters
#' 
#' We complete the specification of the POMP by specifying the measurement model.
#' To obtain the Poisson measurement model described above, we write two Csnippets.
#' The first simulates:
## ----parus-rmeas-defn----------------------------------------------------
rmeas <- Csnippet("pop = rpois(phi*N);")

#' The second computes the likelihood of observing `pop` birds given a true density of `N`:
## ----parus-dmeas-defn----------------------------------------------------
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")

#' [Note the `give_log` argument.
#' When this code is evaluated, `give_log` will be set to 1 if the log likelihood is desired, and 0 else.]
#' We add these into the `pomp` object:
## ----parus-add-meas------------------------------------------------------
pomp(parus,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> parus

#' Now we can simulate the whole POMP.
#' First, let's add some parameters:
## ----ricker-add-params---------------------------------------------------
coef(parus) <- c(N.0=1,e.0=0,r=20,sigma=0.1,phi=200)

## ----ricker-second-sim,results='markup'----------------------------------
library(ggplot2)
sims <- simulate(parus,nsim=3,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=sims,mapping=aes(x=time,y=pop))+geom_line()+
  facet_wrap(~sim,ncol=1,scales="free_y")

#' 
#' #### Adding in the deterministic skeleton
#' 
#' We can add the Ricker model deterministic skeleton to the `parus` `pomp` object.
#' Since the Ricker model is a discrete-time model, its skeleton is a map that takes $N_t$ to $N_{t+1}$ according to the Ricker model equation
#' $$N_{t+1} = r\,N_{t}\,\exp(-N_{t}).$$
#' The following implements this.
## ----parus-skel-defn-----------------------------------------------------
skel <- Csnippet("DN = r*N*exp(-N);")

#' We then add this to the `pomp` object:
## ----parus-add-skel------------------------------------------------------
parus <- pomp(parus,skeleton=map(skel),paramnames=c("r"),statenames=c("N"))

#' Note that we have to inform **pomp** as to which of the variables we've referred to in `skel` is a state variable (`statenames`) and which is a parameter (`paramnames`).
#' In writing a `Csnippet` for the deterministic skeleton, we use `D` to designate the map's value.
#' The `map` call tells **pomp** that the skeleton is a discrete-time dynamical system (a map) rather than a continuous-time system (a vectorfield).
#' 
#' With just the skeleton defined, we are in a position to compute the trajectories of the deterministic skeleton at any point in parameter space.
#' For example, here we compute the trajectory and superimpose it on a plot of one simulation:
## ----parus-first-traj,results='markup'-----------------------------------
traj <- trajectory(parus,params=c(N.0=1,r=12),as.data.frame=TRUE)
plot(N~time,data=sim,type='o')
lines(N~time,data=traj,type='l',col='red')

#' 
#' #### A note on terminology
#' 
#' If we know the state, $x(t_0)$, of the system at time $t_0$, it makes sense to speak about the entire trajectory of the system for all $t>t_0$.
#' This is true whether we are thinking of the system as deterministic or stochastic.
#' Of course, in the former case, the trajectory is uniquely determined by $x(t_0)$, while in the stochastic case, only the probability distribution of $x(t)$, $t>t_0$ is determined.
#' To avoid confusion, we use the term "trajectory" exclusively to refer to *trajectories of a deterministic process*.
#' Thus, the `trajectory` command iterates or integrates the deterministic skeleton forward in time, returning the unique trajectory determined by the specified parameters.
#' When we want to speak about sample paths of a stochastic process, we use the term *simulation*.
#' Accordingly, the `simulate` command always returns individual sample paths from the POMP.
#' In particular, we avoid "simulating a set of differential equations", preferring instead to speak of "integrating" the equations, or "computing trajectories".
#' 
#' ## Exercises
#' 
#' --------------------------
#' 
#' #### Exercise: Ricker model parameters
#' Fiddle with the parameters to try and make the simulations look more like the data.
#' This will help you build some intuition for what the various parameters do.
#' 
#' --------------------------
#' 
#' #### Exercise: Reformulating the Ricker model
#' Reparameterize the Ricker model so that the scaling of $N$ is explicit:
#' $$N_{t+1} = r\,N_{t}\,\exp\left(-\frac{N_{t}}{K}+\varepsilon_t\right).$$
#' 
#' Modify the `pomp` object we created above to reflect this reparameterization.
#' 
#' Modify the measurement model so that
#' $$\mathrm{pop}_t \sim \dist{Negbin}{\phi\,N_t,k},$$
#' i.e., $\mathrm{pop}_t$ is negative-binomially distributed with mean $\phi\,N_t$ and clumping parameter $k$.
#' See `?NegBinomial` for documentation on the negative binomial distribution and [the **R** Extensions Manual section on distribution functions](http://cran.r-project.org/doc/manuals/r-release/R-exts.html#Distribution-functions) for information on how to access these in C.
#' 
#' --------------------------
#' 
#' #### Exercise: Beverton-Holt
#' Construct a `pomp` object for the *Parus major* data and the stochastic Beverton-Holt model
#' $$N_{t+1} = \frac{a\,N_t}{1+b\,N_t}\,\varepsilon_t,$$
#' where $a$ and $b$ are parameters and
#' $$\varepsilon_t \sim \dist{Lognormal}{-\tfrac{1}{2}\sigma^2,\sigma}.$$
#' Assume the same measurement model as before.
#' 
#' ------------------------------
#' 
#' ## [Back to course homepage](http://kingaa.github.io/sbied)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/gh-pages/intro/intro.R)
#' 
#' ----------------------
#' 
#' ## References
