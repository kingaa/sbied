#' ---
#' title: Likelihood for POMP models
#' author: "Aaron A. King and Edward L. Ionides"
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 4
#' bibliography: ../sbied.bib
#' csl: ../ecology.csl
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
#' 
#' --------------------------
#' 
#' [Licensed under the Creative Commons Attribution-NonCommercial license](http://creativecommons.org/licenses/by-nc/4.0/).
#' Please share and remix noncommercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
#' Produced with **R** version `r getRversion()` and **pomp** version `r packageVersion("pomp")`.
#' 
#' --------------------------
#' 
#' 
## ----prelims,cache=FALSE,include=FALSE-----------------------------------
library(pomp)
stopifnot(packageVersion("pomp")>="1.6")
library(ggplot2)
theme_set(theme_bw())
library(plyr)
library(reshape2)
options(stringsAsFactors=FALSE)
set.seed(1221234211)

#' 
#' ## Objectives
#' 
#' Students completing this lesson will:
#' 
#' 1. Gain an understanding of the nature of the problem of likelihood computation for POMP models.
#' 1. Be able to explain the simplest particle filter algorithm.
#' 1. Gain experience in the visualization and exploration of likelihood surfaces.
#' 1. Learn to apply standard optimization methods for the maximization of the likelihood.
#' 
#' ## Theory of the particle filter
#' 
#' ### The likelihood function
#' 
#' - The basis for modern frequentist, Bayesian, and information-theoretic inference.
#' - Method of maximum likelihood introduced by @Fisher1922.
#' - The function itself is a representation of the what the data have to say about the parameters.
#' - A good general reference on likelihood is @Pawitan2001.
#' 
#' #### Definition of the likelihood function
#' 
#' Data are a sequence of $N$ observations, denoted $y_{1:N}^*$.
#' A statistical model is a density function $f(y_{1:N};\theta)$ which defines a probability distribution for each value of a parameter vector $\theta$.
#' To perform statistical inference, we must decide, among other things, for which (if any) values of $\theta$ it is reasonable to model $y^*_{1:N}$ as a random draw from $f(y_{1:N};\theta)$.
#' 
#' The likelihood function is the density function evaluated at the data.
#' It is usually convenient to work with the log likelihood function,
#' $$\loglik(\theta)=\log f(y^*_{1:N};\theta)$$
#' 
#' #### Modeling using discrete and continuous distributions
#' 
#' Recall that the probability distribution $f(y_{1:N};\theta)$ defines a random variable $Y_{1:N}$ for which probabilities can be computed as integrals of $f(y_{1:N};\theta)$.
#' Specifically, for any event $E$ describing a set of possible outcomes of $Y_{1:N}$, 
#' $$P[Y_{1:N} \in E] = \int_E f(y_{1:N};\theta)\, dy_{1:N}.$$ 
#' If the model corresponds to a discrete distribution, then the integral is replaced by a sum and the probability density function is called a *probability mass function*.
#' The definition of the likelihood function remains unchanged.
#' We will use the notation of continuous random variables, but all the methods apply also to discrete models. 
#' 
#' #### Indirect specification of the statistical model via a simulation procedure
#' 
#' - For simple statistical models, we may describe the model by explicitly writing the density function $f(y_{1:N};\theta)$. 
#' One may then ask how to simulate a random variable $Y_{1:N}\sim f(y_{1:N};\theta)$.
#' - For many dynamic models it is much more convenient to define the model via a procedure to simulate the random variable $Y_{1:N}$. 
#' This implicitly defines the corresponding density $f(y_{1:N};\theta)$. 
#' For a complicated simulation procedure, it may be difficult or impossible to write down or even compute $f(y_{1:N};\theta)$ exactly. 
#' - It is important to bear in mind that the likelihood function exists even when we don't know what it is!
#' We can still talk about the likelihood function, and develop numerical methods that take advantage of its statistical properties.
#' 
#' 
#' ### Likelihood for POMP models
#' 
#' *********************
#' The following schematic shows dependence among variables in a POMP model.
#' Measurements, $Y_n$, at time $t_n$ depend on the state, $X_n$, at that time.
#' State variables depend on state variables at the previous timestep.
#' To be more precise, the distribution of the state $X_{n+1}$, conditional on $X_{n}$, is independent of the values of $X_{k}$, $k<n$ and $Y_{k}$, $k\le n$.
#' Moreover, the distribution of the measurement $Y_{n}$, conditional on $X_{n}$, is independent of all other variables.
#' 
#' #### POMP model notation
#' 
#' - Write $X_n=X(t_n)$ and $X_{0:N}=(X_0,\dots,X_N)$.
#' - Let $Y_n$ be a random variable modeling the observation at time $t_n$.
#' - The one-step transition density, $f_{X_n|X_{n-1}}(x_n|x_{n-1};\theta)$, together with the measurement density, $f_{Y_n|X_n}(y_n|x_n;\theta)$ and the initial density, $f_{X_0}(x_0;\theta)$, specify the entire joint density via
#' $$f_{X_{0:N},Y_{1:N}}(x_{0:N},y_{1:N};\theta) = f_{X_0}(x_0;\theta)\,\prod_{n=1}^N\!f_{X_n | X_{n-1}}(x_n|x_{n-1};\theta)\,f_{Y_n|X_n}(y_n|x_n;\theta).$$
#' - The marginal density for sequence of measurements, $Y_{1:N}$, evaluated at the data, $y_{1:N}^*$, is
#' $$\lik(\theta) = f_{Y_{1:N}}(y^*_{1:N};\theta)=\int\!f_{X_{0:N},Y_{1:N}}(x_{0:N},y^*_{1:N};\theta)\, dx_{0:N}.$$
#' 
#' #### Special case: deterministic unobserved state process
#' 
#' Lets' begin with a special case.
#' Suppose that the unobserved state process is deterministic.
#' That is, $X_{n}=x_n(\theta)$ is a known function of $\theta$ for each $n$.
#' What is the likelihood?
#' 
#' Since the probability of the observation, $Y_n$, depends only on $X_n$ and $\theta$, and since, in particular $Y_{m}$ and $Y_{n}$ are independent given $X_{m}$ and $X_{n}$, we have $$\lik(\theta) = \prod_{n} f_{Y_n|X_n}(y_n^*;x_n(\theta),\theta)$$ or $$\ell(\theta) = \log\lik(\theta) = \sum_{n} \log f_{Y_n|X_n}(y_n^*;x_n(\theta),\theta).$$
#' The following diagram illustrates this.
#' 
#' 
#' In this diagram, $\hat y_n$ refers to the model prediction ($\hat y_n = \expect{Y_n \vert X_n=x_n(\theta)}$) and $y_n^*$ to the data.
#' 
#' #### General case: stochastic unobserved state process
#' 
#' For a POMP model, the likelihood takes the form of an integral:
#' $$\lik(\theta)=\prob{y^*_{1:N}|\theta}=\int\!\prod_{n=1}^{N}\!\prob{y^*_n|x_n,\theta}\,\prob{x_n|x_{n-1},\theta}\,dx_{1}\,\cdots\,dx_{N}.$$
#' This integral is high dimensional and, except for the simplest cases, can not be reduced analytically.
#' A general approach to computing integrals, which will be useful here is [Monte Carlo integration](./monteCarlo.html#monte-carlo-integration).
#' 
#' [[See here for a brief introduction to Monte Carlo methods](./monteCarlo.html).]
#' 
#' Let us investigate a Monte Carlo integration scheme for the POMP likelihood.
#' In particular, if we propose trajectories of the unobserved state process according to some probabilistic rule, we can generate a large number of these and approximate $\lik(\theta)$ by its Monte Carlo estimate.
#' Specifically, let us randomly generate $J$ trajectories of length $N$, $x_{n,j}$, $j=1\,\dots,J$, $n=1,\dots,N$.
#' Let $w_j$ denote the probability that we propose trajectory $j$.
#' We compute the likelihood of each trajectory
#' $$\lik{_j}(\theta)=\prod_{n=1}^{N} \prob{y^*_n|x_{n,j},\theta}\,\prob{x_{n,j}|x_{n-1,j},\theta}$$
#' Then by the Monte Carlo theorem, we have 
#' $$\lik(\theta) \approx \frac{1}{J}\,\sum_{j=1}^{J}\!\frac{\lik_j(\theta)}{w_j}.$$
#' 
#' How shall we choose our trajectories?
#' One idea would be to choose them so as to simplify the computation.
#' If we choose them such that
#' $$w_j=\prod_{n=1}^{N} \prob{x_{n,j}|x_{n-1,j},\theta},$$
#' then we have
#' $$\lik(\theta) \approx \frac{1}{J}\,\sum_{j=1}^{J} \frac{\lik_j(\theta)}{w_j} = \frac{1}{J}\,\sum_{j=1}^{J}\!\frac{\prod_{n=1}^{N} \prob{y^*_n|x_{n,j},\theta}\,\prob{x_{n,j}|x_{n-1,j},\theta}}{\prod_{n=1}^{N} \prob{x_{n,j}|x_{n-1,j},\theta}} = \frac{1}{J}\,\sum_{j=1}^{J}\!\prod_{n=1}^{N} \prob{y^*_n|x_{n,j},\theta}$$
#' 
#' This implies that if we generate trajectories by simulation, all we have to do is compute the likelihood of the data with given each trajectory and then average.
#' 
#' Let's go back to the boarding school influenza outbreak to see what this looks like in practice.
#' Let's reconstruct the toy SIR model we were working with.
#' 
## ----flu-construct-------------------------------------------------------
read.table("http://kingaa.github.io/sbied/stochsim/bsflu_data.txt") -> bsflu

rproc <- Csnippet("
  double N = 763;
  double t1 = rbinom(S,1-exp(-Beta*I/N*dt));
  double t2 = rbinom(I,1-exp(-mu_I*dt));
  double t3 = rbinom(R1,1-exp(-mu_R1*dt));
  double t4 = rbinom(R2,1-exp(-mu_R2*dt));
  S  -= t1;
  I  += t1 - t2;
  R1 += t2 - t3;
  R2 += t3 - t4;
")

init <- Csnippet("
  S = 762;
  I = 1;
  R1 = 0;
  R2 = 0;
")

dmeas <- Csnippet("
  lik = dpois(B,rho*R1+1e-6,give_log);
")

rmeas <- Csnippet("
  B = rpois(rho*R1+1e-6);
")

pomp(subset(bsflu,select=-C),
     times="day",t0=0,
     rprocess=euler.sim(rproc,delta.t=1/5),
     initializer=init,rmeasure=rmeas,dmeasure=dmeas,
     statenames=c("S","I","R1","R2"),
     paramnames=c("Beta","mu_I","mu_R1","mu_R2","rho")) -> flu

#' 
#' Let's generate a large number of simulated trajectories at some particular point in parameter space.
## ----bbs-mc-like-2-------------------------------------------------------
simulate(flu,params=c(Beta=3,mu_I=1/2,mu_R1=1/4,mu_R2=1/1.8,rho=0.9),
         nsim=5000,states=TRUE) -> x
matplot(time(flu),t(x["R1",1:50,]),type='l',lty=1,
        xlab="time",ylab=expression(R[1]),bty='l',col='blue')
lines(time(flu),obs(flu,"B"),lwd=2,col='black')

#' 
#' We can use the function `dmeasure` to evaluate the log likelihood of the data given the states, the model, and the parameters:
## ----bbs-mc-like-3,cache=T-----------------------------------------------
ell <- dmeasure(flu,y=obs(flu),x=x,times=time(flu),log=TRUE,
                params=c(Beta=3,mu_I=1/2,mu_R1=1/4,mu_R2=1/1.8,rho=0.9))
dim(ell)

#' According to the equation above, we should sum up the log likelihoods across time:
## ----bbs-mc-like-4-------------------------------------------------------
ell <- apply(ell,1,sum)
summary(exp(ell))

#' The variability in the individual likelihoods is high and therefore the likelihood esitmate is imprecise.
#' We will need many simulations to get an estimate of the likelihood sufficiently precise to be of any use in parameter estimation or model selection.
#' 
#' What is the problem?
#' Essentially, very few of the trajectories pass anywhere near the data and therefore almost all have extremely bad likelihoods.
#' Moreover, once a trajectory diverges from the data, it almost never comes back.
#' While the calculation is "correct" in that it will converge to the true likelihood as the number of simulations tends to $\infty$, we waste a lot of effort investigating trajectories of very low likelihood.
#' *This is a consequence of the fact that we are proposing trajectories in a way that is completely unconditional on the data.*
#' The problem will get much worse with longer data sets.
#' 
#' ### The particle filter
#' 
#' We arrive at a more efficient algorithm by factorizing the likelihood in a different way:
#' $$\lik(\theta)=\prob{y^*_{1:N}|\theta}
#' =\prod_{n}\,\prob{y^*_n|y^*_{1:n-1},\theta} 
#' =\prod_{n}\,\int\!\prob{y^*_n|x_n,\theta}\,\prob{x_n|y^*_{1:n-1},\theta}\,dx_{n}.\tag{1}$$
#' Now, the Markov property gives us the Chapman-Kolmogorov equation,
#' $$\prob{x_n|y^*_{1:n-1},\theta} 
#' = \int\!\prob{x_n|x_{n-1},\theta}\,\prob{x_{n-1}|y^*_{1:n-1},\theta}\,dx_{n-1},\tag{2}$$
#' and Bayes' theorem tells us that
#' $$\prob{x_{n}|y^*_{1:n},\theta} = \prob{x_{n}|y^*_{n},y^*_{1:n-1},\theta} =\frac{\prob{y^*_{n}|x_{n},\theta}\,\prob{x_{n}|y^*_{1:n-1},\theta}}{\displaystyle\int\!\prob{y^*_{n}|x_{n},\theta}\,\prob{x_{n}|y^*_{1:n-1},\theta}\,dx_{n}}.\tag{3}$$
#' 
#' This suggests that we keep track of two key distributions.
#' We'll refer to the distribution of $X_n | y^*_{1:n-1}$ as the *prediction distribution* at time $n$ and
#' the distribution of $X_{n} | y^*_{1:n}$ as the *filtering distribution* at time $n$.
#' 
#' Let's use Monte Carlo techniques to estimate the integrals.
#' Suppose $\left\{x_{n-1,j}^{F}\right\}_{j=1}^J$ is a set of points drawn from the filtering distribution at time $n-1$.
#' Eqn.&nbsp;2 tells us that we obtain a sample $\left\{x_{n,j}^{P}\right\}$ of points drawn from the prediction distribution at time $n$ by simply simulating the process model:
#' $$X_{n,j}^{P} \sim \mathrm{process}(x_{n-1,j}^{F},\theta), \qquad j=1,\dots,J.$$
#' Having obtained $\left\{x_{n,j}^{P}\right\}$, we obtain a sample of points from the filtering distribution at time $n$ by *resampling* from $\left\{x_{n,j}^{P}\right\}$ with weights proportional to $\prob{y^*_{n}|x_{n},\theta}$.
#' The Monte Carlo theorem tells us, too, that the conditional likelihood 
#' $$\lik_n(\theta) = \prob{y^*_n|y^*_{1:n-1},\theta} = \sum_{x_{n}}\,\prob{y^*_{n}|x_{n},\theta}\,\prob{x_{n}|y^*_{1:n-1},\theta} \approx \frac{1}{J}\,\sum_j\,\prob{y^*_{n}|x_{n,j}^{P},\theta}.$$
#' We can iterate this procedure through the data, one step at a time, alternately simulating and resampling, until we reach $n=N$.
#' The full log likelihood is then approximately
#' $$\loglik(\theta) = \log{\lik(\theta)} \approx \sum_n \log{\lik_n(\theta)}.$$
#' It can be shown that this estimate of the likelihood is unbiased.
#' 
#' This is known as the *sequential Monte Carlo* algorithm or the *particle filter*.
#' Key references include @Kitagawa1987, @Arulampalam2002, and the book by @Doucet2001.
#' Pseudocode for the above is provided by @King2016.
#' 
#' 
#' ## Sequential Monte Carlo in **pomp**
#' 
#' Here, we'll get some practical experience with the particle filter, and the likelihood function, in the context of our influenza-outbreak case study.
#' 
#' In **pomp**, the basic particle filter is implemented in the command `pfilter`.
#' We must choose the number of particles to use by setting the `Np` argument.
#' 
## ----flu-pfilter-1,cache=T-----------------------------------------------
pf <- pfilter(flu,Np=5000,params=c(Beta=3,mu_I=1/2,mu_R1=1/4,mu_R2=1/1.8,rho=0.9))
logLik(pf)

#' 
#' We can run a few particle filters to get an estimate of the Monte Carlo variability:
## ----flu-pfilter-2,cache=T-----------------------------------------------
pf <- replicate(10,
                pfilter(flu,Np=5000,
                        params=c(Beta=3,mu_I=1/2,mu_R1=1/4,mu_R2=1/1.8,rho=0.9)))
ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)

#' 
#' ## The graph of the likelihood function: The likelihood surface
#' 
#' It is extremely useful to visualize the geometric surface defined by the likelihood function.
#' 
#' - If $\Theta$ is two-dimensional, then the surface $\ell(\theta)$ has features like a landscape.
#' - Local maxima of $\ell(\theta)$ are peaks.
#' - Local minima are valleys.
#' - Peaks may be separated by a valley or may be joined by a ridge. 
#' If you go along the ridge, you may be able to go from one peak to the other without losing much elevation. 
#' Narrow ridges can be easy to fall off, and hard to get back on to.
#' - In higher dimensions, one can still think of peaks and valleys and ridges. 
#' However, as the dimension increases it quickly becomes hard to imagine the surface.
#' 
#' To get an idea of what the likelihood surface looks like in the neighborhood of a point in parameter space, we can construct some likelihood *slices*.
#' We'll make slices in the $\beta$ and $\mu_I$ directions.
#' Both slices will pass through the central point.
#' 
## ----flu-like-slice,cache=TRUE,results='hide'----------------------------
sliceDesign(
  center=c(Beta=2,mu_I=1,mu_R1=1/4,mu_R2=1/1.8,rho=0.9),
  Beta=rep(seq(from=0.5,to=4,length=40),each=3),
  mu_I=rep(seq(from=0.5,to=2,length=40),each=3)) -> p

library(foreach)
library(doParallel)
registerDoParallel()

set.seed(108028909,kind="L'Ecuyer")

foreach (theta=iter(p,"row"),.combine=rbind,
         .inorder=FALSE,
         .options.multicore=list(set.seed=TRUE)
) %dopar% {
  library(pomp)
  pfilter(flu,params=unlist(theta),Np=5000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> p

#' 
#' Note that we've used the **foreach** package with the parallel backend (**doParallel**) to parallelize these computations.
#' To ensure that we have high-quality random numbers in each parallel *R* session, we use a parallel random number generator (`kind="L'Ecuyer"`, `.options.multicore=list(set.seed=TRUE)`).
#' 
## ----flu-like-slice-plot,cache=FALSE,echo=FALSE--------------------------
library(magrittr)
library(reshape2)
library(ggplot2)
p %>% 
  melt(measure=c("Beta","mu_I")) %>%
  subset(variable==slice) %>%
  ggplot(aes(x=value,y=loglik,color=variable))+
  geom_point()+
  facet_grid(~variable,scales="free_x")+
  guides(color=FALSE)+
  labs(x="parameter value",color="")+
  theme_bw()

#' 
#' ------------------------
#' 
#' #### Exercise: Likelihood slice
#' 
#' Add likelihood slices along the $\rho$ direction.
#' 
#' ------------------------
#' 
#' Slices offer a very limited perspective on the geometry of the likelihood surface.
#' With just two parameters, we can evaluate the likelihood at a grid of points and visualize the surface directly.
## ----flu-grid1-----------------------------------------------------------
bake(file="flu-grid1.rds",seed=421776444,kind="L'Ecuyer",{
  
  expand.grid(Beta=seq(from=1.5,to=5,length=50),
              mu_I=seq(from=0.7,to=4,length=50),
              mu_R1=1/4,mu_R2=1/1.8,
              rho=0.9) -> p
  
  library(foreach)
  library(doParallel)
  registerDoParallel()
  
  ## Now we do the computation
  foreach (theta=iter(p,"row"),.combine=rbind,.inorder=FALSE,
           .options.multicore=list(set.seed=TRUE)
  ) %dopar% 
  {
    library(pomp)
    pfilter(flu,params=unlist(theta),Np=5000) -> pf
    theta$loglik <- logLik(pf)
    theta
  }
  
})-> p

## ----flu-grid1-plot,echo=F,purl=T----------------------------------------
library(magrittr)
library(reshape2)
library(plyr)
p %<>% arrange(Beta,mu_I,mu_R1,mu_R2,rho)
saveRDS(p,file="flu-grid1.rds")
p %>% 
  mutate(loglik=ifelse(loglik>max(loglik)-50,loglik,NA)) %>%
  ggplot(aes(x=Beta,y=mu_I,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient()+
  labs(x=expression(beta),y=expression(mu[I]))

#' 
#' In the above, all points with log likelihoods less than 50 units below the maximum are shown in grey.
#' 
#' ------------------------
#' 
#' #### Exercise: 2D likelihood slice
#' 
#' Compute a slice of the likelihood in the $\beta$-$\rho$ plane.
#' 
#' ------------------------
#' 
#' 
#' ## Maximizing the likelihood
#' 
#' Call the whole parameter space $\Theta$. 
#' Let $\Theta^*$ be a subset of $\Theta$, constraining parameters according to a scientific hypothesis of interest. 
#' For example, in a disease transmission model, $\Theta^*$ could assert that the probability of a case being reported is $\rho=0.8$.
#' 
#' - We define the maximized log likelihoods for $\Theta$ and $\Theta^*$ to be
#' $$\ell_\mathrm{max}=\max\{\ell(\theta):\theta\in\Theta\},\quad
#' \ell^*_\mathrm{max}=\max\{\ell(\theta): \theta\in\Theta^*\}.$$
#' - Intuitively, a model with a higher maximized likelihood should be preferable to a model with a substantially lower maximized likelihood. 
#' - However, since $\Theta^*$ is a subset of $\Theta$, it is mathematically necessary that
#' $$ \ell_\mathrm{max} \ge \ell^*_\mathrm{max}.$$
#' This raises the question of how close $\ell^*_\mathrm{max}$ should be to $\ell_\mathrm{max}$ to make it reasonable to prefer the simpler model $\Theta^*$ over the more complex model $\Theta$.
#' - The principle of parsimony (Occam's razor) advises that we be satisfied with the simpler model unless there is good evidence to do otherwise. 
#' For a formal hypothesis test, we accordingly set our null hypothesis to be $\Theta^*$ and our alternative hypothesis to be $\Theta$.
#' - A likelihood ratio test rejects $\Theta^*$ in favor of $\Theta$ when 
#' $$ \ell_\mathrm{max} -\ell_\mathrm{max}^* > c.$$
#' - An elegant mathematical property (Wilks' theorem) says that, for regular parametric models where $N$ is large and $\Theta$ has $d$ more free parameters than $\Theta^*$, then $2(\ell_\mathrm{max}-\ell^*_\mathrm{max})$ has a chi-square distribution with $d$ degrees of freedom.
#' - For the concrete situation where $\Theta^*$ fixes a single parameter, $d=1$, and we look for a test of size $0.05$, this suggests we reject $\Theta^*$ if
#' $$\ell_\mathrm{max} - \ell^*_\mathrm{max} > 1.92$$
#' since $\prob{\chi^2_1>3.84}=0.05$.
#' - One can carry out a simulation study to assess the actual size of this test, if one is concerned whether the asymptotic property of Wilks is sufficiently accurate. 
#' Fortunately, Wilks' theorem is often a good approximation for many finite-sample problems.
#' - Wilks' theorem gives a convenient, quick scientific interpretation of maximized log likelihood values.
#' One can choose later whether to refine the interpretation via further simulation studies.
#' - Akaike's information criterion (AIC) is defined by
#' $$\mathrm{AIC} = -2(\mbox{maximized log likelihood}) +2(\mbox{# parameters}).$$
#' This criterion makes a slightly different decision, recommending $\Theta$ over $\Theta^*$ if $$\ell_\mathrm{max} -\ell_\mathrm{max}^* > d.$$
#' The justification of AIC is based on minimizing prediction error. AIC tends to prefer larger models than Occam's razor: 
#' heuristically, it values simplicity not for its own sake, but only because unnecessary parameters lead to over-fitting and hence greater out-of-fit forecasting error. 
#' - Wilks' theorem applies only to nested hypotheses (when $\Theta^*$ is a subset of $\Theta$) whereas AIC is applicable to compare non-nested models, which may have entirely different structure. 
#' - Although AIC is not designed to be a formal statistical test, it is a commonly used objective rule for model selection.
#' This rule could be intrepreted as a hypothesis test, with the size and power investigated by simulation, if desired.
#' 
#' --------------------------
#' 
#' #### Exercise: AIC as a formal statistical test
#' 
#' Determine the size of AIC as a hypothesis test for nested hypotheses with $d=1$ in a regular parametric situation. 
#' 
#' --------------------------
#' 
#' ### Point estimates for parameters: The maximum likelihood estimate (MLE)
#' 
#' We define maximum likelihood estimates (MLEs) $\hat\theta$ and $\hat\theta^*$ such that
#' $$\ell(\hat\theta)=\ell_\mathrm{max},\quad \ell(\hat\theta^*)=\ell_\mathrm{max}^*.$$
#' 
#' - If the likelihood function has a flat region, or ridge, at its maximum then the MLE is not unique. 
#' Alternatively, one can talk about a maximum likelihood surface describing the set of parameter values for which $\ell(\hat\theta)=\ell_\mathrm{max}$.
#' - Flat, or nearly flat, ridges in the likelihood surface are not an idle concern. 
#' Many dynamic models have combinations of parameters that are weakly identified: they cannot be well estimated on the basis of the data.
#' 
#' ### Confidence intervals for parameters: Profile likelihood
#' 
#' The likelihood ratio test with $d=1$ gives a good way to construct confidence intervals. Suppose we are interested in a specific parameter, $\theta_k$, and we want to consider whether the data support the possibility that $\theta_k=\theta_k^*$ in the absence of assumptions on the other parameters.
#' We can then take $\Theta^*$ to be the subset of $\Theta$ satisfying $\theta_k=\theta_k^*$.
#' Using the chi-square approximation to the likelihood ratio statistic, a 95% confidence interval for $\theta_k$ consists of all the values $\theta_k^*$ for which
#' $$\ell_\mathrm{max}-\ell_\mathrm{max}^* < 1.92.$$
#' 
#' A way to visualize the information about a specific parameter $\theta_k$ is via the profile likelihood function, defined as 
#' $$\ell_\mathrm{profile}(\theta_k^*) = \max\{\ell(\theta): \theta_k=\theta_k^*\}.$$
#' We then plot $\ell_\mathrm{profile}(\theta_k)$ against $\theta_k$. 
#' 
#' - The set of values of $\theta_k$ for which $\ell_\mathrm{profile}(\theta_k)$ lies above a horizontal line with $y$-axis value $\ell_\mathrm{max}-c$ gives an approximate confidence interval (according to Wilks' theorem) with confidence level given by $\prob{\chi^2_1<2c}$.
#' - The maximum of $\ell_\mathrm{profile}(\theta_k)$ over all values of $\theta_k$ is $\ell_\mathrm{max}$.
#' - Thus, a profile plot allows us to visualize an entire spectrum of confidence intervals.
#' - If the profile plot has two peaks (i.e., $\ell_\mathrm{profile}(\theta_k)$ is bimodal) then a likelihood ratio test helps us to assess whether or not both peaks provide adequate explanations of the data.
#' 
#' ## Maximizing the likelihood using the particle filter
#' 
#' In the toy example we've been working with, the default parameter set is not particularly close to the MLE.
#' One way to find the MLE is to try optimizing the estimated likelihood directly.
#' There are of course many standard optimization algorithms we might use for this.
#' However, three issues arise immediately:
#' 
#' 1. The particle filter gives us a stochastic estimate of the likelihood.
#' We can reduce this variability by making $J$ larger, but we cannot make it go away.
#' If we use a deterministic optimizer (i.e., one that assumes the objective function is evaluated deterministically), then we must control this variability somehow.
#' For example, we can fix the seed of the pseudo-random number generator (RNG).
#' A side effect will be that the objective function becomes jagged, marked by many small local knolls and pits.
#' Alternatively, we can use a stochastic optimization algorithm, with which we will be only be able to obtain estimates of our MLE.
#' This is the trade-off between a rough and a noisy objective function.
#' 1. Because the particle filter gives us just an estimate of the likelihood and no information about the derivative, we must choose an algorithm that is "derivative-free".
#' There are many such, but we can expect less efficiency than would be possible with derivative information.
#' Note that finite differencing is not an especially promising way of constructing derivatives. 
#' The price would be a $n$-fold increase in cpu time, where $n$ is the dimension of the parameter space.
#' Also, since the likelihood is noisily estimated, we would expect the derivative estimates to be even noisier.
#' 1. Finally, the parameters set we must optimize over is not unbounded.
#' We must have $\beta,\mu_I>0$ and $0<\rho<1$.
#' We must therefore select an optimizer that can solve this *constrained maximization problem*, or find some of way of turning it into an unconstrained maximization problem.
#' For example, we can transform the parameters onto a scale on which there are no constraints.
#' 
#' Here, let's opt for deterministic optimization of a rough function.
#' We'll try using `optim`'s default method: Nelder-Mead, fixing the random-number generator seed to make the likelihood calculation deterministic.
#' Since Nelder-Mead is an unconstrained optimizer, we must transform the parameters.
#' The following `Csnippet`s encode an appropriate transformation and its inverse, and introduce them into the `pomp` object.
## ----flu-partrans--------------------------------------------------------
toEst <- Csnippet("
 TBeta = log(Beta);
 Tmu_R1 = log(mu_R1);
 Tmu_I = log(mu_I);
 Trho = logit(rho);
")

fromEst <- Csnippet("
 TBeta = exp(Beta);
 Tmu_I = exp(mu_I);
 Tmu_R1 = exp(mu_R1);
 Trho = expit(rho);
")

pomp(flu,toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     paramnames=c("Beta","mu_I","mu_R1","rho")) -> flu

#' 
#' Let's fix a reference point in parameter space and insert these parameters into the `pomp` object:
## ----flu-ref-params------------------------------------------------------
coef(flu) <- c(Beta=2,mu_I=1,mu_R1=512/sum(bsflu$B),mu_R2=512/sum(bsflu$C),rho=0.9)

#' 
#' The following constructs a function returning the negative log likelihood of the data at a given point in parameter space.
#' The parameters to be estimated are named in the `est` argument.
#' Note how the `freeze` function is used to fix the seed of the RNG.
#' Note too, how this function returns a large (and therefore bad) value when the particle filter encounters and error.
#' This behavior makes the objective function more robust.
#' 
## ----flu-like-optim-1----------------------------------------------------
neg.ll <- function (par, est) {
  allpars <- coef(flu,transform=TRUE)
  allpars[est] <- par
  try(
    freeze(
      pfilter(flu,params=partrans(flu,allpars,dir="fromEst"),
              Np=2000),
      seed=915909831
    )
  ) -> pf
  if (inherits(pf,"try-error")) 1e10 else -logLik(pf)
}

#' 
#' Now we call `optim` to minimize this function:
## ----flu-like-optim-2----------------------------------------------------
## use Nelder-Mead with fixed RNG seed
fit <- optim(
  par=c(log(2), log(1), log(0.9/(1-0.9))),
  est=c("Beta","mu_I","rho"),
  fn=neg.ll,
  method="Nelder-Mead",
  control=list(maxit=400,trace=0)
)

mle <- flu
coef(mle,c("Beta","mu_I","rho"),transform=TRUE) <- fit$par
coef(mle)

fit$val

lls <- replicate(n=5,logLik(pfilter(mle,Np=20000)))
ll <- logmeanexp(lls,se=TRUE); ll

#' 
#' We plot some simulations at these parameters.
## ----flu-sims------------------------------------------------------------
simulate(mle,nsim=10,as.data.frame=TRUE,include.data=TRUE) -> sims

#' The data are shown in blue.
#' The `r max(sims$sim)` simulations are shown in red.
## ----flu-sims-plot,echo=F------------------------------------------------
ggplot(data=sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  guides(color=FALSE)+
  geom_line()

#' 
#' 
#' --------------------------
#' 
#' #### Exercise: More slices
#' 
#' Construct likelihood slices through the MLE we just found.
#' 
#' --------------------------
#' 
#' #### Exercise: Visualizing the likelihood surface
#' 
#' Evaluate the likelihood at points on a grid lying in a 2D slice through the MLE we found above.
#' Each group should choose a different slice.
#' Afterward, we'll compare results across groups.
#' 
#' --------------------------
#' 
#' #### Exercise: Global maximization
#' 
#' The search of parameter space we conducted above was local.
#' It is possible that we found a local maximum, but that other maxima exist with higher likelihoods.
#' Conduct a more thorough search by initializing the Nelder-Mead starting points across a wider region of parameter space.
#' Do you find any other local maxima?
#' 
#' --------------------------
#' 
#' #### Exercise: Modify the measurement model
#' 
#' The Poisson measurement model used here may not seem easy to interpret.
#' Formulate an alternative measurement model and maximize the likelihood to compare the alternative model.
#' 
#' --------------------------
#' 
#' #### Exercise: Fit more parameters.
#' 
#' Try to estimate $\beta$, $\mu_I$, $\rho$, and $\mu_{R1}$ simultaneously.
#' Does your estimate of $\mu_{R1}$ differ from the value we computed from the raw data?
#' How do you interpret the agreement or lack thereof?
#' 
#' --------------------------
#' 
#' ## Biological interpretation of parameter estimates
#' 
#' When we write down a mechanistic model for an epidemiological system, we have some idea of what we intend parameters to mean; a reporting rate, a contact rate between individuals, an immigration rate, a duration of immunity, etc. 
#' 
#' - The data and the parameter estimation procedure do not know about our intended interpretation of the model. 
#' It can and does happen that some parameter estimates, statistically consistent with the data, may be scientifically absurd according to the biological reasoning that went into building the model. 
#' - This can arise as a consequence of weak identifiability. 
#' - It can also be a warning that the data do not agree that our model represents reality in the way we had hoped.
#' This is a signal that more work is needed on model development.
#' - Biologically unreasonable parameter estimates can sometimes be avoided by fixing some parameters at known, reasonable values. 
#' However, this risks suppressing the warning that the data were trying to give about weaknesses in the model, or in the biological interpretation of it.
#' - This issue will be discussed further in connection with the case studies.
#' 
#' --------------------------
#' 
#' ## [Back to course homepage](http://kingaa.github.io/sbied)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/gh-pages/pfilter/pfilter.R)
#' 
#' --------------------------
#' 
#' ## References
