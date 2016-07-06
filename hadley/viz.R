#' ---
#' title: "Data visualization"
#' author: Aaron A. King
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 4
#' bibliography: ../course.bib
#' csl: ../ecology.csl
#' ---
#' 
## ----prelims,include=FALSE,cache=FALSE-----------------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
  )

set.seed(594709947L)
library(ggplot2)
theme_set(theme_bw())

#' 
#' ## How to use this document.
#' 
#' This is an extremely condensed introduction to **R**'s base graphics and---more importantly---the powerful data-visualization package **ggplot2**, developed by Hadley Wickham.
#' Run the codes shown and study the outputs to learn about these tools.
#' When questions are posed, do your best to answer them.
#' 
#' For your convenience, [the **R** codes for this document are provided in an **R** script](http://raw.githubusercontent.com/kingaa/short-course/gh-pages/hadley/viz.R) which you can download, edit, and run.
#' 
#' ## Getting started: **R**'s base graphics
#' 
#' ### Transgenic mosquito experiment
#' 
#' Let's load the data on transgenic mosquito survival time.
#' 
## ------------------------------------------------------------------------
dat <- read.csv("http://kingaa.github.io/short-course/hadley/mosquitoes.csv")

#' 
#' Let's compare the average lifespan of transgenic vs wildtype mosquitoes from this experiment.
#' The following split the data into two subsets, one for each genetic type.
## ------------------------------------------------------------------------
wt <- subset(dat,type=="wildtype",select=lifespan)
tg <- subset(dat,type=="transgenic",select=-type)

#' 
#' Let's try and visualize the data.
## ------------------------------------------------------------------------
plot(dat)
op <- par(mfrow=c(1,2))
hist(tg$lifespan,breaks=seq(0,55,by=5),ylim=c(0,40))
hist(wt$lifespan,breaks=seq(0,55,by=5),ylim=c(0,40))
par(op)

#' 
#' **Question:** What does the second `par` command accomplish?
#' 
#' Another way to visualize a distribution is via the *empirical cumulative distribution plot*.
#' 
## ------------------------------------------------------------------------
plot(sort(dat$lifespan),seq(1,nrow(dat))/nrow(dat),type='n')
lines(sort(wt$lifespan),seq(1,nrow(wt))/nrow(wt),type='s',col='blue')
lines(sort(tg$lifespan),seq(1,nrow(tg))/nrow(tg),type='s',col='red')

#' **Question:** What does `type="n"` do in the first line above?
#' 
#' ### Mammal body and brain sizes
#' 
#' The data on mammal body and brain sizes is included in the **MASS** package:
## ------------------------------------------------------------------------
library(MASS)

plot(mammals)
plot(mammals,log='x')
plot(mammals,log='xy')
plot(mammals$body,mammals$brain,log='xy')
plot(brain~body,data=mammals,log='xy')

#' 
#' ### Oil production
#' 
## ------------------------------------------------------------------------
read.csv("http://kingaa.github.io/short-course/hadley/oil_production.csv",
         comment.char="#") -> oil
head(oil)
summary(oil)
plot(oil)
plot(Gbbl~year,data=oil,subset=region=="North.America",type='l')
lines(Gbbl~year,data=oil,subset=region=="Eurasia",type="l",col='red')

library(reshape2)

dcast(oil,year~region) -> wideOil
names(wideOil)
wideOil$total <- wideOil$Africa+wideOil$Asia+wideOil$Central+wideOil$Eurasia+wideOil$Europe+wideOil$Middle+wideOil$North.America
wideOil$total <- apply(wideOil[,-1],1,sum)
plot(wideOil$year,wideOil$total,type='l')

#' 
#' ## A systematic approach to visualization: the Grammar of Graphics
#' 
#' Parts of a graphic:
#' 
#' 1. ***Data***
#' 1. ***Geometrical object***: point, line, box, bar, density plot, contours, ribbons
#' 1. ***Statistical transformations***: bins, mean, median, quantile, ECDF, identity 
#' 1. ***Aesthetic attributes***: x and y position, color, fill, size, shape, line type, transparency
#' 1. ***Scales***: map the data onto the aesthetic attributes
#' 1. A ***coordinate system***: maps x and y position onto the page
#' 1. A ***faceting system***: multiple plots
#' 
#' You construct a graphical visualization by choosing the constituent parts.
#' This is implemented in the **ggplot2** package.
#' 
#' ### References
#' 
#' - [ggplot2.org](http://ggplot2.org)
#' - [ggplot2 documentation](http://docs.ggplot2.org/)
#' 
#' 
#' ## Examples
#' 
#' ### Energy production
#' 
## ------------------------------------------------------------------------
read.csv("http://kingaa.github.io/short-course/hadley/energy_production.csv",
         comment.char="#") -> energy

library(ggplot2)

ggplot(data=energy,mapping=aes(x=year,y=TJ,color=region,linetype=source))+geom_line()
ggplot(data=energy,mapping=aes(x=year,y=TJ,color=region))+geom_line()+facet_wrap(~source)
ggplot(data=energy,mapping=aes(x=year,y=TJ,color=source))+geom_line()+facet_wrap(~region,ncol=2)

#' 
#' What can you conclude from the above?
#' Try plotting these data on the log scale (`scale_y_log10()`).
#' How does your interpretation change?
#' 
## ------------------------------------------------------------------------
ggplot(data=energy,mapping=aes(x=year,y=TJ))+geom_line()
ggplot(data=energy,mapping=aes(x=year,y=TJ,group=source))+geom_line()

#' 
#' **Question:** How do you account for the appearance of the two plots immediately above?
#' 
## ------------------------------------------------------------------------
ggplot(data=energy,mapping=aes(x=year,y=TJ,group=source:region))+geom_line()

#' 
#' **Question:** What does the `group` aesthetic do?
#' 
#' Let's aggregate across regions by year and source of energy.
## ------------------------------------------------------------------------
library(reshape2)

tot <- dcast(energy,year+source~'TJ',value.var="TJ",fun.aggregate=sum)
ggplot(data=tot,mapping=aes(x=year,y=TJ,color=source))+geom_line()
ggplot(data=tot,mapping=aes(x=year,y=TJ,fill=source))+geom_area()


#' 
#' Now let's aggregate across years by region and source.
#' 
## ------------------------------------------------------------------------
reg <- dcast(energy,region+source~'TJ',value.var="TJ",fun.aggregate=mean)
ggplot(data=reg,mapping=aes(x=region,y=TJ,fill=source))+
   geom_bar(stat="identity")+coord_flip()

#' 
#' An even better way to manipulate the data is to use the **plyr** package.
#' [See the data munging tutorial.](./data_munging.html)
#' 
## ------------------------------------------------------------------------
library(plyr)

ddply(energy,~region+source,summarize,TJ=mean(TJ)) -> x

ggplot(data=x,mapping=aes(x=region,y=TJ,fill=source))+
   geom_bar(stat="identity")+coord_flip()

ddply(x,~region,mutate,frac=TJ/sum(TJ)) -> y

ggplot(data=y,mapping=aes(x=region,y=frac,fill=source))+
   geom_bar(stat="identity")+coord_flip()+labs(x="fraction of production")


#' 
#' In the above, we first average across years for every region and source.
#' Then, for each region, we compute the fraction of the total production due to each source.
#' Finally, we plot the fractions using a barplot.
#' The `coord_flip` coordinate specification gives us horizontal bars instead of the default vertical bars.
#' Fancy!
#' 
#' Let's compare fossil fuel production to renewable.
## ------------------------------------------------------------------------
library(plyr)

mutate(energy,
       source=as.character(source),
       source1=mapvalues(source,
                         from=c("Hydro","Other Renewables","Coal","Oil","Gas"),
                         to=c("Renewable","Renewable","Carbon","Carbon","Carbon"))
       ) -> energy

ddply(energy,~source1+region+year,summarize,TJ=sum(TJ)) -> x

ggplot(data=x,mapping=aes(x=year,y=TJ,fill=source1))+
    geom_area()+
    facet_wrap(~region,scales="free_y",ncol=2)

ddply(energy,~source1+year,summarize,TJ=sum(TJ)) -> x

ggplot(data=x,mapping=aes(x=year,y=TJ,fill=source1))+
    geom_area()

#' --------------------------
#' 
#' ### Exercise
#' 
#' Ask a question regarding one of the datasets shown here and devise a visualization to answer it.
#' 
#' --------------------------
#' 
#' ## [Back to course homepage](http://kingaa.github.io/short-course)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/short-course/gh-pages/hadley/viz.R)
#' 
#' --------------------------
