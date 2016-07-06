#' ---
#' title: "Data munging with **plyr**, **reshape2**, and **magrittr**"
#' author: "Aaron A. King"
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 4
#' bibliography: ../course.bib
#' csl: ../ecology.csl
#' ---
#' 
#' 
#' ## How to use this document.
#' 
#' This is an extremely condensed introduction to the powerful data-munging tools developed by Hadley Wickham and contained in the packages **plyr**, **reshape2**, and **magrittr**.
#' Run the codes shown and study the outputs to learn about these tools.
#' For your convenience, the [**R** codes for this document are provided in a script](http://raw.githubusercontent.com/kingaa/short-course/gh-pages/hadley/viz.R) which you can download, edit, and run.
#' 
#' ## Reshaping data with **reshape2**
#' 
#' The **reshape2** package works with a metaphor of *melting* and *casting*.
#' 
#' ### Melting
#' 
#' Melting takes a wide data frame and makes it long.
#' Multiple columns are combined into one *value* column with a *variable* column keeping track of which column the different values came from.
#' Only the columns containing *measure* variables are reshaped;
#' those containing *identifier* variables are left alone.
#' 
## ------------------------------------------------------------------------
library(reshape2)

x <- data.frame(a=letters[1:10],b=1:10,
                c=sample(LETTERS[1:3],10,replace=TRUE),d=sample(1:10,10,replace=T))
x
melt(x,id.vars=c("a","b"))
melt(x,measure.vars=c("c","d")) -> y; y

#' 
#' ### Casting
#' 
#' Casting turns a long data frame into a wide one.
#' A single column (called the *value* column) is separated into multiple columns according to the specification given.
#' Use `dcast` or `acast` according to whether you want the result as a data frame or an array.
## ------------------------------------------------------------------------
dcast(y,a+b~variable) -> d1; d1
class(d1)
acast(y,b~variable) -> a1; a1
class(a1); dim(a1)
acast(y,a~b~variable) -> a2; a2
class(a2); dim(a2)

#' 
#' ## Split-apply-combine with **plyr**
#' 
#' **plyr** implements a very flexible and intuitive syntax for split-apply-combine computations.
#' That is, it allows you to split data according to a wide range of criteria, apply some operation to each piece, them recombine the pieces back together.
#' 
#' In the following, we first detail the "basic" functions that make up the "apply" piece of split-apply-combine.
#' Then, we discuss the "split" and "combine" pieces.
#' 
#' ### Basic **plyr** functions
#' 
#' The following are the basic functions for manipulating data using **plyr**.
#' 
#' #### `arrange`
#' 
#' `arrange` sorts a data frame according to specifications.
#' 
## ------------------------------------------------------------------------
library(plyr)

x <- data.frame(a=letters[1:10],b=runif(10),c=sample(LETTERS[1:3],10,replace=TRUE))
arrange(x,a,b,c)
arrange(x,b,c,a)
arrange(x,c,b,a)

#' 
## ------------------------------------------------------------------------
read.csv("http://kingaa.github.io/short-course/hadley/energy_production.csv",comment="#") -> energy
arrange(energy,region,source,year)
arrange(energy,-TJ,year)

#' 
#' #### `count`
#' 
#' `count(x)` counts the combinations that occur and returns a data frame.
## ------------------------------------------------------------------------
count(x,~c)
count(x,~a+c)
count(x,vars=c('a','c'))

#' 
## ------------------------------------------------------------------------
count(energy,~source+region)
count(energy,~source+TJ)

#' 
#' #### `summarise` and `summarize`
#' 
#' Given a data frame, `summarise` (synonym `summarize`), produces a new data frame.
## ------------------------------------------------------------------------
summarize(x,mean=mean(b),sd=sd(b),top=c[1])

#' 
## ------------------------------------------------------------------------
summarize(energy,tot=sum(TJ),n=length(TJ))
summarize(energy,range(year))
summarize(energy,min(year),max(year),interval=diff(range(year)))

#' 
#' #### `mutate`
#' Given a data frame, `mutate` modifies, adds, or removes variables.
## ------------------------------------------------------------------------
x <- mutate(x,d=2*b,c=tolower(c),e=b+d,a=NULL); x

#' 
#' #### `subset`
#' 
#' `subset` doesn't belong to **plyr**, but would if it didn't already exist in the **base** package.
#' This function allows you to choose a subset of rows and/or columns.
#' The `subset` argument specifies a logical condition: those rows that satisfy it are chosen.
#' The `select` argument picks out which columns to keep or throw away.
## ------------------------------------------------------------------------
subset(x,d>1.2)
subset(x,select=c(b,c))
subset(x,select=-c(d))
subset(x,d>1.2,select=-e)

## ------------------------------------------------------------------------
subset(energy,year>2010,select=c(source,TJ))
subset(energy,year>2010&source%in%c("Nuclear","Oil"),select=-source)

#' 
#' #### `merge` and `join`
#' 
#' `merge` belongs to the **base** package; 
#' `join` belongs to **plyr**.
#' They both do versions of the database *join* operation.
#' 
## ------------------------------------------------------------------------
x <- expand.grid(a=1:3,b=1:5)
y <- expand.grid(a=1:2,b=1:5,c=factor(c("F","G")))
m1 <- merge(x,y); m1
m2 <- merge(x,y,by='a'); m2
m3 <- merge(x,y,all=TRUE); m3
m4 <- merge(x,y,by='a',all=TRUE); m4

#' 
#' `join` is more general implementing the *database join operations*.
#' It can perform a *left join*, a *right join*, an *inner join*, or a *full join*.
#' Read the documentation (`?join`) for explanations.
#' 
## ------------------------------------------------------------------------
join(x,y,by=c('a','b'),type='left')
join(x,y,by=c('a','b'),type='right')
join(x,y,by=c('a','b'),type='inner')
join(x,y,by=c('a','b'),type='full')
join(x,y,by='a',type='full')
join(x,y,by='a',type='inner')

#' 
#' 
#' ### The `-ply` functions
#' 
#' **plyr** provides a systematic, intuitive, and regular expansion of base **R**'s `apply` family (`apply`, `lapply`, `sapply`, `tapply`, `mapply`) and `replicate`.
#' Collectively, these functions implement the split-apply-combine pattern of computation.
#' They first split the data up according to some criterion, then apply some function, then combine the results.
#' The functions are all named according to the scheme `XYply`, where `X` tells about the class of the source object and `Y` the class of the desired target object.
#' In particular `X` and `Y` can be in `d` (data-frames), `a` (arrays), `l` (lists), `_` (null), and `r` (replicate).
#' 
#' #### `ddply`
#' 
#' This is probably the most useful of the lot.
#' It splits a data frame according to some criterion, conveniently expressed as a formula involving the variables of the data frame, applies a specified function, and combines the results back into a data frame.
#' It is best to use a function that returns a data frame, but if the function returns something else, `ddply` will attempt to coerce the value into a data frame.
#' Here are some examples:
## ------------------------------------------------------------------------
x <- ddply(energy,~region+source,subset,TJ==max(TJ)); x
x <- ddply(energy,~region+source,summarize,TJ=mean(TJ)); x

#' Notice that only combinations of the variables that exist are included in the result by default.
#' 
#' #### `daply`
#' 
#' This one is very similar, except that (as the name implies), the result is returned as an array:
## ------------------------------------------------------------------------
daply(energy,~region,function(df) sum(df$TJ))
daply(energy,~region+source,function(df) sum(df$TJ))

#' 
#' #### `dlply`
#' 
#' This splits the data according to the given specifications, applies the function, and returns each result (as its name implies) as a distinct element of a list.
## ------------------------------------------------------------------------
dlply(energy,~region,summarize,TJ=sum(TJ))

#' 
#' #### `adply`, `aaply`, `alply`
#' 
#' These take arrays and, like the **base** function `apply`, divide the array up into slices along specified directions.
#' They then apply a function to each slice and return the results in the desired form (if possible).
#' As an example, we first create an array from `dat`, then act on it with each of these.
## ------------------------------------------------------------------------
mutate(energy,time=year-min(year)) -> dat
daply(dat,~source+region,function(df) min(df$time)) -> A; A
aaply(A,1,max)

#' 
#' #### Exercise
#' Create some simple arrays and practice using these functions.
#' 
#' 
#' #### `llply`, `laply`, `ldply`
#' 
#' These functions are generalizations of `lapply` and `sapply`.
#' 
#' #### Exercise
#' Create a few simple lists and practice using these functions.
#' 
#' 
#' #### `mlply`, `maply`, `mdply`
#' 
#' These work with multi-argument functions.
#' 
#' #### Exercise
#' Create a simple data frame and practice using these functions.
#' 
#' 
#' ### Other functions
#' 
#' #### `rename`, `revalue`, `mapvalues`
#' 
#' `rename` helps one to change the (column) names of a data frame.
## ------------------------------------------------------------------------
x <- rename(energy,c(TJ='energy',year="time")); head(x)

#' 
#' `revalue` allows you to change one or more of the levels of a factor without worrying about how the factors are coded.
#' 
#' `mapvalues` does the same, but works on vectors of any type.
## ------------------------------------------------------------------------
mutate(energy,region=revalue(region,c(`Asia and Oceania`="Asia",
                                      `Central and South America`="Latin.America"))); 

mutate(energy,source=mapvalues(source,from=c("Coal","Gas","Oil"),
                               to=c("Carbon","Carbon","Carbon")))

#' 
#' ## The **magrittr** syntax
#' 
#' ![ceci n'est pas une pipe](MagrittePipe.png)  
#' Ren&eacute; Magritte, *La Trahison des Images*
#' 
#' 
#' **magrittr** gives a set of "pipe" operators.
#' These allow one to chain operations together.
#' When calculations get complex, it is easier and more natural to view them as a chain of operations instead of using nested function calls or defining intermediate variables.
#' 
#' ### The `%>%` operator
#' 
#' ```
#' f(g(data, a, b, c, ...), d, e, ...)
#' ```
#' 
#' is equivalent to
#' 
#' ```
#' data %>% g(a, b, c, ...) %>% f(d, e, ...)
#' ```
#' 
#' ### The `%<>%` operator
#' 
#' ```
#' x %>% f(a, b, c, ...) -> x
#' ```
#' is equivalent to
#' ```
#' x %<>% f(a, b, c, ...)
#' ```
#' 
## ------------------------------------------------------------------------
library(magrittr)

energy %>% 
  subset(year>=1990) %>%
  ddply(~source+year,summarize,TJ=sum(TJ)) %>%
  ddply(~source,summarize,TJ=mean(TJ))

#' 
#' --------------------------
#' 
#' ## [Back to course homepage](http://kingaa.github.io/short-course)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/short-course/gh-pages/hadley/viz.R)
#' 
#' --------------------------
#' 
#' ## References
#' 
