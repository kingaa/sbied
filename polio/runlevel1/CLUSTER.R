switch(
  system("hostname",intern=TRUE),
  theorygroup01={
    library(doFuture)
    library(iterators)
    cl <- makeClusterMPI(250,verbose=TRUE,logdir="/tmp")
    plan(cluster,workers=cl)
  },
  {
    library(doFuture)
    library(iterators)
    plan(multicore)
  }
  )
