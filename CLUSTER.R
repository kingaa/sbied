switch(
  system("hostname",intern=TRUE),
  theorygroup01={
    library(iterators)
    library(doFuture)
    cl <- makeClusterMPI(250,verbose=TRUE,logdir="/tmp")
    registerDoFuture()
    plan(cluster,workers=cl)
  },
  {
    library(iterators)
    library(doFuture)
    registerDoFuture()
    plan(multisession)
  }
  )
