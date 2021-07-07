switch(
  system("hostname",intern=TRUE),
  theorygroup01={
    library(doMPI)
    cl <- startMPIcluster(250,verbose=TRUE,logdir="/tmp")
    registerDoMPI(cl)
  },
  {
    library(doParallel)
    registerDoParallel()
  }
  )
