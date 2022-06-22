switch(
  system("hostname",intern=TRUE),
  theorygroup01={
    library(doMPI)
    cl <- startMPIcluster(250,verbose=TRUE,logdir="/tmp")
    registerDoMPI(cl)
    cores <- 250
  },
  {
    library(doParallel)
    cores <- as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE"))
    if (is.na(cores)) cores <- detectCores()  
    registerDoParallel(cores)
  }
)
