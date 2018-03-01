# Pairwise FST 
fst.pairwise.dioecious = function(N,i,j,xi,xj,ti,tj,si,sj) {
  
  
  # Define 1-D array of genotype counts for the first and second group of individuals
  # Calculate number of individuals in each group 
  
  if (is.character(N)) {
    dir.res.name <- paste(getwd(),N,sep="/")
    
    # Load file for first group of individuals
    name.file <- paste(dir.res.name,"/N",ti,".RData",sep="")
    load(name.file)
    if (si == "M") {
      Ni <- N_M[,i,xi]
      n_i <- sum(N_M[,i,xi])
    } else {
      Ni <- N_F[,i,xi]
      n_i <- sum(N_F[,i,xi])
    }
    
    # Load file for second group of individuals
    name.file <- paste(dir.res.name,"/N",tj,".RData",sep="")
    load(name.file)
    if (sj == "M") {
      Nj <- N_M[,j,xj]
      n_j  <- sum(N_M[,j,xj])
    } else {
      Nj <- N_F[,j,xj]
      n_j  <- sum(N_F[,j,xj])
    }

    
  } else {
    if (si == "M") {
      Ni <- N[[1]][,i,xi,ti]
      n_i <- sum(N[[1]][,i,xi,ti])
    } else {
      Ni <- N[[2]][,i,xi,ti]
      n_i <- sum(N[[2]][,i,xi,ti])
    }
    
    if (sj == "M") {
      Nj <- N[[1]][,j,xj,tj]
      n_j  <- sum(N[[1]][,j,xj,tj])
    } else {
      Nj <- N[[2]][,j,xj,tj]
      n_j  <- sum(N[[2]][,j,xj,tj])
    }
  }
  
  # Calculate number of individuals in the total "population"
  # Population means the two groups together
  n_T	<- n_i + n_j
  
  # Calculate allele frequencies in each group and in the total "population"
  p_i	<- freq.all(Ni)
  p_j	<- freq.all(Nj)
  p_T	<- ( n_i*p_i + n_j*p_j ) / n_T
  
  # Calculate expected heterozygosities in each group and in the total "population"
  H_S_i	<- 1 - sum(p_i^2)
  H_S_j	<- 1 - sum(p_j^2)
  H_t	<- 1 - sum(p_T^2)
  
  # Calculate fst as ratio of heterozygosities
  fst <- ( H_t - (n_i * H_S_i + n_j * H_S_j)/(n_i+ n_j) ) / H_t
  return(fst)
  
}
