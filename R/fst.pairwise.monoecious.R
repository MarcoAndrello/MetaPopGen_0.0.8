# Pairwise FST 
fst.pairwise.monoecious = function(N,i,j,xi,xj,ti,tj) {
  
  # Define 1-D array of genotype counts for the first and second group of individuals
  # Calculate number of individuals in each group 
  
  if (is.character(N)) {
    dir.res.name <- paste(getwd(),N,sep="/")
    
    # Load file for first group of individuals
    name.file <- paste(dir.res.name,"/N",ti,".RData",sep="")
    load(name.file)
    Ni <- N[,i,xi]
    n_i <- sum(N[,i,xi])
    
    # Load file for second group of individuals
    name.file <- paste(dir.res.name,"/N",tj,".RData",sep="")
    load(name.file)
    Nj <- N[,j,xj]
    n_j  <- sum(N[,j,xj])
    
  } else {
    
    Ni <- N[,i,xi,ti]
    n_i <- sum(N[,i,xi,ti])
    
    Nj <- N[,j,xj,tj]
    n_j  <- sum(N[,j,xj,tj])
    
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
