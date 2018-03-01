# Global fst, dioecious species
fst.global.dioecious <- function(N,s,t) {
  
    if (is.character(N)) {
		
		dir.res.name <- paste(getwd(),N,sep="/")
		name.file <- paste(N,"/N",t,".RData",sep="")
		load(name.file)
		
		if (s == "M") {
			N <- N_M
		} else {
			N <- N_F
		}		
		
    
	} else {
		
		if (s == "M") {
			N <- N[[1]][,,,t]
		} else {
			N <- N[[2]][,,,t]
		}
		
    
  }
  

  
  # Define basic variables
  m <- dim(N)[1]            # Number of genotypes
  l <- (sqrt(1+8*m)-1)/2    # Number of alleles
  n <- dim(N)[2]            # Number of demes
  
  # Genotype arrays and number of individuals for each deme
  N_i <- apply(N,c(1,2),sum)
  n_i <- apply(N,2,sum)
  
  # Genotype arrays and number of individuals for total population
  N_T <- apply(N,1,sum)
  n_T  <- sum(N)
    
  # Calculate allele frequencies in each group and in the total population
  p_i <- array(NA,dim=c(n,l))
  for (i in 1 : n) {
    p_i[i,] <- freq.all(N_i[,i])    
  }
  p_T	<- freq.all(N_T)
  
  # Calculate expected heterozygosities in each group and in the total population
  H_S_i <- array(NA,dim=n)
  for (i in 1 : n) {
    H_S_i[i] <- 1 - sum(p_i[i,]^2)
  }
  H_t	<- 1 - sum(p_T^2)
  
  # Calculate fst as ratio of heterozygosities
  fst <- ( H_t - sum(n_i * H_S_i,na.rm=T)/n_T)  / H_t
  return(fst)
  
}