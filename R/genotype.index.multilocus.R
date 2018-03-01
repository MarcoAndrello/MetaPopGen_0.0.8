genotype.index.multilocus <- function(List_gene){
  
  #Given list of number of alleles per locus, the number of gametes is prod(List_gene)
  l <- prod(List_gene)
  
  y <- array(NA,dim=c(l,l))
  diagy <- array(NA,dim=l)
  diagy[1] <- 1
  for (i in 2 : l) {
    diagy[i] <- diagy[i-1] + (l-(i-2))
  }
  diag(y)<-diagy
  for (i in 1 : (l-1)) {
    y[i,i:l] <- seq(y[i,i],(y[i,i]+l-i))
    y[i:l,i] <- seq(y[i,i],(y[i,i]+l-i))
  }
  
  

  # Set dimension names, for <=26 loci
  nloc <- length(List_gene)
  alleles_list <- list()
  for (loc in 1 : nloc) {
      alleles <- vector("character",List_gene[loc])
      for (i in 1 : List_gene[loc]){
          alleles[i] <- paste0(LETTERS[loc],i)
      }
      alleles_list[[loc]] <- alleles
  }
  alleles_list <- expand.grid(alleles_list)
  alleles_list <- apply(alleles_list,2,as.character)
  
  dn <- apply(alleles_list,1,paste0,collapse="")
  dimnames(y) <- list(dn,dn)
  
  return(y)
}

