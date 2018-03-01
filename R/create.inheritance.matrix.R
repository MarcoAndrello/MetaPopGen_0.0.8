create.inheritance.matrix <- function(index, List_gene, r=0.5, mu) {
  
  g <- length(List_gene) # Number of loci
  m <- max(index)        # Number of multi-locus genotypes
  l <- prod(List_gene)   # Number of multi-locus alleles
  
  # Create GENO matrix to know which allele comes from which genotype - could be made more efficient...
  GENO <- array(NA,dim=c(m,2*g))    # Matrix: for each genotype (lines), gives which allele at each locus
  #(locus 1: 1st and 2nd column; locus 2: rd and 4th column; etc)
  names <- c()
  for(i in 1 : l){
    for(ii in 1 : l){
      if(is.na(GENO[index[i,ii],1])){
        w <- 0
        for (h in 1 : g){
          a <- c(substr(rownames(index)[i],2*h,2*h),substr(colnames(index)[ii],2*h,2*h))
          GENO[index[i,ii],c(h+w,2*h)] <- as.numeric(a)
          w <- w+1
        }
        names <- c(names,paste0(rownames(index)[i],"/",colnames(index)[ii]))
      }
    }
  }
  dimnames(GENO) <- list(names,rep(1:2,g))

  for (j in 1 : m) {
    genotype <- strsplit(rownames(GENO)[j],"/")
    alleles.1 <- strsplit(genotype[[1]][1],"")[[1]][seq(2,(2*g),2)]
    alleles.2 <- strsplit(genotype[[1]][2],"")[[1]][seq(2,(2*g),2)]
    tabl.all <- rbind(as.numeric(alleles.1),as.numeric(alleles.2))
    colnames(tabl.all) <- LETTERS[1:g]
    freq.all.in.genotype <- rep(0,l)
    #### mi fermo qua
    
    for (k in 1 : l) {
      gamete <- rownames(index)[k]
      
      }
    }
  
  
  
  }