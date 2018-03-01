# Creation of the matrix of probability used in reproduction.
# This matrix gives the probabilities for each gamete type as a funciton of the parental genotypes, mutation and recombination rates

create.probability.matrix <- function(index, List_gene, r=0.5, mu){
  
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
  
  # Create RECOMB matrix
  RECOMB<-array(0, dim=c(l,m)) # Matrix giving the probability of production of each gameto-type (lines) for each genotype (columns)
  
  dimnames(RECOMB)<-list(colnames(index),rownames(GENO))
  
  if (length(List_gene)==2){
    
    TYPE<-expand.grid(1:List_gene[1],1:List_gene[2]) # seems to give the gameto-type: which allele for each locus for each line of RECOMB
    colnames(TYPE) <- LETTERS[1:g] ## I added this
    for (k in 1:m){
      genotype<-GENO[k,]
      for(j in 1:l){
        if ((genotype[1]== TYPE[j,1]) & (genotype[3]==TYPE[j,2])){
          RECOMB[j,k]<-RECOMB[j,k]+((1/2)*(1-r))
        }
        if ((genotype[1]== TYPE[j,1]) & (genotype[4]==TYPE[j,2])){
          RECOMB[j,k]<-RECOMB[j,k]+(r/2)
        }
        if ((genotype[2]== TYPE[j,1]) & (genotype[3]==TYPE[j,2])){
          RECOMB[j,k]<-RECOMB[j,k]+(r/2)
        }
        if ((genotype[2]== TYPE[j,1]) & (genotype[4]==TYPE[j,2])){
          RECOMB[j,k]<-RECOMB[j,k]+((1/2)*(1-r))
        }
      }
    }
  } else {
    p<-2^g
    
    for(k in 1:m){
      if(substr(rownames(GENO)[k],1,2*g)==substr(rownames(GENO)[k],2*g+2,4*g+1)){ # If it is a homozygote
        P<-which(colnames(index)==substr(rownames(GENO)[k],1,2*g)) # .. it produces only one gametotype, equal to either gamete forming the genotype
        RECOMB[P,k]<-1
      }
      else {
        alleles_list <- list()
        w<-0
        for (h in 1 : g) {
          alleles <- vector("character",2)
          alleles[1] <- paste0(LETTERS[h],GENO[k,h+w])
          alleles[2]<-paste0(LETTERS[h],GENO[k,h+1+w])
          w<-w+1
          alleles_list[[h]] <- alleles
        }
        alleles_list <- expand.grid(alleles_list)
        alleles_list<-apply(alleles_list,1,paste0,collapse="")
        for (h in 1: length(alleles_list)){
          P<-which(rownames(index)==alleles_list[h])
          RECOMB[P,k]<-RECOMB[P,k]+1/p
        }
      }
    }
  }
  
  
  MU<-array(NA,dim=c(l,l))
  type<-rownames(index)
  dimnames(MU)<-list(genotype_i=type,genotype_j=type)
  
  for(w in 1:l){
    for(ww in 1:l){
      p<-c()
      for (h in 1:g){
        if(substr(type[w],2*h,2*h)==substr(type[ww],2*h,2*h)){
          p<-c(p,1-mu[h])
        }
        else {
          p<-c(p,mu[h])
        }
      }
      MU[w,ww]<-prod(p)
    }
  }
  
  PROBA<-array(0,dim=c(l,m))
  dimnames(PROBA)<-list(type=rownames(index),genotype_parent=rownames(GENO))
  
  for (k in 1:m){
    type<-c()
    p<-c()
    for(kk in 1:l){
      if(RECOMB[kk,k]==0)next
      else{
        type<-c(type,kk)
        p<-c(p,RECOMB[kk,k])
      }
    }
    p_i<-c()
    for( i in 1:length(type)){
      p_i<-p[i]*MU[type[i],]
      PROBA[,k]<-PROBA[,k]+c(p_i)
    }
  }
  return(PROBA)
}

