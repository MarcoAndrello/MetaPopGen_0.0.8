
#Calcul allele and genotypes frenquencies in a population  

Freq.genotype<- function(N){
  m<-length(N)
  P<-array(NA,m)
  ntot<-sum(N)
  for(i in 1:m){
    P[i]<-N[i]/ntot
  }
  return(P)
}

Freq.allele<-function(N,List_gene,index){
  P<-array(NA,dim=sum(List_gene))
  GENO<-create.genotype.matrix(List_gene,index)
  ntot<-sum(N)
  k<-1
  for (kk in 1:length(List_gene)){
    for(kkk in 1:List_gene[kk]){
      pos<-Position(GENO,kk,kkk)
      P[k]<-0.5*sum(N[pos])/ntot
      k<-k+1
    }
  }
  return(P)
}


#gives the position in the Genotype matrix GENO, with the number of the gene and the number
#of the allele 

Position<-function(GENO,Gene,allele){
  i<-Gene
  pos<-c()
  for(k in 1:length(GENO[,1])){
    if(GENO[k,i]==allele){
      pos<-c(pos,k)
    }
    if (GENO[k,i+2]==allele){
      pos<-c(pos,k)
    }
  }
  return(pos)
}







