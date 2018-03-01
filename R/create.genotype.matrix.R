#GENO matrice gives 
create.genotype.matrix<-function(List_gene,index){

  g<-length(List_gene)
  m<-max(index)
  l<-prod(List_gene)
  
  GENO<-array(NA,dim=c(m,2*g))
  
  names<-c()
  for(i in 1:l){
    for(ii in 1:l){
      if(is.na(GENO[index[i,ii],1])){
        w<-0
        for (h in 1:g){
          a<-c(substr(rownames(index)[i],2*h,2*h),substr(colnames(index)[ii],2*h,2*h))
          GENO[index[i,ii],c(h+w,2*h)]<-as.numeric(a)
          w<-w+1
        }
        names<-c(names,paste0(rownames(index)[i],"/",colnames(index)[ii]))
      }
    }
  }
  dimnames(GENO)=list(names,rep(1:2,g))
  return(GENO)
}
