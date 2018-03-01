pop.loci<-function(Proba,List_gene,loci,N){
  m <- length(Proba[1,])
  n<-length(N[1,])
  mLoc<-List_gene[loci]*(List_gene[loci]+1)/2
  
  name_Loc<-c()
  for (k in 1:m){
    name_Loc<-c(name_Loc,paste0(substr(colnames(Proba)[k],loci+loci-1,loci+loci),substr(colnames(Proba)[k],2*length(List_gene)+loci+loci,2*length(List_gene)+loci+loci+1)))
    
  }
  k<-1
  while(k<=length(name_Loc)){
    kk<-1
    r<-c()
  
    while(kk<=length(name_Loc)){
      if(kk!=k){
        if (name_Loc[k]==name_Loc[kk]){
          r<-c(r,kk)
        }
        if (paste0(substr(name_Loc[k],1,2),substr(name_Loc[k],3,4))==paste0(substr(name_Loc[kk],3,4),substr(name_Loc[kk],1,2))){
          rep<-c(rep,paste0(substr(name_Loc[kk],3,4),substr(name_Loc[kk],1,2))) 
        }
      }
      kk<-kk+1
    }
    k<-k+1
    if (length(r)!= 0){
      name_Loc<-name_Loc[-r]
    }
  }
  
  NLoc<-array(0,dim=c(length(name_Loc),n))
  dimnames(NLoc)<-list(name_Loc,c(1:n))
 
  
  
  for (k in 1:length(name_Loc)){
    L<-which(paste0(substr(colnames(Proba),loci+loci-1,loci+loci),substr(colnames(Proba),2*length(List_gene)+loci+loci,2*length(List_gene)+loci+loci+1))==rownames(NLoc)[k] )
    for (l in 1:length(L)){
      NLoc[k,]<-NLoc[k,]+N[L[l],]
    }
  }
  
  k<-1
  while(k<=length(name_Loc)){
    rep<-c()
    kk<-1
    while(kk<=length(name_Loc)){
      if(kk!=k){
        if (paste0(substr(name_Loc[k],1,2),substr(name_Loc[k],3,4))==paste0(substr(name_Loc[kk],3,4),substr(name_Loc[kk],1,2))){
          
          NLoc[k,]<-NLoc[k,]+NLoc[kk,]
          NLoc<-NLoc[-kk,]
          name_Loc<-name_Loc[-kk]
        }
      }
      kk<-kk+1
    }
    k<-k+1
   }
  return(NLoc)
}
