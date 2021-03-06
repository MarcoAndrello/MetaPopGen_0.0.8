# sim.metapopgen.coral.multilocus
# Marco Andrello
# 07/09/2018

# This function is based on sim.metapopgen.monoecious.multilocus It adapts the life cycle for coral colonies (individuals)
# There are five classes: 0+, 1+, 2+, 3+ and dead.
# Individuals can persist in the 3+ class. Individuals in the "dead" class decay at a constant rate
# Recruitment is dependent on the total number of individuals in the deme, including dead individuals.


sim.metapopgen.coral.multilocus <- function(input.type, demographic.data, N1,
                                                 sigma, phi_F, phi_M, List_gene, mu,
                                                 r=0.5, delta,
                                                 kappa0, T_max, save.res=F,
                                                 save.res.T=seq(1:T_max), verbose=F) {
  
  ##########################################################################
  
  # Initial definitions
  
  ##########################################################################
  
  
  if (input.type=="data.frame") {
    
    print("Input type = data.frame")
    
    a <- metapopgen.input.convert.monoecious(demographic.data)
    N1    <- a[[1]]
    sigma <- a[[2]]
    phi_M <- a[[3]]
    phi_F <- a[[4]]
    rm(a)
    
    
  } else {
    if (input.type == "array") {
      
      print("Input type = array")
      
    } else {
      stop("Unknown value for argument input.type. It must be either data.frame, array or txt")
      
    }
  }
  
  
  
  # Define basic variables
  
  m <- dim(N1)[1]                 # Number of genotypes
  l <- prod(List_gene)            # Number of multilocus gametes 
  n <- dim(N1)[2]                 # Number of demes
  z <- dim(N1)[3]                 # Number of age-classes
  
  
  # Define genotype indices matrix, giving the index for each genotype - Can we put it inside Create.probability.matrix ??
  index <- genotype.index.multilocus(List_gene)
  
  print ("creating gamete production probability matrix between parental genotypes and gametes...")
  Proba <- create.probability.matrix(index,List_gene,r,mu)
  print("...done")
  
  ##########################################################################
  # Check if input data are time-dependent or not
  ##########################################################################
  
  # Survival
  if (is.na(dim(sigma)[4])) {
    sigma <- array(rep(sigma,T_max),c(m,n,z,T_max))
  }
  
  # Female fecundity
  if (is.na(dim(phi_F)[4])) {
    phi_F <- array(rep(phi_F,T_max),c(m,n,z,T_max))
  }
  
  # Male fecundity
  if (is.na(dim(phi_M)[4])) {
    phi_M <- array(rep(phi_M,T_max),c(m,n,z,T_max))
  }
  
  # Dispersal
  if (is.na(dim(delta)[3])) {
      is.delta.constant <- T
  } else {
      is.delta.constant <- F
  }
  
  # Carrying capacity
  if (is.vector(kappa0)) {
    kappa0 <- array(rep(kappa0,T_max),c(n,T_max))
  }
  
  
  ##########################################################################
  # Initialize state variables
  ##########################################################################
  
  print("Initializing variables...")
  if (save.res){
    N <- N1
    rm(N1)    
  } else {
    N       <- array(NA,dim=c(m,n,z,T_max))
    N[,,,1] <- N1
    rm(N1)
    L       <- array(NA,dim=c(m,n,T_max))
    S       <- array(0,dim=c(m,n,T_max))
  }
  
  
  
  ##########################################################################
  # Define functions
  ##########################################################################
  
  # Survival
  surv <-
    function(sigma,N) {
      rbinom(1,N,sigma)
    }
  
  # Type gamete define the number of each type of gamete given by a genotype
  Type_gamete<-function(fec,Proba){
    
    l=length(Proba[,1]) 
    m=length(Proba[1,])
    
    Gamete<-array(0,dim=c(l,m))
    for (k in 1:m){
      err<- try (as.vector(rmultinom(1,fec,Proba[,k])),silent=T) #meiose
      if (class(err)=="try-error") {
        prob.mvr <- fec[k] * Proba[,k]				# Vector of means of the multivariate normal distribution
        var.mvr <- fec[k]* Proba[,k] * (1-Proba[,k])	# Vector of variances of the multivariate normal distribution
        sigma.mvr <- diag(var.mvr, l)			# Variance-covariance matrix of the multivariate normal distribution
        
        for (i.mvr in 1 : l){
          for (j.mvr in 1 : l) {
            if (i.mvr == j.mvr) next
            sigma.mvr[i.mvr,j.mvr] <- -fec[k] * Proba[i.mvr,k] * Proba[j.mvr,k]
          }
        }
        err<-as.vector(round(loinorm(1,prob.mvr,sigma.mvr)))
      }
      Gamete[,k]<-err
    }
    type<-array(0,l)
    for (j in 1:l){
      type[j]<-sum(Gamete[j,])
    }
    return(type)
  }
  
  # Reproduction multilocus 
  repr <-
    function(Nprime,phi_F,phi_M,l,m,Proba) {
      
      # Calculate number of female gametes for each gametype
      fecx<- array(0,dim=c(m,z))	# Number of female gametes produced by all the individuals of each genotype in each age class
      fec<- array(0,dim=m)			# Number of female gametes produced by all the individuals of each genotype
      
      for (k in 1 : m) {
        for (x in 1 : z) {
          fecx[k,x] <- sum(as.numeric(rpois(Nprime[k,x],phi_F[k,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
        }
        fec[k] <- sum(fecx[k,])
      }
      
      G_F <-Type_gamete(fec,Proba)
      
      
      # Calculate number of male gametes for each gametype
      
      fecx <- array(0,dim=c(m,z))	# Number of male gametes produced by all the individuals of each genotype in each age class
      fec <- array(0,dim=m)			# Number of male gametes produced by all the individuals of each genotype
      for (k in 1 : m) {
        for (x in 1 : z) {
          fecx[k,x] <- sum(as.numeric(rpois(Nprime[k,x],phi_M[k,x])))	# This is the contribution of variation in reproductive success among individuals to genetic drift
        }
        fec[k] <- sum(fecx[k,])
      }

      G_M <- Type_gamete(fec,Proba)
      
      # Union of gametes to form zygotes
      if (sum(G_F) <= sum(G_M)) {
        
        mat_geno<- array(0,dim=c(l,l))
        Gprime_M<- G_M
        for (j in 1 : l) {
          in_dist<- Gprime_M 
          odds<- array(1,dim=l)
          ndraws<- G_F[j]
          # If the mutlivariate hypergeometric does not work, use the multinomial
          # If the multinomial does not work, use the multivariate normal
          err1<- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
          if (class(err1)=="try-error") {
            
            err2<- try(as.numeric(rmultinom(1,ndraws,in_dist)),silent=T)
            if (class(err2)=="try-error") {
              
              # Use multivariate normal
              prob <- in_dist/sum(in_dist)
              mu.mvr <- ndraws * prob  			               # Vector of means of the multivariate normal distribution
              var.mvr <- ndraws * prob * (1-prob)	       # Vector of variances of the multivariate normal distribution
              sigma.mvr <- diag(var.mvr, l)			             # Variance-covariance matrix of the multivariate normal distribution
              
              for (i.mvr in 1 : l){
                for (j.mvr in 1 : l) {
                  if (i.mvr == j.mvr) next
                  sigma.mvr[i.mvr,j.mvr] <- -ndraws * prob[i.mvr] * prob[j.mvr]
                }
              }
              extr <- as.vector(round(loinorm(1,mu.mvr,sigma.mvr)))
              
            } else {
              extr<- err2                                     # Use multinomial
            }
            
          } else {
            extr <- err1
          }
          
          mat_geno[j,]<- extr
          Gprime_M<- Gprime_M - extr
        }
        mat_geno_l<- mat_geno
        mat_geno_l[upper.tri(mat_geno_l)]<- 0
        mat_geno_u<- mat_geno
        mat_geno_u[lower.tri(mat_geno_u,diag=T)]<- 0
        mat_geno_f<- mat_geno_l + t(mat_geno_u)
        L <- mat_geno_f[lower.tri(mat_geno_f,diag=T)]
        
      } else {
        
        
        mat_geno <- array(0,dim=c(l,l))
        Gprime_F <- G_F
        for (j in 1 : l) {
          in_dist <- Gprime_F 
          odds <- array(1,dim=l)
          ndraws <- G_M[j]
          # If the mutlivariate hypergeometric does not work, use the multinomial
          # If the multinomial does not work, use the multivariate normal
          
          err1 <- try(rMWNCHypergeo(1,in_dist,ndraws,odds),silent=T)
          
          if (class(err1)=="try-error") {
            
            err2 <- try(as.numeric(rmultinom(1,ndraws,in_dist)),silent=T)
            
            if (class(err2)=="try-error") {
              
              # Use multivariate normal
              prob <- in_dist/sum(in_dist)
              
              mu.mvr <- ndraws * prob  			               # Vector of means of the multivariate normal distribution
              var.mvr <- ndraws * prob * (1-prob)	       # Vector of variances of the multivariate normal distribution
              sigma.mvr <- diag(var.mvr, l)			             # Variance-covariance matrix of the multivariate normal distribution
              
              for (i.mvr in 1 : l){
                for (j.mvr in 1 : l) {
                  if (i.mvr == j.mvr) next
                  sigma.mvr[i.mvr,j.mvr] <- -ndraws * prob[i.mvr] * prob[j.mvr]
                }
              }
              
              extr <- as.vector(round(loinorm(1,mu.mvr,sigma.mvr)))
              
            } else {
              extr<- err2                                     # Use multinomial
            }
            
          } else {
            extr <- err1                                       # Use multivariate hypergeometric
          }
          
          mat_geno[j,]<- extr
          Gprime_F<- Gprime_F - extr
        }
        mat_geno_l<- mat_geno
        mat_geno_l[upper.tri(mat_geno_l)]<- 0
        mat_geno_u<- mat_geno
        mat_geno_u[lower.tri(mat_geno_u,diag=T)]<- 0
        mat_geno_f<- mat_geno_l + t(mat_geno_u)
        L<- mat_geno_f[lower.tri(mat_geno_f,diag=T)]
        
      }

      return(L)
      
    }
  
  
  
  
  # Dispersal
  disp <-
    function(L,delta){
      delta_lost <- max(0,1 - sum(delta)) # The maximum function is needed to avoid errors due to precision
      delta_add <- c(delta,delta_lost)
      S <- rmultinom(1,L,delta_add)
      return(S)
    }
  
  
  # Recruitment
  
  # For corals, survival of settlers is dependent on the number of individuals in the colony, including dead individuals, and the maximum carrying capacity
  # If there are already N individuals, then settler survival is dependent on (kappa0-N) = maximum number of individuals that can still be recruited
  # S        Number of settlers of all genotypes. Dimension: m. Corresponds to S[,i] in the main code.
  # N        Number of adults of all genotypes and age-classes. Dimensions: m*z. Corresponds to Nprime[,i,]
  # m        Number of genotypes
  # kappa0   Carrying capacity. Scalar. Corresponds to kappa0[i,t]
  # S[,i],Nprime[,i,],m,kappa0[i,t],recr.dd
  #S <- S[,i]
  #N <- Nprime[,i,]
  #kappa0 <- kappa0[i,t]
  recr <-
      function(S,N,m,kappa0) {
          Ntot <- sum(N) # Total number of individuals in the deme (summed over genotypes)
          Stot <- sum(S)
          Recr <- kappa0 - Ntot
          
          # Compute recruitment probability (settler survival)
          if (Recr <= 0){
              sigma0 <- 0
          } else {
              sigma0 <- Recr / Stot
          }
          
          if (sigma0 > 1) sigma0 <- 1
          
          # Use recruitment probability to calculate the number of recruits
          R <- array(0,dim=m)
          for (k in 1 : m) {
              R[k] <- rbinom(1,S[k],sigma0)
          }
          return(R)
      }
  
  
  
  ##########################################################################
  # Simulate metapopulation genetics
  ##########################################################################
  if (save.res){
    dir.res.name <- paste(getwd(),format(Sys.time(), "%Y-%b-%d-%H.%M.%S"),sep="/")
    dir.create(dir.res.name)
    
    if (1 %in% save.res.T) {
      file.name <- "N1.RData"
      save(N,file=paste(dir.res.name,file.name,sep="/"))
    }
  }
  
  
  print("Running simulation...")
  for (t in 1 : (T_max-1)) {
    
    print(t)
    
    
    # At each time-step, redefine variable Nprime
    # If save.res, redefine also larval and settlers numbers
    if (save.res) {
      Nprime  <- array(NA,dim=c(m,n,z))
      L       <- array(NA,dim=c(m,n))
      S       <- array(0,dim=c(m,n))
    } else {
      Nprime  <- array(NA,dim=c(m,n,z))
      #L       <- array(NA,dim=c(m,n,T_max))
      #S       <- array(0,dim=c(m,n,T_max))
    }
    
    
    ### Survival
    
    # If there is only one age-class, we must force the third dimension. (This should not be necessary for corals)
    if (length(dim(sigma))==2) dim(sigma)[3] <- 1
    
    if (verbose) print("Apply survival function")
    for (i in 1 : n) {
      for (x in 1 : z) {
        for (k in 1 : m) {
          
          if (save.res){
            Nprime[k,i,x] = surv(sigma[k,i,x,t],N[k,i,x])
          } else {
            Nprime[k,i,x] = surv(sigma[k,i,x,t],N[k,i,x,t])
          }
        }
      }
    }
    
    
    
    if (verbose) print("Apply reproduction function")
    # If there is only one age-class, we must force the third dimension. (This should not be necessary for corals)
    if (length(dim(phi_F))==2) dim(phi_F)[3] <- 1
    if (length(dim(phi_M))==2) dim(phi_M)[3] <- 1
    
    for (i in 1 : n) {
      if (save.res) {
        if (sum(Nprime[,i,])==0) { # To save computing time
          L[,i] = 0
          next
        } else {
          L[,i] <- repr(Nprime[,i,],phi_F[,i,,t],phi_M[,i,,t],l,m,Proba)
        }
      } else {
        if (sum(Nprime[,i,])==0) { # To save computing time
          L[,i,t] = 0
          next
        } else {
          L[,i,t] <- repr(Nprime[,i,],phi_F[,i,,t],phi_M [,i,,t],l,m,Proba)
        }
      }
    }
    
    
    if (verbose) print("Apply dispersal function")
    for (i in 1 : n) {
        if (verbose) cat("Deme",i,"\n")
        for (k in 1 : m) {
            if (is.delta.constant){ 
                delta.1 <- delta[,i]
            } else {
                delta.1 <- delta[,i,t]
            }
            if (save.res) {
                y = disp(L[k,i],delta.1)
                S[k,] <- S[k,] + y[1:n]       
            } else {
                y = disp(L[k,i,t],delta.1)
                S[k,,t] <- S[k,,t] + y[1:n]
            }
        }
    } 
    
    # Aging
    # Calculates the number of individuals dead in this time step.
    # dead[k,i,x] is the number of dead individuals for genotype, deme and stage.
    # It is the difference between N (number of individuals before survival phase) and Nprime (number of individuals after survival phase)
    if (save.res) {
        dead <- N - Nprime
    } else {
        dead <- N[,,,t] - Nprime
    }
    dead <- dead[,,1:4] # we do not count the difference between individuals of the fifth stage class (already dead): this difference would give the number of destroyed coral skeletons
    
    if (verbose) print("Aging")
    for (i in 1 : n){
        for (x in 1 : z) {
            
            if (x == 1) { # For x==1 (stage 0+), we set it to 0 and calculate it later, after calculating recruitment
                if (save.res) {
                    N[,i,x] <- 0
                } else {
                    N[,i,x,t+1] <- 0
                }
                next
            }
            
            if (x == 4) { # For x==4 (stage 3+), includes the 2+' and 3+'
                for (k in 1 : m) {
                    if (save.res) {
                        N[k,i,x] <- Nprime[k,i,x-1] + Nprime[k,i,x]
                    } else {
                        N[k,i,x,t+1] <- Nprime[k,i,x-1] + Nprime[k,i,x]
                    }
                }
                next
            }
            
            if (x == 5) { # For x==5 (stage "dead"), includes the dead' and all the new dead (sum on stage-classes for that genotype and that deme)
                for (k in 1 : m) {
                    if (save.res) {
                        N[k,i,x] <- Nprime[k,i,x] + sum(dead[k,i,])
                    } else {
                        N[k,i,x,t+1] <- Nprime[k,i,x] + sum(dead[k,i,])
                    }
                }
                next
            }
            
            # For x==2 or x==3 (stages 1+ and 2+), includes the 0+' or the 1+', respectively.
            for (k in 1 : m) {
                if (save.res) {
                    N[k,i,x] <- Nprime[k,i,x-1]
                } else {
                    N[k,i,x,t+1] <- Nprime[k,i,x-1]
                }
            }
        }
    }
    
    if (verbose) print("Apply recruitment function")
    # To calculate recruitment, we need to know how many (alive+dead) individuals are present in the deme
    # Loop on the demes
    for (i in 1 : n) {
        if (save.res) {
            # Apply the recruitment function to find the number of recruits from the number of settlers
            # Arguments passed to the recruitment function:
            # S[,i]             the vector of settler genotype numbers for this deme i
            # m                 the number of genotypes
            # kappa0[i,t]       the carrying capacity for this deme and this time
            # sum(Nprime[,i,]) the number of individuals already present in this deme in all stage-classes (including the dead individuals)
            # This sums all alive and all dead individuals remaining after aging (decay of dead)
            N[,i,1] <- recr(S[,i],sum(N[,i,]),m,kappa0[i,t])
            
        } else {
            
            # Apply the recruitment function to find the number of recruits from the number of settlers
            # Arguments passed to the recruitment function:
            # S[,i,t]           the vector of settler genotype numbers for this deme i and this time t
            # other parameters as above
            N[,i,1,t+1] <- recr(S[,i,t],sum(N[,i,]),m,kappa0[i,t])
            
        }
    }
    
    
    # Save results if save.res=T
    if (save.res){
      if ((t+1) %in% save.res.T) {
        file.name <- paste("N",(t+1),".RData",sep="")
        save(N,file=paste(dir.res.name,file.name,sep="/"))
      }
    }
    
    
  }
  print(T_max)
  print("...done")
  if (save.res==F)
    return(N)
}
