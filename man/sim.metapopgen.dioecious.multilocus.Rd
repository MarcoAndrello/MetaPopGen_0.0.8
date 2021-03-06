\name{sim.metapopgen.dioecious.multilocus}
\alias{sim.metapopgen.dioecious.multilocus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Main function to simulate metapopulation genetics (dioecious species).with several loci. Check the aid of sim.metapopgen.dioecious for more information. 
}
\description{
  Simulates the metapopulation genetics. The function will calculate automatically from its arguments the number of gamete type \code{l}, the number of genotypes \code{m}, the number of demes \code{n} and the number of age-classes \code{z}. It is therefore important that all the arguments have the correct dimensions and size.
}
\usage{
  sim.metapopgen.dioecious(input.type,demographic.data,
                           N1_M,N1_F,sigma_M,sigma_F,
                           phi_M,phi_F,List_gene,List_mu,r, 
                           delta,recr.dd="settlers",
                           kappa0,T_max,save.res=F,
                           save.res.T=seq(1,T_max),verbose=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.type}{
    Character vector defining whether the input data are giving as arrays (\code{"array"}) or as a data frame (\code{"data.frame"})
  }
  \item{demographic.data}{
    If the input is given as a data.frame (\code{input.type = "data.frame"}), a data frame giving 1) survival probabilities and fecundities for each genotype, deme, age-class, sex and time and 2) number of individuals for each genotype, deme, sex and age-class at the beginning of the simulation
  }
  \item{N1_M}{
    Number of male individuals at the beginning of the simulation. 3-D array giving the number of individuals of each genotype in each deme and each age-class at the beginning of the simulation. The first dimension is the genotype, the second dimension is the deme and the third dimension is the age-class. To be used when \code{input.type = "array"}.
  }
  \item{N1_F}{
    Number of female individuals at the beginning of the simulation. The same as \code{N1_M}.
  }
  \item{sigma_M}{
    Survival probabilities for males. 4-D array of length c(m,n,z,T_max). \code{sigma[j,i,x,t]} is the survival probability for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{sigma_M} can be a 3-D array of length c(m,n,z): this will be interpreted as male survival probabilities constant in time. All the elements of \code{sigma} are bounded between [0,1]. To be used when \code{input.type = "array"}.
  }
  \item{sigma_F}{
    Survival probabilities for females. The same as \code{sigma_M}.
  }
  \item{phi_F}{
    Mean female fecundity. 4-D array of length c(m,n,z,T_max). \code{phi_F[j,i,x,t]} is the female fecundity for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{phi_F} can be 
    a 3-D array of length c(m,n,z): this will be interpreted as female fecundities constant in time. The female fecundity is the per capita number of female gametes (eggs for animals, embrio sacs for plants). All the elements of \code{phi_F} must be natural numbers. To be used when \code{input.type = "array"}.
    Variation in fecundity among individuals follows a Poisson distribution.
  }
  \item{phi_M}{
    Mean male fecundity. 4-D array of length c(m,n,z,T_max). \code{phi_M[j,i,x,t]} is the male fecundity for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. Alternatively, \code{phi_M} can be 
    a 3-D array of length c(m,n,z): this will be interpreted as male fecundities constant in time. The male fecundity is the per capita number of male gametes (sperms for animals, pollen grains for plants). All the elements of \code{phi_F} must be natural numbers. To be used when \code{input.type = "array"}.
    Variation in fecundity among individuals follows a Poisson distribution.
  }
  \item{List_gene}{
    List with the number of allele for each loci.
  }
  \item{List_mu}{
   List with the mutation rate of each loci. 
  }
  \item{r}{
  Parameter defining the recombination probability in simulations with two loci.      \code{r} can take values between 0 and 0.5, with \code{r}=0.5 (default) for two      unlinked loci and  0<\code{r}<0.5 for linked loci. With more than two loci, loci are   considered unlinked. 
  }
  \item{delta}{
    Dispersal probabilities. Square matrix of size \code{n*n}. \code{delta[i,j]} is the dispersal probability from deme \code{j} to deme \code{i}. Alternatively, \code{delta} can be a 2-D array of length c(n,n): this will be interpreted as dispersal constant in time.
  }
  \item{recr.dd}{
    Defines how settler survival probabilities depend on density. If "settlers" (default), settler survival probabilities depends on the number of settlers, through the function \code{\link[MetaPopGen:settler.survival]{settler.survival}}. If "adults", settler survival probabilities depends on the number of adults (all individuals already present in the deme): only up to \code{kappa0[i,t]} (see next argument) individuals are allowed in a deme.
  }
  \item{kappa0}{
    2-D array of length c(n,T_max). Carrying capacity used to calculate the density-dependent settler survival probabilities. If the argument \code{recr.dd = "settlers"}, \code{kappa0[i,t]} is the settler carrying capacity for deme i at time t (maximum number of settlers allowed to recruit), and is used in the function \code{\link[MetaPopGen:settler.survival]{settler.survival}}. If \code{recr.dd = "adults"}, \code{kappa0[i,t]} is the maximum number of individuals allowed for deme i at time t. \code{kappa0} can also be a 1-D array of length n: this will be interpreted as a carrying capacity constant in time.
  }
  \item{T_max}{
    Number of years of simulations.
  }
  \item{save.res}{
    If FALSE (default), the result of the simulation is stored in memory and returned to the user. If TRUE, the function does not return anything, but save the result in different files on disk. See "value" below.
  }
  \item{save.res.T}{
    Vector. If save.res=TRUE, this argument is used to define the time steps for which the state of the metapopulation will be saved on disk. It is used to save disk space in the case of long simulations with many demes or genotypes. It deafults to \code{seq(1,T_max)}, that is, the results are saved every time step.
  }
  \item{verbose}{
    Boolean. If TRUE, the function will print much information, mainly used for debugging purposes. Default to FALSE.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  If \code{save.res} is FALSE: a named list with two elements, \code{N_M} and \code{N_F}, both 4-D arrays (ex: \code{N_M[k,i,x,t] N_F[k,i,x,t]}) of length c(m,n,z,tmax) giving the number of male and female individuals of genotype \code{k}, age \code{x} in deme \code{i} at time \code{t}. \cr
  If \code{save.res} is TRUE: the function does not return anything. Results are saved in files on disk. The function will create a folder named with the date and time of your simulation and create one RData file per year of simulation (but see argument \code{save.res.T} for defining the time steps for which the results are saved). Each file will contain the variables N_M and N_F, both 3-D arrays (\code{N_M[k,i,x] N_F[k,i,x]}) of length c(m,n,z) giving the number of male and fenale individuals of genotype \code{k}, age \code{x} in deme \code{i} at that time.
}
\references{
}
\author{
}
\note{
}
\seealso{ For more information, have a look in the tutoriel 3 that you can find on: https://sites.google.com/site/marcoandrello/metapopgen 
}
\examples{
  data("MultiDioecious")
  N<-sim.metapopgen.dioecious.multilocus(input.type = "array",N1_M=N1_M,N1_F=N1_F,sigma_M = sigma_M,sigma_F = sigma_F,phi_M = phi_M,phi_F = phi_F,List_gene = List_gene,List_mu=List_mu,delta=delta,kappa0=kappa0,T_max=5)
  
  # This will save results on disk:
  sim.metapopgen.dioecious.multilocus(input.type = "array",N1_M=N1_M,N1_F=N1_F,sigma_M = sigma_M,sigma_F = sigma_F,phi_M = phi_M,phi_F = phi_F,List_gene = List_gene,List_mu=List=mu, delta=delta,kappa0=kappa0,T_max=5,save.res=T)
  
  #with a Data.frame
  data("DataFrameDioeciousMulti")
  sim.metapopgen.dioecious.multilocus(input.type = "data.frame",demographic.data = demographic.data,List_gene = List_gene,List_mu=List_mu, delta=delta,kappa0=kappa0,T_max=5)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
