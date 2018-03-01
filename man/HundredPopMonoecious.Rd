\name{HundredPopMonoecious}
\alias{HundredPopMonoecious}
\docType{data}
\title{
Example datasets for a hundred-deme monoecious metapopulation
}
\description{
Datasets of survival, fecundity, mutation, dispersal and initial population size for a hundred-deme monoecious metapopulation.
}
\usage{data(HundredPop)}
\format{
This dataset contains the objects \code{N1, sigma, phi_F, phi_M, mu, delta, kappa0} and \code{T_max}, as in the dataset \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}}. The help page of the \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} dataset contains a detailed explanation for each object.\cr 
}
\details{
The \code{HundredPopMonoecious} dataset has one locus with 2 alleles (therefore, 3 genotypes), 100 demes and 5 age-classes. Initial demes are composed of 1000 individuals per age class, all of the same genotype. Survival and fecundities do not differ among demes, years or genotypes, but change with age. Survival increases from 0.4 at age 1 to 0.6 at ages 3-5. Fecundity is zero from age 1 to age 3, is moderate at age 4 and is maximal at age 5.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(HundredPopMonoecious)

# Shows the genotypic composition of demes
colours <- rainbow(3)
barplot(N1[,,1],xlab="Deme",ylab="Number of age-1 individuals",col=colours)
par(xpd=NA)
lg.txt <- c("A1A1", "A1A2", "A2A2")
legend(0.2,1200,lg.txt,fill=colours,ncol=3,bty="n")

# Run the simulation
N <- sim.metapopgen.monoecious(input.type="array",N1=N1, sigma=sigma, phi_F=phi_F, phi_M=phi_M, mu=mu, delta=delta, kappa0=kappa0, T_max=T_max)

# Shows the genotypic composition of demes at the end of the simulation
# Compare with that at the beginning!
colours <- rainbow(3)
barplot(N[,,1,dim(N)[4]],xlab="Deme",ylab="Number of age-1 individuals",col=colours)
lg.txt <- c("A1A1", "A1A2", "A2A2")
legend(0.2,1200,lg.txt,fill=colours,ncol=3,bty="n")

# Temporal evolution of fst. Pairwise comparison for deme 1 vs 2, 1 vs 50 and 1 vs 100
fst <- array(NA,dim=c(3,100))
vdeme <- c(2,50,100)
for (t in 1 : 100) {
for (i in 1 : 3){
deme <- vdeme[i]
fst[i,t] <- fst.pairwise.monoecious(N,1,deme,1,1,t,t)
}
}
plot(c(1:100),fst[1,],type="l",col="blue",xlab="Time",ylab="Pairwise Fst",ylim=c(0,0.4))
matplot(c(1:100),fst[2,],type="l",col="darkgreen",add=T)
matplot(c(1:100),fst[3,],type="l",col="red",add=T)
legend(30,0.4,c("Deme 1 - Deme 2","Deme 1 - Deme 50","Deme 1 - Deme 100"),col=c("blue","darkgreen","red"),lty=1)
}

\keyword{datasets}
