\name{FourPopMonoecious}
\alias{FourPopMonoecious}
\docType{data}
\title{
Example datasets for a four-deme monoecious metapopulation
}
\description{
Datasets of survival, fecundity, mutation, dispersal, initial population size and number of years of simulation for a four-deme monoecious metapopulation.
}
\usage{data(FourPopMonoecious)
data(FourPopMonoecious.selection)}
\format{
This dataset contains the objects \code{N1, sigma, phi_F, phi_M, mu, kappa0, delta} and \code{T_max}.\cr 
\cr 
\code{N1} is a 3-D array giving the number of individuals of each genotype in each deme and each age-class at the beginning of the simulation. The first dimension is the genotype, the second dimension is the deme and the third dimension is the age-class, so \code{N1[j,i,x]} is the number of individuals of genotype \code{j},in deme \code{i}, of age \code{x} at the beginning of the simulation. It is important to understand the order of genotypes; it will be easier to explain through some examples. If there are two alleles, say A1 and A2, there will be three genotypes; the order in which they are stored in N1 is as follows: A1A1, A1A2 and A2A2. If there are three alleles (A1, A2 and A3), the order is: A1A1, A1A2, A1A3, A2A2, A2A3, A3A3. With four alleles, the order is: A1A1, A1A2, A1A3, A1A4, A2A2, A2A3, A2A4, A3A3, A3A4, A4A4. So with four alleles \code{N[6,1,3]} gives the number of A2A4 individuals in deme 1 that are of age 3. \cr
\cr
\code{sigma}, \code{phi_F} and \code{phi_M} give the survival probability, female fecundity and male fecundity, respectively. The female fecundity is the per capita number of female gametes (eggs for animals, embrio sacs for plants) and the male fecundity is the per capita number of male gametes (sperms for animals, pollen grains for plants). All the elements of \code{phi_F} and \code{phi_M} must be natural numbers. All the elements of \code{sigma} must be rational numbers comprised between [0,1]. \code{sigma}, \code{phi_F} and \code{phi_M} are 4-D arrays and share the same structure: \code{sigma[j,i,x,t]} is is the survival probability for individuals of genotype \code{j}, age \code{x} in deme \code{i} at time \code{t}. The order of genotypes is of course the same as that in \code{N1}. \cr
\cr
\code{mu} is a square matrix of size \code{l*l}, with \code{l} = number of alleles. \code{mu[i,j]} is the probability that allele \code{j} mutates into allele \code{i}. The elements of the diagonal of \code{mu} are the probabilities that the allele does not mutate, i.e. it must be \code{colSums(mu)} = 1. It is assumed that the locus follows a K-allele mutation model (the number of alleles is fixed). This is appropriate for certain type of markers, such as single-nucleotide-polymorphisms (SNP). \cr \cr
\code{kappa0} is the carrying capacity. See the function \code{\link[MetaPopGen:sim.metapopgen.monoecious]{sim.metapopgen.monoecious}} for further details on this parameter. \cr \cr
\code{delta} is a square matrix of size \code{n*n}, with \code{n} = number of demes. \code{delta[i,j]} is the dispersal probability from deme \code{j} to deme \code{i}. Dispersal probabilities must therefore be obtained before running the model. For marine organisms, they can be obtained through biophysical simulations of larval dispersal.\cr
\cr
\code{T_max} Number of years of simulation.\cr
\cr
The dependency of survival and fecundity on genotype means that one can model selection through its effect on the survival and/or fecundity of a given genotype. Setting the same values for all genotype is equivalent to simulate a neutral locus. Changing values with genotype mimicks the effects of selection on survival and/or fecundity. The \code{FourPop.selection} dataset is the same as \code{FourPop} except that, in deme number 2, the first genotype (A1A1) has survival = 0  and the second, third and fourth genotypes (A1A2, A1A3 and A1A4) have survival = 0.18 (this corresponds to a selection coefficient s = 1 and dominance coefficient h = 0.7).

The dependency on deme and time means that the user is free to explore virtually any scenario of time-varying selection, spatial environmental variation, environmental gradients, etc. \cr
}

\details{
The \code{FourPopMonoecious} dataset has one locus with 4 alleles (therefore, 10 genotypes), 4 demes and 6 age-classes. The locus is fixed for a different allele in different demes. Survival and fecundities do not differ among demes, years or genotypes, but change with age. Survival increases from 0.4 at age 1 to 0.6 at ages 3-6. Fecundity is zero from age 1 to age 4, is moderate at age 5 and is maximal at age 6.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(FourPopMonoecious)
# Shows the survival and fecundity functions
plot(sigma[1,1,,1],xlab="Age",ylab="Survival probability")
plot(phi_F[1,1,,1],xlab="Age",ylab="Female fecundity")

# Shows the genotypic composition of demes
par(xpd=NA)
colours <- rainbow(10)
barplot(N1[,,1],names.arg=c("1","2","3","4"),xlab="Deme",ylab="Number of age-1 individuals",col=colours)
lg.txt <- c("A1A1", "A1A2", "A1A3", "A1A4", "A2A2", "A2A3", "A2A4", "A3A3", "A3A4", "A4A4")
legend(-1,1500,lg.txt,fill=colours,ncol=5,bty="n")

# Shows the dispersal matrix
dim(delta)
require(lattice)
levelplot(delta[,,1])

# Shows the mutation matrix
dim(mu)
levelplot(mu)

# Run the simulation
N <- sim.metapopgen.monoecious(input.type="array",N1=N1, sigma=sigma, phi_F=phi_F, phi_M=phi_M, mu=mu, delta=delta, kappa0=kappa0)

# Shows the genotypic composition of demes at the end of the simulation
# Compare with that at the beginning!
colours <- rainbow(10)
barplot(N[,,1,dim(N)[4]],names.arg=c("1","2","3","4"),xlab="Deme",ylab="Number of age-1 individuals",col=colours)
lg.txt <- c("A1A1", "A1A2", "A1A3", "A1A4", "A2A2", "A2A3", "A2A4", "A3A3", "A3A4", "A4A4")
legend(-1,1500,lg.txt,fill=colours,ncol=5,bty="n")

# Simulation with selection against the first allele in the second deme
data(FourPopMonoecious.selection)
N <- sim.metapopgen.monoecious(input.type="array",N1=N1, sigma=sigma, phi_F=phi_F, phi_M=phi_M, mu=mu, delta=delta, kappa0=kappa0)
par(xpd=NA)
colours <- rainbow(10)
barplot(N[,,2,dim(N)[4]],names.arg=c("1","2","3","4"),xlab="Deme",ylab="Number of age-2 individuals",col=colours)
lg.txt <- c("A1A1", "A1A2", "A1A3", "A1A4", "A2A2", "A2A3", "A2A4", "A3A3", "A3A4", "A4A4")
legend(0.2,600,lg.txt,fill=colours,ncol=5,bty="n")



}
\keyword{datasets}
