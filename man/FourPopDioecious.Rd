\name{FourPopDioecious}
\alias{FourPopDioecious}
\docType{data}
\title{
Example datasets for a four-deme dioecious metapopulation.
}
\description{
Datasets of survival, fecundity, mutation, dispersal and initial population size for a four-deme metapopulation. Dioecious species.
}
\usage{data(FourPopDioecious)}
\format{
This dataset contains the objects \code{N1_M, N1_F, sigma_M, sigma_F, phi_M, phi_F, mu, kappa0, delta} and \code{T_max}.\cr 
\cr 
\code{N1_M} and \code{N1_F} are a 3-D arrays giving the number of male and female individuals, respectively, of each genotype in each deme and each age-class at the beginning of the simulation.  See \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description.\cr
\cr
\code{sigma_M}, \code{sigma_F}, \code{phi_M} and \code{phi_F} give the male survival probability, female survival probability, male fecundity and female fecundity, respectively.  See \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description.\cr
\cr
\code{mu} is is the mutation matrix. See \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description.\cr
\cr
\code{delta} is the dispersal matrix. See \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description.\cr
 \cr
\code{kappa0} is the carrying capacity. See the function \code{\link[MetaPopGen:sim.metapopgen.dioecious]{sim.metapopgen.dioecious}} for further details on this parameter.\cr
\cr
\code{T_max} is the number of years of simulation.
}

\details{
The \code{FourPopDioecious} dataset has one locus with 2 alleles (therefore, 3 genotypes), 4 demes and 6 age-classes. Survival and fecundities do not differ among demes, years or genotypes, but change with age. Survival probabilities differ with sex.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(FourPopDioecious)

# Run the simulation
N <- sim.metapopgen.dioecious(input.type="array",N1_M = N1_M, N1_F=N1_F, sigma_M=sigma_M, sigma_F=sigma_F, phi_F=phi_F, phi_M=phi_M, mu=mu, delta=delta, kappa0=kappa0, T_max=T_max)
}
\keyword{datasets}
