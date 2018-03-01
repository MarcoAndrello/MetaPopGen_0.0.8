\name{FourPopDataFrameMonoecious}
\alias{FourPopDataFrameMonoecious}
\docType{data}
\title{
Example datasets for a four-deme monoecious metapopulation, with input given as a data frame.
}
\description{
Dataset of survival, fecundity, mutation, dispersal, initial population size and number of years of simulation for a four-deme metapopulation. 
}
\usage{data(FourPopDataFrameMonoecious)}
\format{
This dataset contains the objects \code{demographic.data, mu, kappa0, delta} and \code{T_max}.\cr 
\cr 
\code{demographic.data} is a data frame containing survival probabilities and fecundities for each genotype, deme, age-class and time and the initial number of individuals for time = 1. 
\cr \cr
\code{mu} is is the mutation matrix. See \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description.\cr
\cr
\code{delta} is the dispersal matrix. See \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description.\cr
 \cr
\code{kappa0} is the carrying capacity. See the function \code{\link[MetaPopGen:sim.metapopgen.monoecious]{sim.metapopgen.monoecious}} for further details on this parameter.\cr
\cr
\code{T_max} Number of years of simulation.\cr
} 

\details{
The \code{FourPopDataFrameMonoecious} dataset has one locus with 2 alleles (therefore, 3 genotypes), 4 demes and 6 age-classes. Survival and fecundities do not differ among demes, years or genotypes, but change with age. Survival increases from 0.4 at age 1 to 0.6 at ages 3-6. Fecundity is zero from age 1 to age 4, is moderate at age 5 and is maximal at age 6.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(FourPopDataFrameMonoecious)

# Run the simulation
N <- sim.metapopgen.monoecious(input.type="data.frame",demographic.data=demographic.data, mu=mu, delta=delta, kappa0=kappa0)
}
\keyword{datasets}
