\name{FourPopDataFrameDioecious}
\alias{FourPopDataFrameDioecious}
\docType{data}
\title{
Example datasets for a four-deme dioecious metapopulation, with input given as a data frame.
}
\description{
Datasets of survival, fecundity, mutation, dispersal and initial population size for a four-deme metapopulation. 
}
\usage{data(FourPopDataFrameDioecious)}
\format{
This dataset contains the objects \code{demographic.data, mu, delta} and \code{T_max}.\cr 
\cr 
\code{demographic.data} is a data frame containing survival probabilities and fecundities for each genotype, deme, age-class and time and the initial number of individuals for time = 1. 
\cr \cr
\code{mu} is is the mutation matrix. See \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description.\cr
\cr
\code{delta} is the dispersal matrix. See \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}} for a detailed description.\cr
\cr
\code{kappa0} is the carrying capacity. See the function \code{\link[MetaPopGen:sim.metapopgen.dioecious]{sim.metapopgen.dioecious}} for further details on this parameter.\cr
\cr
\code{T_max} is the number of years of simulation.
}

\details{
The \code{FourPopDataFrameDioecious} dataset has one locus with 2 alleles (therefore, 3 genotypes), 4 demes and 6 age-classes. Survival and fecundities do not differ among demes, years or genotypes, but change with age.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(FourPopDataFrameDioecious)

# Run the simulation
N <- sim.metapopgen.dioecious(input.type="data.frame",demographic.data=demographic.data, mu=mu, delta=delta, kappa0=kappa0, T_max=T_max)
}
\keyword{datasets}
