\name{MetaPopGen-package}
\alias{MetaPopGen-package}
\alias{MetaPopGen}
\docType{package}
\title{
Simulation of metapopulation genetics
}
\description{
Functions and code to simulate metapopulation genetics. Includes genetic drift, mutation, selection and migration. The metapopulation has age-structure. Migration is limited to the first age class (like seeds in plants and larvae in marine organisms). Functions are available to simulate one-locus genetics or for multi-locus simulations.
}
\details{
\tabular{ll}{
Package: \tab MetaPopGen\cr
Type: \tab Package\cr
Version: \tab 0.0.8\cr
Date: \tab 2017-12-01\cr
License: \tab GPL-2\cr
}
This package simulates metapopulation genetics accounting for genetic drift, migration, selection and mutation. 
}
\author{
Marco Andrello <marco.andrello@gmail.com>, Christelle Noirot
}
\references{
None
}

\keyword{ package }
\seealso{
The main functions for simulating one-locus metapopulation genetics are \code{\link[MetaPopGen:sim.metapopgen.monoecious]{sim.metapopgen.monoecious}} and \code{\link[MetaPopGen:sim.metapopgen.dioecious]{sim.metapopgen.dioecious}}. There are several example datasets, such as \code{\link[MetaPopGen:FourPopMonoecious]{FourPopMonoecious}}.
For multilocus simulations, see functions \code{\link[MetaPopGen:sim.metapopgen.monoecious.multilocus]{sim.metapopgen.monoecious.multilocus}} and \code{\link[MetaPopGen:sim.metapopgen.dioecious.multilocus]{sim.metapopgen.dioecious.multilocus}}
}
\examples{

}
