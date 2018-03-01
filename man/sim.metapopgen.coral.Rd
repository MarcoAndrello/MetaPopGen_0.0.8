\name{sim.metapopgen.coral}
\alias{sim.metapopgen.coral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate metapopulation genetics of coral colonies
}
\description{
This function is built upon \code{\link[MetaPopGen:sim.metapopgen.monoecious]{sim.metapopgen.monoecious}} to simulate the genetics of coral metapopulations. Each coral colony is considered as an individual. \cr
There are five stage-classes, corresponding to age 0+, age 1+, age 2+, age 3+ and dead individuals. Individuals can persist in the 3+ class. Individuals in the "dead" class are coral skeletons. The "survival" probability of the dead class indicates whether a coral skeleton persists or is destroyed.
Recruitment is dependent on the total number of individuals in the deme, including dead individuals, because they occupy space. \cr
The parameter \code{kappa0} is the maximum number of individuals in the deme (carrying capacity).
}
\usage{
sim.metapopgen.coral(N1, sigma, phi_F, phi_M, kappa0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N1}{
Number of individuals at the beginning of the simulation. 3-D array giving the number of individuals of each genotype in each deme and each of the five stage-classes (0+, 1+, 2+, 3+ and dead) at the beginning of the simulation. The first dimension is the genotype, the second dimension is the deme and the third dimension is the stage-class. To be used when \code{input.type = "array"}.
}
  \item{sigma, phi_F, phi_M}{
They have the same meaning as in \code{\link[MetaPopGen:sim.metapopgen.monoecious]{sim.metapopgen.monoecious}}, the only difference being the meaning of the fifth stage class (dead individuals). The "survival" probability of the dead class, \code{sigma[i,j,5]}, indicates whether a coral skeleton persists or is destroyed: for example, if \code{sigma[i,j,5] = 0.85}, coral skeletons are destroyed with a 0.15 probability each time-step. The fecundities of the fifth stage class (dead individuals) must be set to 0.
}
  \item{kappa0}{
  2-D array of length c(n,T_max). Carrying capacity used to calculate the density-dependent settler survival probabilities. \code{kappa0[i,t]} is the maximum number of individuals allowed for deme i at time t. If the total number of individuals Ntot(i,t) present in a deme (including dead individuals: coral skeletons) is larger than \code{kappa0[i,t]}, then recruitment probability is set to 0. Otherwise, recruitment probability is set to [\code{kappa0[i,t]} - Ntot(i,t)] / Stot(i,t), where Stot(i,t) is the total number of settlers in deme i at time t. \code{kappa0} can also be a 1-D array of length n: this will be interpreted as a carrying capacity constant in time.
}
  \item{...}{
All the other parameters to be passed to the function, as in \code{\link[MetaPopGen:sim.metapopgen.monoecious]{sim.metapopgen.monoecious}}
}
}
  
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The same as in \code{\link[MetaPopGen:sim.metapopgen.monoecious]{sim.metapopgen.monoecious}}.
}
\references{
}
\author{
}
\note{
}
\seealso{
}
\examples{
# Simple example with 2 demes, 2 alleles and 5 years of simulations
data(Dataset.coral)
sim.metapopgen.coral(input.type="array", N1= N1, sigma = sigma, phi_M = phi_M, phi_F = phi_F, mu = mu, delta = delta, kappa0 = kappa0, T_max = T_max, save.res = save.res)

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
