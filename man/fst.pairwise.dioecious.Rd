\name{fst.pairwise.dioecious}
\alias{fst.pairwise.dioecious}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates Fst between two groups of dioecious individuals 
}
\description{
}
\usage{
fst.pairwise.dioecious(N,i,j,xi,xj,ti,tj,si,sj)
fst.pairwise.dioecious(folder.name,i,j,xi,xj,ti,tj,si,sj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
A 4-D array (ex: \code{N[k,i,x,t]}) of length c(m,n,z,T_max) giving the number of individuals of genotype \code{k}, age \code{x} in deme \code{i} at time \code{t}. Typically, the output of the function \code{\link[MetaPopGen:sim.metapopgen.dioecious]{sim.metapopgen.dioecious}}
}
  \item{folder.name}{
Character string giving the name of the folder where the results of \code{\link[MetaPopGen:sim.metapopgen.dioecious]{sim.metapopgen.dioecious}} are saved, when this function is used with the argument  \code{save.res} set to TRUE.
}
  \item{i, xi, ti}{
Identifiers of the first group of individuals: it specifies the individuals of deme \code{i}, age \code{xi} at time \code{ti}
}
  \item{i, xi, ti}{
Identifiers of the second group of individuals: as above. 
}
  \item{si, sj}{
Identifiers of the sex of the first and second group of individuals. Character string: \code{"M"} for males and \code{"F"} for females.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Pairwise genetic differentiation index Fst calculated as ratio of heterozygosities: \code{(Ht - Hs) / Ht}. \code{Ht} is the expected heterozygosity over the population formed by the two groups. \code{Hs} is the average expected heterozygosity over the two groups. The average is weighted by the number of individuals in each group: \code{ Hs = (n_i * H_S_i + n_j * H_S_j)/(n_i + n_j)} where \code{n_i} and \code{n_j} are the number of individuals in the first and second group, respectively, and \code{H_s_i} and \code{H_s_j} are the expected heterozygosities in the first and second group, respectively.
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

# Load the example data and simulates metapopulation genetics
data(FourPopDioecious)
N <- sim.metapopgen.dioecious(input.type="array", N1_M=N1_M, N1_F=N1_F, sigma_M=sigma_M, sigma_F=sigma_F, phi_F=phi_F, phi_M=phi_M, mu=mu, delta=delta, settler.surv=settler.surv)

# Calculates fst between various pairs of groups of individuals
fst.pairwise.dioecious(N,1,2,1,1,5,5,"M","M") # Between one-year old male individuals of deme 1 and 2 at time 5
fst.pairwise.dioecious(N,1,1,1,1,5,5,"M","F") # Between male and female one-year old individuals of deme 1 at time 5


# Example with results from sim.metapopgen.dioecious saved on disk
data(FourPopDioecious)
N <- sim.metapopgen.dioecious(input.type="array", N1_M=N1_M, N1_F=N1_F, sigma_M=sigma_M, sigma_F=sigma_F, phi_F=phi_F, phi_M=phi_M, mu=mu, delta=delta, settler.surv=settler.surv, save.res=T)
# Let us assume that results are stored in a folder named "2014-avr.-11-12.53"
fst.pairwise.dioecious("2014-avr.-11-12.53",1,1,1,1,5,5,"M","F") # Between male and female one-year old individuals of deme 1 at time 5
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
