\name{fst.pairwise.monoecious}
\alias{fst.pairwise.monoecious}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates Fst between two groups of monoecious individuals 
}
\description{
}
\usage{
fst.pairwise.monoecious(N,i,j,xi,xj,ti,tj)
fst.pairwise.monoecious(folder.name,i,j,xi,xj,ti,tj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
A 4-D array (ex: \code{N[k,i,x,t]}) of length c(m,n,z,T_max) giving the number of individuals of genotype \code{k}, age \code{x} in deme \code{i} at time \code{t}. Typically, the output of the function \code{\link[MetaPopGen:sim.metapopgen]{sim.metapopgen.monoecious}}
}
  \item{folder.name}{
Character string giving the name of the folder where the results of \code{\link[MetaPopGen:sim.metapopgen]{sim.metapopgen.monoecious}} are saved, when this function is used with the argument  \code{save.res} set to TRUE.
}
  \item{i, xi, ti}{
Identifiers of the first group of individuals: it specifies the individuals of deme \code{i}, age \code{xi} at time \code{ti}
}
  \item{i, xi, ti}{
Identifiers of the second group of individuals: as above. 
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
data(FourPopMonoecious)
N <- sim.metapopgen.monoecious(input.type="array",N1=N1, sigma=sigma, phi_F=phi_F, phi_M=phi_M, mu=mu, delta=delta, settler.surv=settler.surv)

# Calculates fst between various pairs of groups of individuals
fst.pairwise.monoecious(N,1,2,1,1,50,50) # Between one-year old individuals of deme 1 and 2 at time 50
fst.pairwise.monoecious(N,1,1,1,1,1,50) # Between one-year old individuals of deme 1 at time 50 and at time 1 (temporal fst)

# Temporal evolution of fst
fst <- array(NA,dim=c(3,50))
for (t in 1 : 50) {
for (deme in 1 : 3){
fst[deme,t] <- fst.pairwise.monoecious(N,1,(deme+1),1,1,t,t)
}
}
plot(c(1:50),fst[1,],type="l",col="blue",xlab="Year",ylab="Fst")
matplot(c(1:50),fst[2,],type="l",col="darkgreen",add=T)
matplot(c(1:50),fst[3,],type="l",col="red",add=T)
legend(30,1,c("1 vs 2","1 vs 3","1 vs 4"),col=c("blue","darkgreen","red"),lty=1)

# Example with results from sim.metapopgen.monoecious saved on disk
# Let us assume that results are stored in a folder named "2013-nov-13-12.04"
fst.pairwise.monoecious("2013-nov-13-12.04",1,1,1,1,1,50) # Between one-year old individuals of deme 1 at time 50 and at time 1 (temporal fst)


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
