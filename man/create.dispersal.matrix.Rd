\name{create.probability.matrix}
\alias{create.probability.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create probability matrix
}
\description{
Function to create a matrix givin the probability to have a type of gamete knowing the genotype of the individual before the meiose.

}
\usage{
create.probability.matrix(index,List_gene,r=0.5,List_mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{index}{
 is a matrix calculated with the function genotype.index.multilocus, indicates the row of each genotype in the simulation. 
}
  \item{List_gene}{
is a list giving the number of loci and allele for each loci. The length of the list correspond to the number of loci, and each number indecates the number of alleles for this loci. 
}
\item{r}{
 defines the recombinasion probability if the simulation use only two loci. If \code{r}=0.5, loci are independant, if 0<\code{r}<0.5 loci are linked. For three and more loci in the simulation, loci are concidered independant. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A square symmetric matrix whose [i,j] element gives the probability of dispersal from deme j to deme i.
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
lon <- c(12.3, 13.3, 14.2, 16.9, 16.1, 17.5, 18.1, 19.7, 20.3)
lat <- c(45.4, 46.5, 48.9, 47.1, 45.8, 46.1, 47.0, 46.8, 47.2)
x <- cbind(lon,lat)
create.dispersal.matrix(x,0.01,longlat=T)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
