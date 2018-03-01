\name{create.probability.matrix}
\alias{create.probability.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create probability matrix
}
\description{
Function giving the probability of each gamete type for each parental genotype.

}
\usage{
create.probability.matrix(genotype.index, List_gene, r=0.5, List_mu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype.index}{
 The output of the function genotype.index.multilocus, giving the index (order) of each genotype in the simulation. 
}
  \item{List_gene}{
A vector of length \code{g} (number of loci), giving the number of alleles at each locus. 
}
\item{r}{
 Parameter defining the recombination probability in simulations with two loci. \code{r} can take values between 0 and 0.5, with \code{r}=0.5 (default) for two unlinked loci and  0<\code{r}<0.5 for linked loci. With more than two loci, loci are considered unlinked. 
}
\item{List_mu}{
A vector giving the mutation probability of one allele into an other existing allele for each locus
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A square symmetric matrix whose [i,j] element gives the probability of having a gamete of type i knowing that the parent is of genotype j.
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
# With two loci with two alleles each; linked loci
# There are 4 possible gametes: A1B1, A2B1, A1B2, A2B2
# There are 10 possible genotypes:
# A1B1/A1B1, A1B1/A2B1, A1B1/A1B2, A1B1/A2B2,
#            A2B1/A2B1, A2B1/A1B2, A2B1/A2B2,
#                       A1B2/A1B2, A1B2/A2B2,
#                                  A2B2/A2B2
List_gene <- c(2,2)
genotype.index <- genotype.index.multilocus(List_gene)
List_mu <- c(1e-06,1e-06)
r <- 0.2 
create.probability.matrix(genotype.index, List_gene, r, List_mu)

# With 4 loci of respectively 2,3,2 and 5 alleles
List_gene <- c(2,3,2,5)
genotype.index <- genotype.index.multilocus(List_gene)
List_mu <- c(1e-06,1e-06)
create.probability.matrix(genotype.index, List_gene, List_mu=List_mu)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
