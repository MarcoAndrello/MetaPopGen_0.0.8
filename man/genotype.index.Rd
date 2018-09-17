\name{genotype.index}
\alias{genotype.index.multilocus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotype indices
}
\description{
Functions giving genotype indices as a function of the number of alleles
}
\usage{
genotype.index(l)
genotype.index.multilocus(List_gene)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{
Number of alleles
}
  \item{List_gene}{
Vector giving the number of alleles at each locus 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For \code{genotype.index}, a square matrix whose [i,j] element gives the index of the genotype formed by the i and j alleles. For \code{genotype.index.multilocus}, a square matrix whose [i,j] element gives the index of the multilocus genotype formed by the i and j gamete types.
}
\references{
}
\author{
Marco Andrello
}
\note{
}
\seealso{
}
\examples{
# Monolocus case:
# To know the index of the A5A9 genotype if there are 13 alleles:
y <- genotype.index(13)
y[5,9]

# Multilocus case:
genotype.index.multilocus(c(2,3))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
