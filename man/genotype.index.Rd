\name{genotype.index}
\alias{genotype.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper function to get genotype indices as a function of the number of alleles
}
\description{
}
\usage{
genotype.index(l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{
Number of alleles
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A square matrix whose [i,j] element gives the index of the genotype formed by the i and j allele
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
# If you Want to know the index of the A5A9 genotype if there are 13 alleles
y <- genotype.index(13)
y[5,9]
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
