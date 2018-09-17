\name{pop.loci}
\alias{pop.loci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Genotype frequencies per locus
}
\description{
}
\usage{
    pop.loci(Proba,List_gene,loci,N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Proba}{
    Matrix giving the probability that a given genotype (column) produces a given gamete type (row). It is calculated with the function \code{\link[MetaPopGen:create.probability.matrix]{create.probability.matrix}}. 
    }
    \item{List_gene}{
        Vector of numbers of alleles per locus.
    }
    \item{loci}{
        Index of the locus for which the genotype frequencies are sought.
    }
    \item{N}{
        Array of 2 dimension giving the number of individuals of each genotype in each deme: the first dimension is the genotype and the second dimension is the deme.  
    }
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    A square matrix giving the number of individuals of each genotype (row) in each deme (column).
}
\references{
}
\author{Christelle Noirot}
\note{
}
\seealso{
}
\examples{
   
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
