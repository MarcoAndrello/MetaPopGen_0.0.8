\name{pop.loci}
\alias{pop.loci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Funcion to separate the population in genotype as if there where only one loci concidered 
}
\description{
}
\usage{
    pop.loci(Proba,List_gene,loci,N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Proba}{
    Matrix with the names of all the genotypes. This array is calculated with the function {create.probability.matrix}
    }
    \item{List_gene}{
        List with the number of alleles for each loci
    }
    \item{loci}{
        number giving the place in List_gene of the loci to concider
    }
    \item{N}{
        Array of 2 dimension giving the size of the population for each genotype for each deme, first dimension is the genotypes and second the deme.  
    }
}
\details{
    %%  ~~ If necessary, more details than the description above ~~
}
\value{
    A square matrix with the size of the popultion for each genotype of the concidered loci in each deme.
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
The index with two loci each with two allele
List_gene<-c(2,2)
genotype.index.multilocus(List_gene)
    
             A1B1 A2B1 A1B2 A2B2
    A1B1    1    2    3    4
    A2B1    2    5    6    7
    A1B2    3    6    8    9
    A2B2    4    7    9   10
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
