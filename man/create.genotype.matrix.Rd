\name{create.genotype.matrix}
\alias{create.genotype.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create genotype matrix
}
\description{
  This function returns a matrix giving all the possible multilocus genotype. In column is given the different genotypes and in row the alleles. 
}
\usage{
  create.genotype.matrice(List_gene)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List_gene}{
    A list giving the number of loci and allele for each loci. The length of the list correspond to the number of loci, and each number indecates the number of alleles for this loci. 
  }
  \item{index}{
  is a matrix calculated with the function genotype.index.multilocus, indicates the row of each genotype in the simulation. 
}
}

\examples{
  List_gene<-c(2,2) # example of two loci each with two alleles 
  index<-genotype.index.multilocus(List_gene)
  create.genotype.matrix(List_gene,index)
  
            1 2 1 2
A1B1/A1B1 1 1 1 1
A1B1/A2B1 1 2 1 1
A1B1/A1B2 1 1 1 2
A1B1/A2B2 1 2 1 2
A2B1/A2B1 2 2 1 1
A2B1/A1B2 2 1 1 2
A2B1/A2B2 2 2 1 2
A1B2/A1B2 1 1 2 2
A1B2/A2B2 1 2 2 2
A2B2/A2B2 2 2 2 2
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
