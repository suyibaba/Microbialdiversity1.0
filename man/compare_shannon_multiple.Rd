% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shannon_diversity.R
\name{compare_shannon_multiple}
\alias{compare_shannon_multiple}
\title{Compare Shannon Diversity Between Multiple Groups}
\usage{
compare_shannon_multiple(abundance_matrix, shannon_base = exp(1))
}
\arguments{
\item{abundance_matrix}{Matrix where rows are species and columns are samples}

\item{shannon_base}{Base for logarithm calculation (default: exp(1))}
}
\value{
Matrix of p-values for pairwise comparisons
}
\description{
Compare Shannon Diversity Between Multiple Groups
}
\examples{
communities <- matrix(c(10,20,30,40, 40,30,20,10, 25,25,25,25), 
                     nrow=4, ncol=3)
compare_shannon_multiple(communities)
}
\author{
Akinsuyi Oluwamayowa Samuel
}
