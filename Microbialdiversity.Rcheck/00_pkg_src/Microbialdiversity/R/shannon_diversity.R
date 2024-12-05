#' Shannon Diversity Functions
#' @name shannon_diversity
#' @import ecolTest
#' @importFrom utils install.packages
NULL

##' Calculate Shannon Diversity Index
##'
##' This function calculates the Shannon diversity index for a vector of species counts or proportions.
##'
##' @param abundances A numeric vector of species counts or proportions.
##' @return The Shannon diversity index.
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' # Example with equal abundances
##' calculate_shannon(c(1, 1, 1, 1))  # Should return ~1.386
##'
##' # Example with unequal abundances
##' calculate_shannon(c(10, 20, 30, 40))  # Should return ~1.279
##' @export
calculate_shannon <- function(abundances) {
  if (any(abundances < 0)) {
    stop("Abundances must be non-negative.")
  }
  # Normalize to proportions if not already
  proportions <- abundances / sum(abundances)
  
  # Calculate Shannon diversity
  shannon_index <- -sum(proportions * log(proportions))
  
  return(shannon_index)
}

##' Compare Shannon Diversity Between Two Groups
##'
##' @param x Numeric vector of abundances for first community
##' @param y Numeric vector of abundances for second community
##' @param shannon_base Base for logarithm calculation (default: exp(1))
##' @return List containing test statistics and p-value
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' sample1 <- c(10, 20, 30, 40)
##' sample2 <- c(40, 30, 20, 10)
##' compare_shannon_two(sample1, sample2)
##' @importFrom ecolTest Hutcheson_t_test
##' @export
compare_shannon_two <- function(x, y, shannon_base = exp(1)) {
  if (!requireNamespace("ecolTest", quietly = TRUE)) {
    install.packages("ecolTest")
  }
  library(ecolTest)
  result <- Hutcheson_t_test(x, y, shannon.base = shannon_base)
  return(result)
}

##' Compare Shannon Diversity Between Multiple Groups
##'
##' @param abundance_matrix Matrix where rows are species and columns are samples
##' @param shannon_base Base for logarithm calculation (default: exp(1))
##' @return Matrix of p-values for pairwise comparisons
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' communities <- matrix(c(10,20,30,40, 40,30,20,10, 25,25,25,25), 
##'                      nrow=4, ncol=3)
##' compare_shannon_multiple(communities)
##' @importFrom ecolTest multiple_Hutcheson_t_test
##' @export
compare_shannon_multiple <- function(abundance_matrix, shannon_base = exp(1)) {
  if (!requireNamespace("ecolTest", quietly = TRUE)) {
    install.packages("ecolTest")
  }
  library(ecolTest)
  result <- multiple_Hutcheson_t_test(abundance_matrix, 
                                      shannon.base = shannon_base)
  return(result)
}