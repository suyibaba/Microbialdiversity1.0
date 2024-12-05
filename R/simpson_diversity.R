#' Simpson Diversity Functions
#' @name simpson_diversity
#' @import ecolTest
NULL

#' Calculate Simpson Diversity Index
#' This function calculates the Simpson diversity index (1-D) for a vector of species counts or proportions.
#' The Simpson index ranges from 0 (low diversity) to 1 (high diversity).
#'
#' @param abundances A numeric vector of species counts or proportions.
#' @return The Simpson diversity index (1-D).
#' @author
#' Akinsuyi Oluwamayowa Samuel
#' @examples
#' # Example with equal abundances
#' calculate_simpson(c(1, 1, 1, 1))  # Should return 0.75
#'
#' # Example with unequal abundances
#' calculate_simpson(c(10, 20, 30, 40))  # Should return ~0.72
#' @export
calculate_simpson <- function(abundances) {
  if (any(abundances < 0)) {
    stop("Abundances must be non-negative.")
  }
  proportions <- abundances / sum(abundances)
  simpson_index <- 1 - sum(proportions^2)
  return(simpson_index)
}

#' Compare Simpson Diversity Between Two Groups
#'
#' @param x Numeric vector of abundances for first community
#' @param y Numeric vector of abundances for second community
#' @return List containing test statistics and p-value
#' @author
#' Akinsuyi Oluwamayowa Samuel
#' @examples
#' sample1 <- c(10, 20, 30, 40)
#' sample2 <- c(40, 30, 20, 10)
#' compare_simpson_two(sample1, sample2)
#' @export
compare_simpson_two <- function(x, y) {
  # Calculate Simpson indices
  simpson_x <- calculate_simpson(x)
  simpson_y <- calculate_simpson(y)
  
  # Bootstrap comparison
  n_boot <- 1000
  boot_diff <- numeric(n_boot)
  
  for(i in 1:n_boot) {
    # Resample with replacement
    boot_x <- sample(x, length(x), replace = TRUE)
    boot_y <- sample(y, length(y), replace = TRUE)
    
    # Calculate difference in Simpson indices
    boot_diff[i] <- calculate_simpson(boot_x) - calculate_simpson(boot_y)
  }
  
  # Calculate p-value
  observed_diff <- simpson_x - simpson_y
  p_value <- mean(abs(boot_diff) >= abs(observed_diff))
  
  return(list(
    simpson_x = simpson_x,
    simpson_y = simpson_y,
    difference = observed_diff,
    p_value = p_value
  ))
}

#' Compare Simpson Diversity Between Multiple Groups
#'
#' @param abundance_matrix Matrix where rows are species and columns are samples
#' @return Matrix of p-values for pairwise comparisons
#' @author
#' Akinsuyi Oluwamayowa Samuel
#' @examples
#' communities <- matrix(c(10,20,30,40, 40,30,20,10, 25,25,25,25), 
#'                      nrow=4, ncol=3)
#' compare_simpson_multiple(communities)
#' @export
compare_simpson_multiple <- function(abundance_matrix) {
  n_communities <- ncol(abundance_matrix)
  p_values <- matrix(NA, n_communities, n_communities)
  
  for(i in 1:(n_communities-1)) {
    for(j in (i+1):n_communities) {
      result <- compare_simpson_two(abundance_matrix[,i], 
                                    abundance_matrix[,j])
      p_values[i,j] <- result$p_value
      p_values[j,i] <- result$p_value
    }
  }
  
  diag(p_values) <- 1
  return(p_values)
}
