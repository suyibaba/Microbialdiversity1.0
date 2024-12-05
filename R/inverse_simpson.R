##' Calculate Inverse Simpson Diversity Index
##'
##' This function calculates the Inverse Simpson diversity index for a vector of species counts or proportions.
##' The Inverse Simpson index ranges from 1 (low diversity) to the number of species (high diversity).
##'
##' @param abundances A numeric vector of species counts or proportions.
##' @return The Inverse Simpson diversity index.
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' # Example with equal abundances
##' calculate_inverse_simpson(c(1, 1, 1, 1))  # Should return 4
##'
##' # Example with unequal abundances
##' calculate_inverse_simpson(c(10, 20, 30, 40))  # Should return ~3.57
##' @export
calculate_inverse_simpson <- function(abundances) {
  if (any(abundances < 0)) {
    stop("Abundances must be non-negative.")
  }
  proportions <- abundances / sum(abundances)
  inverse_simpson_index <- 1 / sum(proportions^2)
  return(inverse_simpson_index)
}

##' Compare Inverse Simpson Diversity Between Two Groups
##'
##' @param x Numeric vector of abundances for first community
##' @param y Numeric vector of abundances for second community
##' @return List containing test statistics and p-value
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' sample1 <- c(10, 20, 30, 40)
##' sample2 <- c(40, 30, 20, 10)
##' compare_inverse_simpson_two(sample1, sample2)
##' @export
compare_inverse_simpson_two <- function(x, y) {
  # Calculate Inverse Simpson indices
  simpson_x <- calculate_inverse_simpson(x)
  simpson_y <- calculate_inverse_simpson(y)
  
  # Bootstrap comparison
  n_boot <- 1000
  boot_diff <- numeric(n_boot)
  
  for(i in 1:n_boot) {
    # Resample with replacement
    boot_x <- sample(x, length(x), replace = TRUE)
    boot_y <- sample(y, length(y), replace = TRUE)
    
    # Calculate difference in Inverse Simpson indices
    boot_diff[i] <- calculate_inverse_simpson(boot_x) - calculate_inverse_simpson(boot_y)
  }
  
  # Calculate p-value
  observed_diff <- simpson_x - simpson_y
  p_value <- mean(abs(boot_diff) >= abs(observed_diff))
  
  return(list(
    inverse_simpson_x = simpson_x,
    inverse_simpson_y = simpson_y,
    difference = observed_diff,
    p_value = p_value
  ))
}

##' Compare Inverse Simpson Diversity Between Multiple Groups
##'
##' @param abundance_matrix Matrix where rows are species and columns are samples
##' @return Matrix of p-values for pairwise comparisons
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' communities <- matrix(c(10,20,30,40, 40,30,20,10, 25,25,25,25), 
##'                      nrow=4, ncol=3)
##' compare_inverse_simpson_multiple(communities)
##' @export
compare_inverse_simpson_multiple <- function(abundance_matrix) {
  n_communities <- ncol(abundance_matrix)
  p_values <- matrix(NA, n_communities, n_communities)
  
  for(i in 1:(n_communities-1)) {
    for(j in (i+1):n_communities) {
      result <- compare_inverse_simpson_two(abundance_matrix[,i], 
                                            abundance_matrix[,j])
      p_values[i,j] <- result$p_value
      p_values[j,i] <- result$p_value
    }
  }
  
  diag(p_values) <- 1
  return(p_values)
}

# Example usage:
# Simple diversity calculation
community1 <- c(10, 10, 10, 10)
print(calculate_inverse_simpson(community1))  # Should return 4

# Comparing two communities
community_A <- c(35, 19, 11)  # Uneven distribution
community_B <- c(20, 20, 25)  # More even distribution
result_two <- compare_inverse_simpson_two(community_A, community_B)
print(result_two)

# Comparing multiple communities
communities_matrix <- matrix(
  c(35, 19, 11,   # Community A
    20, 20, 25,   # Community B
    30, 15, 20),  # Community C
  nrow = 3, 
  ncol = 3, 
  byrow = TRUE
)
result_multiple <- compare_inverse_simpson_multiple(communities_matrix)
print(result_multiple)