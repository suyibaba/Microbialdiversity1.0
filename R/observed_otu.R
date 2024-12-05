##' Calculate Observed OTUs
##'
##' This function calculates the number of observed OTUs (Operational Taxonomic Units) 
##' or ASVs (Amplicon Sequence Variants) in a sample by counting non-zero abundances.
##'
##' @param abundances A numeric vector of species counts or proportions.
##' @return The number of observed OTUs/ASVs (species richness).
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' # Example with some zero abundances
##' calculate_observed_otus(c(1, 0, 3, 0, 5))  # Should return 3
##' @export
calculate_observed_otus <- function(abundances) {
  if (any(abundances < 0)) {
    stop("Abundances must be non-negative.")
  }
  observed <- sum(abundances > 0)
  return(observed)
}

##' Compare Observed OTUs Between Two Groups
##'
##' @param x Numeric vector of abundances for first community
##' @param y Numeric vector of abundances for second community
##' @return List containing test statistics and p-value
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' sample1 <- c(10, 20, 30, 0, 40)
##' sample2 <- c(40, 0, 30, 20, 10)
##' compare_observed_otus_two(sample1, sample2)
##' @export
compare_observed_otus_two <- function(x, y) {
  # Calculate observed OTUs
  obs_x <- calculate_observed_otus(x)
  obs_y <- calculate_observed_otus(y)
  
  # Bootstrap comparison
  n_boot <- 1000
  boot_diff <- numeric(n_boot)
  
  for(i in 1:n_boot) {
    # Resample with replacement
    boot_x <- sample(x, length(x), replace = TRUE)
    boot_y <- sample(y, length(y), replace = TRUE)
    
    # Calculate difference in observed OTUs
    boot_diff[i] <- calculate_observed_otus(boot_x) - calculate_observed_otus(boot_y)
  }
  
  # Calculate p-value
  observed_diff <- obs_x - obs_y
  p_value <- mean(abs(boot_diff) >= abs(observed_diff))
  
  return(list(
    observed_otus_x = obs_x,
    observed_otus_y = obs_y,
    difference = observed_diff,
    p_value = p_value
  ))
}

##' Compare Observed OTUs Between Multiple Groups
##'
##' @param abundance_matrix Matrix where rows are species and columns are samples
##' @return Matrix of p-values for pairwise comparisons
##' @author
##' Akinsuyi Oluwamayowa Samuel
##' @examples
##' communities <- matrix(c(10,20,0,40, 40,30,20,0, 25,0,25,25), 
##'                      nrow=4, ncol=3)
##' compare_observed_otus_multiple(communities)
##' @export
compare_observed_otus_multiple <- function(abundance_matrix) {
  n_communities <- ncol(abundance_matrix)
  p_values <- matrix(NA, n_communities, n_communities)
  
  for(i in 1:(n_communities-1)) {
    for(j in (i+1):n_communities) {
      result <- compare_observed_otus_two(abundance_matrix[,i], 
                                          abundance_matrix[,j])
      p_values[i,j] <- result$p_value
      p_values[j,i] <- result$p_value
    }
  }
  
  diag(p_values) <- 1
  return(p_values)
}

# Example usage:
# Simple OTU calculation
community1 <- c(10, 0, 10, 10)
print(calculate_observed_otus(community1))  # Should return 3

# Comparing two communities
community_A <- c(35, 0, 11, 0, 15)  # Has 3 OTUs
community_B <- c(20, 20, 0, 25, 0)  # Has 3 OTUs
result_two <- compare_observed_otus_two(community_A, community_B)
print(result_two)

# Comparing multiple communities
communities_matrix <- matrix(
  c(35, 0, 11, 0, 15,    # Community A
    20, 20, 0, 25, 0,    # Community B
    30, 15, 20, 0, 0),   # Community C
  nrow = 5, 
  ncol = 3, 
  byrow = FALSE
)
result_multiple <- compare_observed_otus_multiple(communities_matrix)
print(result_multiple)