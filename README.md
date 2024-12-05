# Test Microbialdiversity1.0
This is a nice package use for comparing microbial diversity between two or more groups

## How to install the R package
devtools::install_github("YOUR_GITHUB_USERNAME/Microbialdiversity1.0", build_vignettes = TRUE)


## Example
data(microbiome_data)

# Calculate diversity indices for a single sample
sample1 <- microbiome_data$abundances[,1]

# Calculate different diversity measures
shannon <- calculate_shannon(sample1)
simpson <- calculate_simpson(sample1)
inv_simpson <- calculate_inverse_simpson(sample1)
otus <- calculate_observed_otus(sample1)

# Display results
results <- data.frame(
  Measure = c("Shannon", "Simpson", "Inverse Simpson", "Observed OTUs"),
  Value = c(shannon, simpson, inv_simpson, otus)
)
print(results)

# Compare diversity between two groups
control_samples <- microbiome_data$abundances[, microbiome_data$metadata$condition == "Control"]
treatment_samples <- microbiome_data$abundances[, microbiome_data$metadata$condition == "Treatment"]

# Perform statistical comparison
diversity_comparison <- compare_shannon_two(control_samples, treatment_samples)
print("Shannon diversity comparison between groups:")
print(diversity_comparison)
