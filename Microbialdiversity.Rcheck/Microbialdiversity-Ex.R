pkgname <- "Microbialdiversity"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Microbialdiversity')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("calculate_inverse_simpson")
### * calculate_inverse_simpson

flush(stderr()); flush(stdout())

### Name: calculate_inverse_simpson
### Title: Calculate Inverse Simpson Diversity Index
### Aliases: calculate_inverse_simpson

### ** Examples

# Example with equal abundances
calculate_inverse_simpson(c(1, 1, 1, 1))  # Should return 4

# Example with unequal abundances
calculate_inverse_simpson(c(10, 20, 30, 40))  # Should return ~3.57



cleanEx()
nameEx("calculate_observed_otus")
### * calculate_observed_otus

flush(stderr()); flush(stdout())

### Name: calculate_observed_otus
### Title: Calculate Observed OTUs
### Aliases: calculate_observed_otus

### ** Examples

# Example with some zero abundances
calculate_observed_otus(c(1, 0, 3, 0, 5))  # Should return 3



cleanEx()
nameEx("calculate_shannon")
### * calculate_shannon

flush(stderr()); flush(stdout())

### Name: calculate_shannon
### Title: Calculate Shannon Diversity Index
### Aliases: calculate_shannon

### ** Examples

# Example with equal abundances
calculate_shannon(c(1, 1, 1, 1))  # Should return ~1.386

# Example with unequal abundances
calculate_shannon(c(10, 20, 30, 40))  # Should return ~1.279



cleanEx()
nameEx("calculate_simpson")
### * calculate_simpson

flush(stderr()); flush(stdout())

### Name: calculate_simpson
### Title: Calculate Simpson Diversity Index This function calculates the
###   Simpson diversity index (1-D) for a vector of species counts or
###   proportions. The Simpson index ranges from 0 (low diversity) to 1
###   (high diversity).
### Aliases: calculate_simpson

### ** Examples

# Example with equal abundances
calculate_simpson(c(1, 1, 1, 1))  # Should return 0.75

# Example with unequal abundances
calculate_simpson(c(10, 20, 30, 40))  # Should return ~0.72



cleanEx()
nameEx("compare_inverse_simpson_multiple")
### * compare_inverse_simpson_multiple

flush(stderr()); flush(stdout())

### Name: compare_inverse_simpson_multiple
### Title: Compare Inverse Simpson Diversity Between Multiple Groups
### Aliases: compare_inverse_simpson_multiple

### ** Examples

communities <- matrix(c(10,20,30,40, 40,30,20,10, 25,25,25,25), 
                     nrow=4, ncol=3)
compare_inverse_simpson_multiple(communities)



cleanEx()
nameEx("compare_inverse_simpson_two")
### * compare_inverse_simpson_two

flush(stderr()); flush(stdout())

### Name: compare_inverse_simpson_two
### Title: Compare Inverse Simpson Diversity Between Two Groups
### Aliases: compare_inverse_simpson_two

### ** Examples

sample1 <- c(10, 20, 30, 40)
sample2 <- c(40, 30, 20, 10)
compare_inverse_simpson_two(sample1, sample2)



cleanEx()
nameEx("compare_observed_otus_multiple")
### * compare_observed_otus_multiple

flush(stderr()); flush(stdout())

### Name: compare_observed_otus_multiple
### Title: Compare Observed OTUs Between Multiple Groups
### Aliases: compare_observed_otus_multiple

### ** Examples

communities <- matrix(c(10,20,0,40, 40,30,20,0, 25,0,25,25), 
                     nrow=4, ncol=3)
compare_observed_otus_multiple(communities)



cleanEx()
nameEx("compare_observed_otus_two")
### * compare_observed_otus_two

flush(stderr()); flush(stdout())

### Name: compare_observed_otus_two
### Title: Compare Observed OTUs Between Two Groups
### Aliases: compare_observed_otus_two

### ** Examples

sample1 <- c(10, 20, 30, 0, 40)
sample2 <- c(40, 0, 30, 20, 10)
compare_observed_otus_two(sample1, sample2)



cleanEx()
nameEx("compare_shannon_multiple")
### * compare_shannon_multiple

flush(stderr()); flush(stdout())

### Name: compare_shannon_multiple
### Title: Compare Shannon Diversity Between Multiple Groups
### Aliases: compare_shannon_multiple

### ** Examples

communities <- matrix(c(10,20,30,40, 40,30,20,10, 25,25,25,25), 
                     nrow=4, ncol=3)
compare_shannon_multiple(communities)



cleanEx()
nameEx("compare_shannon_two")
### * compare_shannon_two

flush(stderr()); flush(stdout())

### Name: compare_shannon_two
### Title: Compare Shannon Diversity Between Two Groups
### Aliases: compare_shannon_two

### ** Examples

sample1 <- c(10, 20, 30, 40)
sample2 <- c(40, 30, 20, 10)
compare_shannon_two(sample1, sample2)



cleanEx()
nameEx("compare_simpson_multiple")
### * compare_simpson_multiple

flush(stderr()); flush(stdout())

### Name: compare_simpson_multiple
### Title: Compare Simpson Diversity Between Multiple Groups
### Aliases: compare_simpson_multiple

### ** Examples

communities <- matrix(c(10,20,30,40, 40,30,20,10, 25,25,25,25), 
                     nrow=4, ncol=3)
compare_simpson_multiple(communities)



cleanEx()
nameEx("compare_simpson_two")
### * compare_simpson_two

flush(stderr()); flush(stdout())

### Name: compare_simpson_two
### Title: Compare Simpson Diversity Between Two Groups
### Aliases: compare_simpson_two

### ** Examples

sample1 <- c(10, 20, 30, 40)
sample2 <- c(40, 30, 20, 10)
compare_simpson_two(sample1, sample2)



cleanEx()
nameEx("microbiome_data")
### * microbiome_data

flush(stderr()); flush(stdout())

### Name: microbiome_data
### Title: Example Microbiome Dataset
### Aliases: microbiome_data
### Keywords: datasets

### ** Examples

# Load the dataset
data(microbiome_data)

# Calculate Shannon diversity for each sample
shannon_values <- apply(microbiome_data$abundances, 2, calculate_shannon)

# Compare diversity between conditions
control_samples <- microbiome_data$abundances[, 
                  microbiome_data$metadata$condition == "Control"]
treatment_samples <- microbiome_data$abundances[, 
                    microbiome_data$metadata$condition == "Treatment"]
diversity_comparison <- compare_shannon_two(control_samples, treatment_samples)

# View first few rows of abundance data
head(microbiome_data$abundances)

# View sample metadata
print(microbiome_data$metadata)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
