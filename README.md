# finalsize

`R` package to calculate final size for SIR epidemic in a heterogenous population.

## Installation

The easiest way to install the development version of `finalsize` is to use the `devtools` package:

```r
# install.packages("devtools")
library(devtools)
install_github("epiverse-trace/finalsize")
library(finalsize)

```

## Quick start

Use `final_size` to calculate the final size of an epidemic using social contact data from the `socialmixr` package. In the below example, we consider three age groups: 0-19, 20-39 and 40+, with R0 = 2. We also assume that under 20 age group is fully susceptible, whereas only 50% of individuals age 20+ are.

```
# Load socialmixr package

# install.packages("socialmixr")
library(socialmixr)

# Load data from POLYMOD
data(polymod)
contact_data <- contact_matrix(polymod, countries = "United Kingdom", age.limits = c(0,20,40))

c_matrix <-t(contact_data$matrix) # Define contact matrix (entry ij is contacts in group i reported by group j)
d_vector <- contact_data$participants$proportion # Define proportion in each age group
p_suscep <- c(1,0.5,0.5) # Define proportion of age group that is susceptible to infection

# Run final size model
final_size(r0=2,contact_matrix = c_matrix, demography_vector = d_vector, prop_suscep = p_suscep)
```


## Citation

Kucharski AJ, Kwok KO, Wei VW, Cowling BJ, Read JM, Lessler J, Cummings DA, Riley S. [The contribution of social behaviour to the transmission of influenza A in a human population](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1004206). PLOS Pathogens 2014;10(6):e1004206 PMID: 24968312