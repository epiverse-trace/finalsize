
# *finalsize*: Calculate the final size of an epidemic <img src="man/figures/logo.png" align="right" width="130"/>

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/epiverse-trace/finalsize/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/epiverse-trace/finalsize/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/epiverse-trace/finalsize/branch/main/graph/badge.svg)](https://app.codecov.io/gh/epiverse-trace/finalsize?branch=main)
<!-- badges: end -->

*finalsize* provides quick back-of-the-envelope calculations for the
final size of an epidemic outbreak in a population with demographic
variation in contact patterns, and variation within and between age
groups in their susceptibility to disease.

*finalsize* implements methods outlined in [Andreasen
(2011)](https://link.springer.com/article/10.1007/s11538-010-9623-3),
[Miller (2012)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3506030/),
[Kucharski et al.
(2014)](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1004206),
and [Bidari et al.
(2016)](https://www.sciencedirect.com/science/article/pii/S0025556416302802).

*finalsize* can help provide rough estimates of the effectiveness of
pharmaceutical interventions in the form of immunisation programmes, or
the effect of naturally acquired immunity through previous infection
(see the vignette).

*finalsize* relies on [Eigen](https://gitlab.com/libeigen/eigen) via
[RcppEigen](https://github.com/RcppCore/RcppEigen) for fast matrix
algebra, and is developed at the [Centre for the Mathematical Modelling
of Infectious
Diseases](https://www.lshtm.ac.uk/research/centres/centre-mathematical-modelling-infectious-diseases)
at the London School of Hygiene and Tropical Medicine as part of the
[Epiverse Initiative](https://data.org/initiatives/epiverse/).

## Installation

The package can be installed using

``` r
install.packages("finalsize")
```

The current development version of *finalsize* can be installed from
[Github](https://github.com/epiverse-trace/finalsize) using the
`remotes` package.

``` r
# install.packages("devtools")
devtools::install_github("epiverse-trace/finalsize")
```

## Quick start

*finalsize* provides the single function `final_size()`, to calculate
the final size of an epidemic.

Here, an example using social contact data from the *socialmixr* package
investigates the final size of an epidemic when the disease has an
R<sub>0</sub> of 2.0, and given three age groups of interest — 0-19,
20-39 and 40+. The under-20 age group is assumed to be fully susceptible
to the disease, whereas individuals aged over 20 are only half as
susceptible.

``` r
# Load socialmixr package
library(socialmixr)
#> 
#> Attaching package: 'socialmixr'
#> The following object is masked from 'package:utils':
#> 
#>     cite

# load finalsize
library(finalsize)

# Load data from POLYMOD
data(polymod)
contact_data <- contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)
#> Using POLYMOD social contact data. To cite this in a publication, use the 'cite' function
#> Removing participants that have contacts without age information. To change this behaviour, set the 'missing.contact.age' option

# Define contact matrix (entry {ij} is contacts in group i reported by group j)
contact_matrix <- t(contact_data$matrix)

# Define population in each age group
demography_vector <- contact_data$demography$population

# Scale the contact matrix to ensure its largest eigenvalue is 1.0
contact_matrix <- contact_matrix / max(eigen(contact_matrix)$values)

# Divide the contact matrix by the demography
# Each i-th row divided by the corresponding i-th element of demography
contact_matrix <- contact_matrix / demography_vector

# Define susceptibility of each group
susceptibility <- matrix(
  data = c(1.0, 0.5, 0.5),
  nrow = length(demography_vector),
  ncol = 1
)

# Assume uniform susceptibility within age groups
p_susceptibility <- matrix(
  data = 1.0,
  nrow = length(demography_vector),
  ncol = 1
)

# R0 of the disease
r0 <- 1.3 # seasonal influenza

# Calculate the proportion of individuals infected
final_size(
  r0,
  contact_matrix,
  demography_vector,
  p_susceptibility,
  susceptibility
)
#>   demo_grp   susc_grp susceptibility  p_infected
#> 1   [0,20) susc_grp_1            1.0 0.036325712
#> 2  [20,40) susc_grp_1            0.5 0.009726187
#> 3      40+ susc_grp_1            0.5 0.006240606
```

## Package vignettes

More details on how to use *finalsize* can be found in the [online
documentation as package
vignettes](https://epiverse-trace.github.io/finalsize/), under
“Articles”.

## Help

To report a bug please open an
[issue](https://github.com/epiverse-trace/finalsize/issues/new/choose).

## Contribute

Contributions to *finalsize* are welcomed. Please follow the [package
contributing
guide](https://github.com/epiverse-trace/finalsize/blob/main/.github/CONTRIBUTING.md).

## Code of conduct

Please note that the *finalsize* project is released with a [Contributor
Code of
Conduct](https://github.com/epiverse-trace/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## Citing this package

``` r
citation("finalsize")
#> 
#> To cite package 'finalsize' in publications use:
#> 
#>   Gupte P, Kucharski A (2022). _finalsize: Calculates Final Size for
#>   SIR Epidemic in Heterogeneous Population_. R package version 0.1.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {finalsize: Calculates Final Size for SIR Epidemic in Heterogeneous
#> Population},
#>     author = {Pratik Gupte and Adam Kucharski},
#>     year = {2022},
#>     note = {R package version 0.1},
#>   }
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-andreasen2011" class="csl-entry">

Andreasen, Viggo. 2011. “The Final Size of an Epidemic and Its Relation
to the Basic Reproduction Number.” *Bulletin of Mathematical Biology* 73
(10): 2305–21. <https://doi.org/10.1007/s11538-010-9623-3>.

</div>

<div id="ref-bidari2016" class="csl-entry">

Bidari, Subekshya, Xinying Chen, Daniel Peters, Dylanger Pittman, and
Péter L. Simon. 2016. “Solvability of Implicit Final Size Equations for
SIR Epidemic Models.” *Mathematical Biosciences* 282 (December): 181–90.
<https://doi.org/10.1016/j.mbs.2016.10.012>.

</div>

<div id="ref-kucharski2014" class="csl-entry">

Kucharski, Adam J., Kin O. Kwok, Vivian W. I. Wei, Benjamin J. Cowling,
Jonathan M. Read, Justin Lessler, Derek A. Cummings, and Steven Riley.
2014. “The Contribution of Social Behaviour to the Transmission of
Influenza A in a Human Population.” *PLoS Pathogens* 10 (6): e1004206.
<https://doi.org/10.1371/journal.ppat.1004206>.

</div>

<div id="ref-miller2012" class="csl-entry">

Miller, Joel C. 2012. “A Note on the Derivation of Epidemic Final
Sizes.” *Bulletin of Mathematical Biology* 74 (9): 2125–41.
<https://doi.org/10.1007/s11538-012-9749-6>.

</div>

</div>
