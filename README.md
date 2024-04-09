
# finalsize: Calculate the final size of an epidemic <img src="man/figures/logo.svg" align="right" width="130"/>

<!-- badges: start -->

<a href="https://app.digitalpublicgoods.net/a/10557"><img src="https://digitalpublicgoods.net/registry/dpgicon.svg" alt="Digital Public Good" height="25">
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/license/mit)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/epiverse-trace/finalsize/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/epiverse-trace/finalsize/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/epiverse-trace/finalsize/branch/main/graph/badge.svg)](https://app.codecov.io/gh/epiverse-trace/finalsize?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/finalsize)](https://CRAN.R-project.org/package=finalsize)
<!-- badges: end -->

*finalsize* is an R package to calculate the final size of a SIR
epidemic in populations with heterogeneity in social contacts and
infection susceptibility.

*finalsize* provides estimates for the total proportion of a population
infected over the course of an epidemic, and can account for a
demographic distribution (such as age groups) and demography-specific
contact patterns, as well as for heterogeneous susceptibility to
infection between groups (such as due to age-group specific immune
responses) and within groups (such as due to immunisation programs). An
advantage of this approach is that it requires fewer parameters to be
defined compared to a model that simulates the full transmission
dynamics over time, such as models in the [*epidemics*
package](https://epiverse-trace.github.io/epidemics/articles/epidemics.html).

*finalsize* implements methods outlined in Andreasen
([2011](#ref-andreasen2011)), Miller ([2012](#ref-miller2012)),
Kucharski et al. ([2014](#ref-kucharski2014)), and Bidari et al.
([2016](#ref-bidari2016)).

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
[Epiverse-TRACE](https://data.org/initiatives/epiverse/).

## Installation

The package can be installed from CRAN using

``` r
install.packages("finalsize")
```

### Development version

The current development version of *finalsize* can be installed from
[Github](https://github.com/epiverse-trace/finalsize) using the
`remotes` package. The development version documentation can be found
[here](https://epiverse-trace.github.io/finalsize/dev/).

``` r
if(!require("pak")) install.packages("pak")
pak::pak("epiverse-trace/finalsize")
```

## Quick start

The main function in *finalsize* is `final_size()`, which calculates the
final size of an epidemic given the $R_0$.

``` r
# load finalsize
library(finalsize)

final_size(1.5)
#>     demo_grp   susc_grp susceptibility p_infected
#> 1 demo_grp_1 susc_grp_1              1  0.5828132
```

Optionally, `final_size()` can estimate the epidemic size for
populations with differences among demographic groups in their social
contact patterns, in their susceptibility to infection.

We can use social contact data (here, from the *socialmixr* package) to
estimate the final size of an epidemic when the disease has an
R<sub>0</sub> of 1.5, and given three age groups of interest — 0-19,
20-39 and 40+. The under-20 age group is assumed to be fully susceptible
to the disease, whereas individuals aged over 20 are only half as
susceptible as those under 20.

``` r
# Load example POLYMOD social contacts data included with the package
data(polymod_uk)

# Define contact matrix (entry {ij} is contacts in group i reported by group j)
contact_matrix <- polymod_uk$contact_matrix

# Define population in each age group
demography_vector <- polymod_uk$demography_vector

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
r0 <- 1.5 # assumed for pandemic influenza

# Calculate the proportion of individuals infected in each age group
final_size(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)
#>   demo_grp   susc_grp susceptibility p_infected
#> 1   [0,20) susc_grp_1            1.0 0.32849966
#> 2  [20,40) susc_grp_1            0.5 0.10532481
#> 3      40+ susc_grp_1            0.5 0.06995193
```

Helper functions included in *finalsize* are provided to calculate the
effective $R_0$, called $R_{eff}$, from demographic and susceptibility
distribution data, while other helpers can convert between $R_0$ and the
transmission rate $\lambda$.

``` r
# calculate the effective R0 using `r_eff()`
r_eff(
  r0 = r0,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)
#> [1] 1.171758
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
#> To cite package 'finalsize' in publications use:
#> 
#>   Gupte P, Van Leeuwen E, Kucharski A (2024). _finalsize: Calculate the
#>   Final Size of an Epidemic_. R package version 0.2.1,
#>   https://epiverse-trace.github.io/finalsize/,
#>   <https://github.com/epiverse-trace/finalsize>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {finalsize: Calculate the Final Size of an Epidemic},
#>     author = {Pratik Gupte and Edwin {Van Leeuwen} and Adam Kucharski},
#>     year = {2024},
#>     note = {R package version 0.2.1, 
#> https://epiverse-trace.github.io/finalsize/},
#>     url = {https://github.com/epiverse-trace/finalsize},
#>   }
```

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

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
2014. “The contribution of social behaviour to the transmission of
influenza A in a human population.” *PLoS pathogens* 10 (6): e1004206.
<https://doi.org/10.1371/journal.ppat.1004206>.

</div>

<div id="ref-miller2012" class="csl-entry">

Miller, Joel C. 2012. “A Note on the Derivation of Epidemic Final
Sizes.” *Bulletin of Mathematical Biology* 74 (9): 2125–41.
<https://doi.org/10.1007/s11538-012-9749-6>.

</div>

</div>
