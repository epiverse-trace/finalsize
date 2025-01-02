# finalsize (development version)

Maintainer is changing to @rozeggo (#212).

1. Updated all GitHub Actions workflows in line with {epiverse-trace/packagetemplate} (#212).

2. Updated DESCRIPTION and license files with new maintainer and new copyright year (#212).

3. Added `R/dev-utils.R` for extra release issue bullet points, and added `tools/check.env` for global environment checks (#212).

4. Corrected internal article links in vignettes (#212).

5. Added continuous benchmarking workflows using {touchstone} following the pattern of {epiforecasts/EpiNow2} (#212).

6. Updated `final_size()` to return the demography-susceptibility group sizes and the absolute value of individuals infected as columns.

# finalsize 0.2.1

This patch adds:

1. Default values for the contact matrix, susceptibility matrix, and susceptibility distribution matrix to allow quick estimates of final size for a given $R_0$ in a population with homogeneous mixing and full susceptibility (#180);

2. Combines the documentation for the functions `r_eff()`, `r0_to_lambda()`, and `lambda_to_r0()` under the name 'r0_conversions';

3. Adds vignettes on the theoretical background and on projecting re-emergence of a disease due to demographic turnover (#193, #194 by @adamkucharski);

4. Streamlines the examples in the Readme;

5. Replaces the PNG logo in the Readme with an SVG version;

6. Updates the GitHub Actions workflow YAMLs, the `.lintr` config file, `_pkgdown.yml`, `.Rbuildignore`, `CITATION.cff` and `CRAN-SUBMISSION`; workflows and configuration files bring the package in line with {packagetemplate} (#181, #186, #188).

7. Standardises the `DESCRIPTION` file.

8. Adds package-level documentation.

9. Adds spellchecking and the `WORDLIST` file.

10. Adds tools to check that the package does not modify the global state, code to flag any partial matching of function arguments,.

11. Updates references cited in the vignettes.

# finalsize 0.2.0

This is the second release of _finalsize_, and includes:

1. The package's C++ code for the solver algorithms has been moved from source files under `src/` into headers under `inst/include/`
2. The package includes a package header under `inst/include` called `finalsize.h`, which allows other Rcpp packages to link to _finalsize_ and reuse the solver algorithms provided here.
3. Three helper functions have been added to:
    - Calculate the effective basic reproductive number $R_{eff}$ in a population with heterogeneous social contacts and susceptibility to infection
    - Convert from the basic reproductive number provided by the user $R_0$ to the transmission rate (denoted by $\lambda$), given the population's social contacts and susceptibility structure
    - Convert from a user-provided transmission rate to the basic reproductive number given the population's social contacts and susceptibility structure
4. Two new vignettes have been added which are intended to serve as contextual background information for users:
    - A guide to constructing susceptibility matrices
    - A comparison with a simple susceptible-infectious-recovered epidemic model
5. Package infrastructure has been modified:
    - The package now specifies C++17 as standard
    - The `Readme.Rmd` uses auto-rendering of the package and repository name to `Readme.md` via an updated `render-readme` Github Actions workflow
    - The Cpplint workflow now includes linting and checking using Cppcheck for header files
6. Updated `NEWS.md` file to track changes to the package.

# finalsize 0.1

Initial release of _finalsize_, an R package to calculate the final size of an epidemic in a population with demographic variation in social contacts and in susceptibility to infection.

This release includes:

1. A choice of equation solver functions.
2. 100% code coverage,
3. A basic usage vignette, and two advanced vignettes,
4. Example data from the POLYMOD dataset obtained using the _socialmixr_ R package,
5. Workflows to render the vignettes and README as a website.
