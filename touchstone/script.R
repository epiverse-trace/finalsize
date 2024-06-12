# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# TODO OPTIONAL Add directories you want to be available in this file or during
# the benchmarks.
# touchstone::pin_assets("some/dir")

# installs branches to benchmark
touchstone::branch_install()

#### Simple final size calculation on uniform population ####
# run for 1000 values of R0
# using iterative solver
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup.R")
  },
  default_iterative = {
    lapply(r0_samples, final_size, solver = "iterative")
  },
  n = 100
)

# using newton solver
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup.R")
  },
  default_newton = {
    lapply(r0_samples, final_size, solver = "newton")
  },
  n = 100
)

#### Final size calculation on heterogeneous population ####
# run for 1000 values of R0
# using iterative solver
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup.R")
  },
  complex_iterative = {
    lapply(
      r0_samples, final_size,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = susceptibility,
      p_susceptibility = p_susceptibility,
      solver = "iterative"
    )
  },
  n = 100
)

# using newton solver
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup.R")
  },
  complex_newton = {
    lapply(
      r0_samples, final_size,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      susceptibility = susceptibility,
      p_susceptibility = p_susceptibility,
      solver = "newton"
    )
  },
  n = 100
)

# create artifacts used downstream in the GitHub Action
touchstone::benchmark_analyze()
