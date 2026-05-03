###############################################################################
# install_packages.R
#
# Pinned dependency installer for CodeThesis.R. Run this once before sourcing
# the analysis script:
#
#   source("install_packages.R")
#   source("CodeThesis.R")
#
# The script uses `remotes::install_version()` so that future CRAN updates do
# not silently change the rendered output of the thesis pipeline.
###############################################################################

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

# Packages actually loaded by CodeThesis.R, with the versions used to produce
# the rendered results. Update only after re-running and re-rendering the
# whole pipeline.
pinned <- list(
  lpSolve = "5.6.21",   # LP construction and solver (lp() builds + solves)
  beepr   = "2.0",      # Audible solver-status notifications
  ggplot2 = "3.5.1"     # qplot() for the budget-sensitivity charts in §5.4
)

repo <- "https://cloud.r-project.org"

for (pkg in names(pinned)) {
  ver <- pinned[[pkg]]
  installed <- tryCatch(
    as.character(utils::packageVersion(pkg)),
    error = function(e) NA_character_
  )
  if (is.na(installed) || installed != ver) {
    message(sprintf("Installing %s %s ...", pkg, ver))
    remotes::install_version(
      package  = pkg,
      version  = ver,
      repos    = repo,
      upgrade  = "never",
      quiet    = FALSE
    )
  } else {
    message(sprintf("%s %s already installed.", pkg, ver))
  }
}

invisible(NULL)
