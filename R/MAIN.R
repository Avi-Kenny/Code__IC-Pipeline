# Main config
cfg <- list(
  run_analysis = T
)

# Secondary config
source("R/config.R", local=T)

# Analysis pipeline
if (cfg$run_analysis) { source("R/analysis.R", local=T) }
