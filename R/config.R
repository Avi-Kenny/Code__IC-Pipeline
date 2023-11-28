# Change library path
.libPaths(c("/home/akenny/R_lib", .libPaths()))

# Set packages
cfg$pkgs <- c(
  "splines", "survival", "SuperLearner", "magrittr", "ggplot2", "dplyr"
)
cfg$pkgs_nocluster <- c()

# Set cluster config
if (Sys.getenv("HOME")=="/home/akenny") {
  # Bionic
  cluster_config <- list(
    js = "slurm",
    dir = paste0("/home/akenny/", Sys.getenv("proj"),
                 "/Code__", Sys.getenv("proj"))
  )
} else if (Sys.getenv("HOME")=="/home/users/avikenny") {
  # Bayes
  cluster_config <- list(
    js = "ge",
    dir = paste0("/home/users/avikenny/Desktop/", Sys.getenv("proj"),
                 "/Code__", Sys.getenv("proj"))
  )
} else {
  cluster_config <- list(js="", dir="")
}

# Load packages
for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
  suppressMessages({ do.call("library", list(pkg)) })
}
