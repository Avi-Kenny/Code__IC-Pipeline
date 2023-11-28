# Change library path
.libPaths(c("/home/akenny/R_lib", .libPaths()))

# Set packages
cfg$pkgs <- c("splines", "survival", "SuperLearner")
cfg$pkgs_nocluster <- c("ggplot2")

# Set cluster config
if (Sys.getenv("HOME")=="/home/akenny") {
  # Bionic
  cluster_config <- list(
    js = "slurm",
    dir = paste0("/home/akenny/", Sys.getenv("project"))
  )
} else if (Sys.getenv("HOME")=="/home/users/avikenny") {
  # Bayes
  cluster_config <- list(
    js = "ge",
    dir = paste0("/home/users/avikenny/Desktop/", Sys.getenv("project"))
  )
} else {
  cluster_config <- list(js="", dir="")
}

# Load packages
for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
  suppressMessages({ do.call("library", list(pkg)) })
}
