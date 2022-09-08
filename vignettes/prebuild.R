# Precompile vignettes from the vignettes/prebuild/ folder by running
#   source("vignettes/prebuild.R")
# Then git commit the generated files.
#
# Vignettes in prebuild should set knitr option
#   fig.path = "figure/example/"
# where "example" is the root name of the vignette source file.
# This script will automatically move generated figures to the
# vignettes/figure/ folder.

local({
  setwd("vignettes")
  on.exit({setwd("..")})
  if (!dir.exists("figure")) {
    dir.create("figure", recursive = TRUE)
  }
  files <- list.files("prebuild", pattern = ".*\\.Rmd$")
  for (file in files) {
    file_root <- sub("\\.Rmd$", "", file)
    knitr::knit(file.path("prebuild", file), file)
  }
})
