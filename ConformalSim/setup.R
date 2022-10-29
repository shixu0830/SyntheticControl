if (!require("gmm", lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"),
             quietly = TRUE))
  install.packages("gmm", repos = "https://cloud.r-project.org",
                   lib = "/home/qijunli/R_libraries",
                   dependencies = TRUE)

if (!require("sandwich", lib.loc = c(.libPaths(), "/home/qijunli/R_libraries"),
             quietly = TRUE))
  install.packages("sandwich", repos = "https://cloud.r-project.org",
                   lib = "/home/qijunli/R_libraries",
                   dependencies = TRUE)

print(parallel::detectCores())