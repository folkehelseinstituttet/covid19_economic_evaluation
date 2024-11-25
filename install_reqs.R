
options(
  repos = structure(c(
 #   SPLVERSE  = "https://docs.sykdomspulsen.no/drat/",
    CRAN      = "https://cran.rstudio.com"
  ))
)

drat:::add("ncov-ic")
install.packages(c("abind", "tidyverse", "data.table", "odin", "odin.dust", "dust", "mcstate", "GA", "tryCatchLog"))

install.packages("spldata")


devtools::install_github(repo="https://github.com/Gulfa/metapop")
devtools::install_github(repo="https://github.com/Gulfa/metapopnorge")
devtools::install_github(repo="https://github.com/Gulfa/metapopdeterministic")


