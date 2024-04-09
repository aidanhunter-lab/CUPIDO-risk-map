#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to run shiny app from a single 'source' command issued to R session in terminal
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load R shiny, installing if needed
pkg <- 'shiny'
# shiny.version <- '1.7.4'
# pkg_ <- paste0('https://cran.r-project.org/package=', pkg, '&version=', shiny.version)
got.pkg <- library(pkg, character.only = TRUE, logical.return = TRUE)
#if(got.pkg) right.version <- packageVersion(pkg) == shiny.version
if(!got.pkg){
  install.packages(pkgs = pkg)
  library(pkg, character.only = TRUE)
}

# Run the shiny app that is saved in the same directory as this file.
runApp()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Commented code below was supposed to automatically install the selected version
# of the shiny package, but I didn't figure out how to automatically install dependencies
# for the install of old versions...

# pkg <- 'shiny'
# shiny.version <- '1.7.4'
# pkg_ <- paste0('https://cran.r-project.org/package=', pkg, '&version=', shiny.version)
# got.pkg <- require(pkg, character.only = TRUE)
# if(got.pkg) right.version <- packageVersion(pkg) == shiny.version
# if(!got.pkg){
#   install.packages(pkgs = pkg_)# , repos = NULL)
#   library(pkg, character.only = TRUE)
# }else{
#   if(!right.version){
#     remove.packages(pkg)
#     install.packages(pkgs = pkg_)# , repos = NULL)
#     library(pkg, character.only = TRUE)
#   }
# }

# # Run the shiny app that is saved in the same directory as this file.
# runApp()
