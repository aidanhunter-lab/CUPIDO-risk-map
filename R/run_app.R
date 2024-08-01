#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to run shiny app from a single 'source' command issued to R session in terminal
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load R shiny, installing if needed
pkg <- 'shiny'
got.pkg <- library(pkg, character.only = TRUE, logical.return = TRUE)
if(!got.pkg){
  install.packages(pkgs = pkg)
  library(pkg, character.only = TRUE)
}

# Run the shiny app that is saved in the same directory as this file.
runApp()
