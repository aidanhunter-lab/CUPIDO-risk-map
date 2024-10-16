#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Utility script to run R Shiny interactive map app via an R session run in a
# terminal using the command "source('run_app.R')". Necessary packages are
# automatically installed if possible, but note that some packages (sf) require
# manual installation of other software.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dr Aidan Hunter, Ecosystems group, British Antarctic Survey (2024).
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load R shiny, installing if needed
pkg <- 'shiny'
got.pkg <- library(pkg, character.only = TRUE, logical.return = TRUE)
if(!got.pkg){
  install.packages(pkgs = pkg)
  library(pkg, character.only = TRUE)
}

# Run the shiny app that is saved in the same directory as this file.
runApp()
