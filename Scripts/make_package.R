## This script is used for building and maintaining the CAMM package
#devtools::install_github("hadley/staticdocs")

# Load packages needed to build CTSim
library(devtools)
library(roxygen2)
library(staticdocs)

setwd('C:/Users/jrcoyle/Documents/Research/CAMM/GitHub')

#create('CAMM')

current_code = as.package('CAMM')

# Load functions
load_all(current_code)

# Update documentation
document(current_code)

# Add Imports and Suggests to DESCRIPTION
# NOTE: NOT WORKING, added manually
setwd('./CAMM')
use_package('MASS')
use_package('sads')
use_package('vegan')
use_package('parallel','Suggests')


# Check the package
setwd('../')
check('CAMM')

# Build the package
build('CAMM')
build_win('CAMM')


# Check install
install.packages('CAMM_0.2.3.zip', repos=NULL)
library(CAMM)

# Make static html documentation
build_site('CAMM', 'HTML')