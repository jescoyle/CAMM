## This script is used for building and maintaining the CAMM package
#devtools::install_github("hadley/staticdocs")

# Load packages needed to build CTSim
library(devtools)
library(roxygen2)
library(staticdocs)

setwd('C:/Users/jrcoyle/Documents/Research/CAMM/GitHub/')

current_code = as.package('CAMM')

# Load functions
load_all(current_code)

# Update documentation
document(current_code)

# Add Imports and Suggests to DESCRIPTION
setwd('./CTSim')
use_package('abind')
use_package('fdrtool')
use_package('gstat')
use_package('poweRlaw')
use_package('raster')
use_package('reshape2')
use_package('sads')
use_package('sp')
use_package('doParallel','Suggests')
use_package('foreach','Suggests')

# Check the package
setwd('../')
check('CTSim')

# Build the package
build('CTSim')
#build_win('CTSim')


# Check install
install.packages('CTSim_0.1.3.zip', repos=NULL)
library(CTSim)

# Make static html documentation
build_site('CTSim', 'HTML')