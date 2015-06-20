Description
=============

SpatialDelayDiff
* Is an R package for simulation, estimating, visualizing, and evaluating the performance of (1) delay-difference models, (2) mean-weight models, and (3) age-structured surplus production models, each with spatial and nonspatial versions.
* Contains a spatial simulator for generating spatial data, and for reducing it down to nonspatial indices of average weight and abundance


Instructions
=============
First, please install TMB (Template Model Builder) here: 
https://github.com/kaskr/adcomp

Next, please use R version >=3.1.1 and install the package:


    # Install package
    install.packages("devtools")
    library("devtools")
    install_github("James-Thorson/spatial_delay-difference")
    # Load package
    library(SpatialDelayDiff)

Please see examples folder for an example of how to run the model:
https://github.com/James-Thorson/spatial_delay-difference/blob/master/examples/Simulation_example.R

Further reading
=============

For more details regarding development and testing of this software please see:
* In press.  Thorson, J., Ianelli, J., Munch, S., Ono, K.‡, and Spencer, P. Spatial delay-difference models for estimating spatiotemporal variation in juvenile production and population abundance.  Canadian Journal of Fisheries and Aquatic Sciences
