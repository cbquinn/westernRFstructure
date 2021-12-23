# westernRFstructure

This folder contains all input files and scripts necessary to create Figure 5 in "Contrasting genetic trajectories of endangered and expanding red fox populations in the western U.S."

Figure 5 is a visualization of genetic diversity summary statistics and effective population size (Ne) as a spatially explicit continuous surface

Visualization relies heavily on sGD package (with funcitons modified slighltly to handle haplotype data (mitochondrial Y chromosome gene diversity)

See **Shirk A, Cushman S (2011). sGD: software for estimating spatially explicit indices of genetic diversity. Molecular Ecology Resources 11(5): 922-934**

## Description of overall approach
Estimates of genetic diversity were calculated using two approaches based on results of genetic structure analyses: (a) populations with strong genetic structure that were determined to be discrete were calculated in the traditional manner in which one summary statistic is estimated based on all genotypes located within a minimum convex polygon, or (b) populations with relatively continuous genetic structure, in which an estimate is calculated for every individual based on genotypes that fall within a pre-defined radius of its location. For the latter, we used the sGD package (Shirk & Cushman 2011). 

To create plots, we assigned every sample a genetic diversity estimate based on either its membership to discrete population (a, above) or its genetic neighborhood (b, above). We  used inverse weighting to interpolate a spatially explicit continuous surface of genetic diversity in regions in 100 km buffer around the data.  

## Workflow

To create plot, run script `VisualizeSpatialPatternsOfNeandDiversity.R`

This script generates plots for autsomal diversity and calls on other scripts to generate plots for mitochondrial and y-chromosome haplotype. It also uses mutliple helper scripts with functions (e.g., for interpolation, buffering, etc).
