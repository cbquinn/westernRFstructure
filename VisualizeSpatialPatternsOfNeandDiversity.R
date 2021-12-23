
# Script to create Figure 5 in western red fox genetic structure paper
# Plot visualizes genetic diversity summary statistics and Ne as a spatially explicit genetic surface

# See Shirk A, Cushman S (2011). sGD: software for estimating spatially explicit indices of genetic diversity. Molecular Ecology Resources 11(5): 922-934
# https://github.com/Andrew-Shirk/sGD

# Estimates of genetic diversity were calculated using two approaches based on results of genetic structure analyses: (a) populations with strong genetic structure that were determined to be discrete were calculated in the traditional manner in which one summary statistic is estimated based on all genotypes located within a minimum convex polygon, or (b) populations with relatively continuous genetic structure, in which an estimate is calculated for every individual based on genotypes that fall within a pre-defined radius of its location. For the latter, we used the sGD package (Shirk & Cushman 2011). To create plots, we assigned every sample a genetic diversity estimate that pertained either to its discrete population (a) or its genetic neighborhood (b) and then used inverse weighting to interpolate a spatially explicit continuous surface of genetic diversity.  



# load libraries
library(sp)
library(sGD) 
library(tidyverse)
library(rgdal)
library(hierfstat)
library(adegenet)
library(raster)
library(rgeos)
library(RColorBrewer)
library(cowplot)


# Calls on functions in 3 sub-scripts"

# Slight modifications made to original sGD function in package because I kept having errors thrown (version incompatabilties?)
source("mysGDfunction.R")

# Helper script that writes genepop file from adegenet object. Used to create input for NeEstimator2
source("write_genepop_function.R")

# Function to interpolate with inverse distance weighting, using idw function in package gstat
source("myInterpolateFunction.R")

###########################
#### Load mapping data ####
###########################

# state outlines for ggplot
states <- map_data("state")
west <- subset(states, region %in% c("california", "oregon", "washington", "nevada", "idaho", "montana", "utah", "wyoming", "arizona", "new mexico", "colorado"))

# read in polygons
discrete_polys <- readOGR("MergedDiscrete_polygon.shp")

# set color palette
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"), space="Lab")

################################################
#### Load genotypes and lat/lon coordinates ####
################################################

# read in genotypes of "Discrete" populuations
d_discrete <- read.structure("discreteRFgenotypes.stru", n.ind = 301, n.loc = 31, onerowperind = T, col.lab = 1, col.pop = 2, col.others = 0, row.marknames = 1, NA.char = "-9", ask = TRUE, quiet = T)
# read in genotypes of "Continuous" populations
d_continuous <- read.structure("continuousRFgenotypes.stru", n.ind = 341, n.loc = 31, onerowperind = T, col.lab = 1, col.pop = 0, col.others = 0, row.marknames = 1, NA.char = "-9", ask = TRUE, quiet = T)
pop(d_continuous) <- rep("Intermountain", nInd(d_continuous))

# read in Lat/Lon coordinates of all samples
coord <- read_tsv("samplelocations_autosomalgenotypes.txt")

# convet to planar coordinates
splatlon <- SpatialPointsDataFrame(coord, coords = coord[, c("lon", "lat")], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
spAlbers <- spTransform(splatlon, CRS("+init=epsg:5070"))
sp_coord <- as.data.frame(spAlbers) %>%
  rename(albersX=lon.1, albersY=lat.1)

# IDs of individuals in discrete vs. continuous
discrete_ids <- row.names(as.data.frame(d_discrete))
continuous_ids <- row.names(as.data.frame(d_continuous))

# spatial df of coords - discrete and continuous
spcoord_discrete <- subset(spAlbers, sampleID %in% discrete_ids)
spcoord_continuous <- subset(spAlbers, sampleID %in% continuous_ids)

# nonspatial df of coords - discrete and continuous

coord_discrete <- sp_coord %>%
  filter(sampleID %in% row.names(as.data.frame(d_discrete))) 
coord_continuous <- sp_coord %>%
  filter(sampleID %in% row.names(as.data.frame(d_continuous))) 

# albers only (planar)
coord_discrete_albers <- sp_coord %>%
  filter(sampleID %in% row.names(as.data.frame(d_discrete))) %>%
  dplyr::select(sampleID, albersX, albersY)
coord_continuous_albers <- sp_coord %>%
  filter(sampleID %in% row.names(as.data.frame(d_continuous))) %>%
  dplyr::select(sampleID, albersX, albersY)

#############################################
#### Continuous (neighborhood) estimates ####
#############################################

# directory of NeEstimator program
NeEstimator_dir <- "/Applications/NeEstimator_64bit_191125"

# define neighborhoold
radius <- 85000 # units are in meters 
# minimum neighborhood size (no. of genotypes)
min_N <- 10

# create pairwise distance matrix (units in meters)
IBD_dist <- distmat(spcoord_continuous, method="ed")

# calculate using sGD function from package sGD
# (modified sGD function slightly because errors were thrown due to package incompatabilities, but no significant changes)
# produces a dataframe with summary statistics for each sample location based on the genotypes in its radius-defined neighborhood
IBD_continuous <- mysGD(genind_obj=d_continuous, xy=coord_continuous_albers, dist.mat=IBD_dist, NH_radius=radius, min_N=min_N, metrics=c("GD", "NS"), NeEstimator_dir=NeEstimator_dir)

# combine with lat/lon data
data_cont <- bind_cols(IBD_continuous, lat=coord_continuous$lat, lon=coord_continuous$lon) %>%
  mutate(Ne=NS_ex0.05) %>%
  rename(sampleID=NH_ID, albersX=X, albersY=Y)

##########################################
#### Discrete (traditional) estimates ####
##########################################

# function to calculate standard error
se <- function(x){
  sd(x)/(length(x)^0.5)
} 

# For all "discrete" populations calculate summary statistics in traditional manner

# Allelic Richness in hierfstat 
ar <- allelic.richness(d_discrete) # per locus
ar$min.all # the number of alleles used for rarefaction
ARmean <- apply(ar$Ar, 2, mean) # mean per population
ARse <- apply(ar$Ar, 2, se)

# Other basic stats
basicstats <- basic.stats(d_discrete)
Homean <- apply(basicstats$Ho,2, mean)
Hose <- apply(basicstats$Ho,2, se)
Hsmean <- apply(basicstats$Hs, 2, mean)
Hsse <- apply(basicstats$Hs, 2, se)
Nind <- apply(basicstats$n.ind.samp, 2, max, na.rm=TRUE)
FIS <- apply(basicstats$Fis, 2, mean)

# massage into same format as output for sGD neighborhood-based estimates
data_discrete <- bind_cols(sampleID=indNames(d_discrete), pop=pop(d_discrete)) %>%
  left_join(coord_discrete) %>%
  left_join(bind_cols(pop=names(Homean), N=Nind, Ar=ARmean, Ho=Homean, Hs=Hsmean ))

# Estimate Ne for discrete populations using LD method in NeEstimator
write.genepop(d_discrete, "discretepops", "populations with discrete population structure in the western US")

# Ran NeEstimator2 through GUI, using LD method with random mating and 0.05 pcrit
# Read in output after deleting all headers and only keeping the tab delimited output
t <- read_tsv("discretepops_genepopLDxLD_edited.txt", col_names = c("pop_firstsamp", "n", "mean_weighted", "ind_alleles", "r2", "r2_exp", "Ne", "CIlow_parametric", "CIhigh_parametric", "CIlow_jackknife", "CIhigh_jackknife", "df_eff"))

# match NeEstimator output ID with popID
# extract sampleID from pop-identifier
t$sampleID <- sub(".*:", "", t$pop_firstsamp)  
# and match to pop
t$pop <- data_discrete$pop[match(t$sampleID, data_discrete$sampleID)]

# add to data frame of diversity results
# final df is in same format as sGD output, except every sample within a discrete population has an identical estimate 
data_discrete <- t %>%
  dplyr::select(pop, Ne) %>%
  right_join(data_discrete)

##########################################
##### Combine discrete and continuous ####
##########################################

data <- bind_rows(data_cont, data_discrete)

######################
##### Interpolate ####
######################

# default extent of study system
myextent <- c(left = -124, bottom = 36, right = -108, top = 48)

# execute interpolation for each metric
idw.Ne <- myInterPolate(df=data, var="Ne", lat="lat", lon="lon", myextent=myextent)
idw.Ar <- myInterPolate(df=data, var="Ar", lat="lat", lon="lon", myextent=myextent)
idw.He <- myInterPolate(df=data, var="Hs", lat="lat", lon="lon", myextent=myextent)
idw.Ho <- myInterPolate(df=data, var="Ho", lat="lat", lon="lon", myextent=myextent)
idw.FIS <- myInterPolate(df=data, var="FIS", lat="lat", lon="lon", myextent=myextent)
idw.Nsamp <- myInterPolate(df=data, var="N", lat="lat", lon="lon")

# interpolated surface resorts to global average in areas far away from where the data are, which is not informative
# clip surface to retain only area within 100 km of a data point where there is a diversity estimate

r.Ne <- BufferGrid(idw.Ne)
r.Ar <- BufferGrid(idw.Ar)
r.He <- BufferGrid(idw.He)
r.Ho <- BufferGrid(idw.Ho)
r.FIS <- BufferGrid(idw.FIS)

###############
#### Plots ####
###############

legtextsz <- 6
legtitsz <- 6
axistxtsz <- 6
ptsize <- 1.25

#### Ne ####

# adjust color scale (want the most contrast when Ne < 100)
thismax <- max(data$Ne, na.rm = TRUE)
mybreaks <- c(0,10,20,30,50,100,120, round(thismax,0)-1)
mylabels <- c(0,10,"","","50",100,"", round(thismax,0)-1)
mylabels[8] <- 350
myscale <- mybreaks/thismax

plot.Ne <- ggplot() +
  geom_tile(data = r.Ne, alpha=0.7, aes(x = x, y = y, fill=var1.pred)) +
  geom_polygon(data = discrete_polys, aes(x=long, y = lat, group = group), 
               color="gray35", fill="transparent", linetype="dashed", size=0.75) + 
  geom_polygon(data = west, aes(x=long, y = lat, group = group), 
               color="black", fill="transparent") + 
  geom_point(data=data, aes(x=lon, y=lat, fill=Ne), shape=21, size=ptsize) + 
  coord_fixed(1.3) +
  scale_fill_gradientn(colors=myPalette(100), na.value="white", name=expression(N[e]),
                    values=myscale, breaks=mybreaks, labels=mylabels) +
  scale_color_manual() +
  coord_cartesian(ylim = c(34, 47), xlim=c(-123, -109)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=axistxtsz),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(.85, .23), 
    legend.text = element_text(size=legtextsz),
    legend.title = element_text(size=legtitsz),
    aspect.ratio=1,
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(3.5, "mm"))
plot.Ne

#### He ####

plot.He <- ggplot() +
  geom_tile(data = r.He, alpha=0.7, aes(x = x, y = y, fill=var1.pred)) +
  geom_polygon(data = discrete_polys, aes(x=long, y = lat, group = group), 
               color="gray35", fill="transparent", linetype="dashed", size=0.75) + 
  geom_polygon(data = west, aes(x=long, y = lat, group = group), 
               color="black", fill="transparent") + 
  geom_point(data=data, aes(x=lon, y=lat, fill=Hs), shape=21, size=ptsize) + 
  coord_fixed(1.3) +
  scale_fill_gradientn(colors=myPalette(100), na.value="white", name=expression(H[E])) +
  scale_color_manual() +
  coord_cartesian(ylim = c(34, 47), xlim=c(-123, -109)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=axistxtsz),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(.85, .23), 
    legend.text = element_text(size=legtextsz),
    legend.title = element_text(size=legtitsz),
    aspect.ratio=1,
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(3.5, "mm"))
plot.He

#### Ho ####

plot.Ho <- ggplot() +
  geom_tile(data = r.Ho, alpha=0.7, aes(x = x, y = y, fill=var1.pred)) +
  geom_polygon(data = discrete_polys, aes(x=long, y = lat, group = group), 
               color="gray35", fill="transparent", linetype="dashed", size=0.75) + 
  geom_polygon(data = west, aes(x=long, y = lat, group = group), 
               color="black", fill="transparent") + 
  geom_point(data=data, aes(x=lon, y=lat, fill=Ho), shape=21, size=ptsize) + 
  coord_fixed(1.3) +
  scale_fill_gradientn(colors=myPalette(100), na.value="white", name=expression(H[O]))+
  scale_color_manual() +
  coord_cartesian(ylim = c(34, 47), xlim=c(-123, -109)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=axistxtsz),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(.85, .23), 
    legend.text = element_text(size=legtextsz),
    legend.title = element_text(size=legtitsz),
    aspect.ratio=1,
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(3.5, "mm"))
plot.Ho

#### AR ####

plot.Ar <- ggplot() +
  geom_tile(data = r.Ar, alpha=0.7, aes(x = x, y = y, fill=var1.pred)) +
  geom_polygon(data = discrete_polys, aes(x=long, y = lat, group = group), 
               color="gray35", fill="transparent", linetype="dashed", size=0.75) + 
  geom_polygon(data = west, aes(x=long, y = lat, group = group), 
               color="black", fill="transparent") + 
  geom_point(data=data, aes(x=lon, y=lat, fill=Ar), shape=21, size=ptsize) + 
  coord_fixed(1.3) +
  scale_fill_gradientn(colors=myPalette(100), na.value="white", name="AR") +
  scale_color_manual() +
  coord_cartesian(ylim = c(34, 47), xlim=c(-123, -109)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=axistxtsz),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(.85, .23), 
    legend.text = element_text(size=legtextsz),
    legend.title = element_text(size=legtitsz),
    aspect.ratio=1,
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(3.5, "mm"))
  plot.Ar

#### mtDNA & Ysats ####

source("plotNeighborhoodmtDNADiversity.R")
# creates mtDNA figure called plot.mt
source("plotNeighborhoodYsatDiversity.R")
# creates Ysats diversity figure called plot.y

# adjust plotting features to match other plots
plot.mt <- plot.mt + 
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=axistxtsz),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(.85, .23), 
    legend.text = element_text(size=legtextsz),
    legend.title = element_text(size=legtitsz),
    aspect.ratio=1,
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(3.5, "mm"))
plot.y <- plot.y +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size=axistxtsz),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(.85, .23), 
    legend.text = element_text(size=legtextsz),
    legend.title = element_text(size=legtitsz),
    aspect.ratio=1,
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(3.5, "mm"),)

#######################
#### Arrange plots ####
#######################

#### main figure (Fig. 5 ###

divfig <- plot_grid(plot.mt, plot.y, plot.He, plot.Ne, 
          labels="AUTO")
divfig

# 175 mm = 2 col = 6.89 inches
save_plot("spatialDiversity_vvmc_combined.pdf",
  divfig,
  base_height = 6.89,
  base_width = 6.89,
)

#### SI plot ####

sifig <- plot_grid(plot.Ho, plot.Ar,
                   labels="AUTO")
sifig
save_plot("spatialDiversity_SI.pdf",
          sifig,
          base_height = 6.89/2,
          base_width = 6.89,
)
