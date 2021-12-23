# Visualize mitochondrial diversity as a spatially explicit surface, following methods used for autosomal microsatellite genotypes

library(sp)
library(sGD) 
library(tidyverse)
library(rgdal)
library(RColorBrewer)
library(readxl)
library(raster)
library(rgeos)

source("geneDiversity.R") # has modified sgd function for mitochondrial DNA
source("myInterpolateFunction.R")

# read in mtDNA haplotype data
haplos <- read_tsv("mitohaplos_cytb_dloop_vvmc.tsv")

##########################################
#### Discrete (traditional) estimates #### 
##########################################

# grps of discrete populations
mydiscrete <- haplos %>%
  filter(discrete=="yes") %>%
  distinct(grp) %>%
  pull(grp)

# calculate gene diversity for each discrete grp
gene_D <- rep(NA,length(mydiscrete)) #empty vector
names(gene_D) <- mydiscrete
for (i in mydiscrete) {
  df <- haplos %>%
    group_by(grp) %>%
    filter(grp == i) %>%
    summarise(geneD=geneDiversity(CAT_cytbdloop)) 
  gene_D[i] <- df$geneD
}

# calculate variance for gene diversity
gene_SE <- rep(NA,length(mydiscrete)) 
names(gene_SE) <- mydiscrete
for (i in mydiscrete) {
  df <- haplos %>%
    group_by(grp) %>%
    filter(grp == i) %>%
    summarise(geneD=geneDiversity_se(CAT_cytbdloop)) 
  gene_SE[i] <- df$geneD
}

mtDNA_geneD_discrete <- bind_rows(gene_D, gene_SE)
mtDNA_geneD_discrete$estimate <- c("geneD", "SE")

write_tsv(mtDNA_geneD_discrete, "mtDNA_geneD_discrete.txt")

#############################################
#### Continuous (neighborhood) estimates #### 
#############################################

continuous <- haplos %>%
  filter(discrete == "no") 

mtHaplos <- continuous %>% dplyr::select(one_of("sampleID", "CAT_cytbdloop"))

coord <- dplyr::select(continuous, one_of(c("lon", "lat")))
#convert to USA_Contiguous_Albers_Equal_Area_Conic_USGS_version
splatlon <- SpatialPoints(coords = coord, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
spAlbers <- spTransform(splatlon, CRS("+init=epsg:5070"))

# calculate distance  matrix
IBD_dist <- distmat(spAlbers,method="ed")

radius <- 85000 # units are in meters 
min_N <- 10

results <- mysGD_mtDNA(mtHaplos=mtHaplos, xy=as.data.frame(spAlbers), dist.mat=IBD_dist, NH_radius=radius, min_N=min_N) 

###################################################
#### Combine Discrete and Continuous estimates ####
###################################################

# join neighborhood estimates with lat/lon data
mydf <- full_join(results, haplos, by=c("NH_ID"="sampleID"))

# fill in N and geneD for the discrete populations
for (i in mydiscrete) {
  mydf$geneD[mydf$grp==i] <- gene_D[i]
  mydf$N[mydf$grp==i] <- sum(mydf$grp==i)
}

##############
#### PLOT ####
##############

states <- map_data("state")
west <- subset(states, region %in% c("california", "oregon", "washington", "nevada", "idaho", "montana", "utah", "wyoming", "arizona", "new mexico", "colorado"))

# set color palette
myPalette <- colorRampPalette(brewer.pal(11, "Spectral"), space="Lab")

# interpolate
idw.output <- myInterPolate(df=mydf, var="geneD", lat="lat", lon="lon", myextent=myextent)

# read in polygons
discrete_polys <- readOGR("MergedDiscrete_polygon.shp")

# convert grid to raster
r <- rasterFromXYZ(idw.output, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# buffer around points (using albers)
# but first, remove points that didn't have enough samples nearby for a neighborhood (NAs)
df <- haplos %>%
  filter(sampleID %in% (mydf$NH_ID[which(!is.na(mydf$geneD))]))
coord <- dplyr::select(df, one_of(c("lon", "lat")))
#convert to USA_Contiguous_Albers_Equal_Area_Conic_USGS_version
splatlon <- SpatialPoints(coords = coord, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
spAlbers <- spTransform(splatlon, CRS("+init=epsg:5070"))
mybuffers_albers <- gBuffer(spAlbers, width=100000)
# convert back to lat/lon
mybuffers_wgs84 <- spTransform(mybuffers_albers, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# intersect with raster
rbuff <- mask(r, mybuffers_wgs84)

# make it plottable in ggplot2
rbuff_spdf <- as(rbuff, "SpatialPixelsDataFrame")
rbuff_df <- as.data.frame(rbuff_spdf)
names(rbuff_df) <- c("var1.pred", "x", "y")

# plot 
plot.mt <- ggplot() +
  geom_tile(data = rbuff_df, alpha=0.7, aes(x = x, y = y, fill=var1.pred)) +
  geom_polygon(data = discrete_polys, aes(x=long, y = lat, group = group), 
               color="gray35", fill="transparent", linetype="dashed", size=0.75) + 
  geom_polygon(data = west, aes(x=long, y = lat, group = group), 
               color="black", fill="transparent") + 
  geom_point(data=mydf, aes(x=lon, y=lat, fill=geneD), shape=21, size=1.25) + 
  coord_fixed(1.3) +
  scale_fill_gradientn(colors=myPalette(100), na.value="white", name="mtDNA\nDiversity") +
  coord_cartesian(ylim = c(34, 47), xlim=c(-123, -109)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(.85, .23), 
    aspect.ratio=1)
plot.mt


