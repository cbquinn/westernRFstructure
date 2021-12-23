# Two helper functions used to create Figure 5 (visualization of genetic diversity in western red foxes) and SI figures

# default extent of study system
myextent <- c(left = -124, bottom = 36, right = -108, top = 48)

# Interpolate per-sample point estimates of genetic diversity to a continuous surface (uses idw function from gstat)
myInterPolate <- function(df, var, lat="lat", lon="lon",myextent) {
  require(gstat)
  # set up data
  mycol <- c(var, lat, lon)
  tempdata <- df %>%
    dplyr::select(one_of(mycol))
  tempdata <- na.omit(tempdata)
  tempdata$x <- tempdata$lon
  tempdata$y <- tempdata$lat
  coordinates(tempdata) <- ~x + y
  # create a blank raster to extent of system
  grd <- expand.grid(x=seq(from=myextent[1], myextent[3], 0.1),
                     y=seq(from=32, myextent[4], 0.1)) # 
  coordinates(grd) <- ~x + y
  # inverse distance weighting
  idw <- idw(formula= get(var) ~ 1, locations=tempdata, newdata=grd )
  idw.output = as.data.frame(idw)
  names(idw.output)[1:3] <- c("long", "lat", "var1.pred") # give names to the modeled variables
  idw.output <- dplyr::select(idw.output, -var1.var)
  return(idw.output)
}

# function to clip csurface to retain only area within 100 km of a data point where there is a diversity estimate
# (function assumes existence of a lot of objects in the workspace)

BufferGrid <- function(idw) {
  # grid to raster
  r <- rasterFromXYZ(idw, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  # remove points with NAs
  # remove NAs that don't have a large enough sample size
  buffercoords <- data %>%
    filter(!is.na(Ne)) %>%
    dplyr::select(albersX, albersY) 
  spbuffercoords <- SpatialPoints(coords = buffercoords, proj4string = CRS("+init=epsg:5070"))
  # buffer
  mybuffers_albers <- gBuffer(spbuffercoords, width=100000)
  # convert back to lat/lon
  mybuffers_wgs84 <- spTransform(mybuffers_albers, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  # intersect with raster
  rbuff <- mask(r, mybuffers_wgs84)
  # make it plottable in ggplot2
  rbuff_spdf <- as(rbuff, "SpatialPixelsDataFrame")
  rbuff_df <- as.data.frame(rbuff_spdf)
  names(rbuff_df) <- c("var1.pred", "x", "y")
  rbuff_df
}
