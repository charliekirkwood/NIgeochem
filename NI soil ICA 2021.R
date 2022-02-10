library(data.table)
library(ggplot2)
library(ica)
library(compositions)
library(Rcpp)
library(raster)
library(ranger)

# data wrangling ####
# load NI soil xrf data as pointdata
pointdata <- fread("data/regionalsoilsaxrf.csv")
names(pointdata)[2:3] <- c("Easting", "Northing")
pointdata

# exclude pH (not part of composition) and LOI 
# (we'll just use unmeasured remainder as a proxy, to avoid having two 'LOI's)
pointdata <- pointdata[,1:55]

# Exclude NA values (they are not compositionally valid)
pointdata <- na.omit(pointdata)

# Convert oxide percentages to parts per million, like all other elements
pointdata[, names(pointdata)[grepl( "O" , names(pointdata))] := lapply(.SD, function(x) x*10000), .SDcols = grepl( "O" , names(pointdata))]
pointdata

# Set zeros to half limit of detection
pointdata[, 4:55 := lapply(.SD, function(x) replace(x, x == 0, min(x[x > 0], na.rm = TRUE) / 2)), .SDcols = c(4:55)]
pointdata

# Create 'R' component, which is the unmeasured remainder
pointdata[, R := 1000000-rowSums(pointdata[,4:55]), ]

# Create ILR matrix for the element selection (all of them in this case) ####
# choose elements to include 
elementsel <- names(pointdata)[-(1:3)]

# create isometric log-ratio matrix
pointdata.ilr <- ilr(acomp(pointdata[,..elementsel]))
pointdata.ilr

# run ICA on ILR matrix (to find independent compositional components)
dataPCA <- icafast(pointdata.ilr, nc=8, alg="def", fun="exp", center = TRUE)
pairs(dataPCA$S, col=rgb(0,0,0,5,maxColorValue=255), upper.panel = NULL, labels=c(paste("Comp.", seq(1,8,1))), cex.labels = 1.2)

# plot proportion of variance explained by each component
ggplot(data = data.frame(PC = 1:length(dataPCA$vafs), propVariance = dataPCA$vafs*100), aes(x= PC, y=propVariance)) + geom_point() + scale_y_continuous(limits=c(0,max(dataPCA$vafs*100*1.01))) + 
  xlab("Component") + ylab("Proportion of Variance (%)") + theme_bw()

# plot a map of the observations, coloured by chosen component (Vx)
pointdataica <- cbind(pointdata, dataPCA$Y)
ggplot(data = pointdataica) + geom_point(aes(x = Easting, y = Northing, col = V1)) + theme_bw()

# gridded maps are better, so a quick random forest for interpolation ####
# first load NI elevation grid, which will double as our interpolation grid too
grid <- fread("data/OSNI_OpenData_50m_DTM.csv")
grid

gridrast <- rasterFromXYZ(grid)
gridrastcoarse <- raster::aggregate(gridrast, fact = 5)
length(gridrastcoarse)

pointdataica[, Elevation := raster::extract(gridrast, pointdata[, c("Easting", "Northing")], method = "bilinear"), ]
pointdataica <- na.omit(pointdataica)

predgrid <- as.data.table(raster::as.data.frame(gridrastcoarse, xy = TRUE, na.rm = TRUE))
names(predgrid) <- c("Easting", "Northing", "Elevation")

fidelity <- 6
for(a in seq(0,90,fidelity)){
  predgrid[, paste0("N_", a) := Easting*cos(a) + Northing*sin(a)]
  predgrid[, paste0("E_", a) := -Easting*sin(a) + Northing*cos(a)]
  
  pointdataica[, paste0("N_", a) := Easting*cos(a) + Northing*sin(a)]
  pointdataica[, paste0("E_", a) := -Easting*sin(a) + Northing*cos(a)]
}

models <- list()
for(i in paste0("V",1:8)){
  models[[i]] <- ranger(as.formula(paste0(i," ~ Elevation + ", paste0("E_", seq(0,90,fidelity), collapse=" + "), " + ", paste0("N_", seq(0,90,fidelity), collapse=" + "))), 
                        data = pointdataica, importance = "none", num.trees = 250, sample.fraction = 0.5, min.node.size = 15, replace = FALSE)
  predgrid[, paste0(i) := predict(models[[i]], predgrid[, -(1:2)])$predictions, ]
}

# List model correlation to out-of-bag data
unlist(lapply(models, '[[', "r.squared")) 

comp <- "V5"
for(comp in paste0("V",1:8)){
ggplot(predgrid) + geom_raster(aes_string(x = "Easting", y = "Northing", fill = comp)) + coord_equal() + theme_bw() +
  geom_point(data = pointdataica, aes_string(x = "Easting", y = "Northing", col = comp), shape = 16, size = 0.5) + theme(legend.position = c(0.05,0.8)) +
  scale_color_continuous(limits = c(min(c(predgrid[, get(comp)], pointdataica[, get(comp)])), max(c(predgrid[, get(comp)], pointdataica[, get(comp)])))) +
  scale_fill_continuous(limits = c(min(c(predgrid[, get(comp)], pointdataica[, get(comp)])), max(c(predgrid[, get(comp)], pointdataica[, get(comp)]))))
  ggsave(paste0("NI_soil_XRF_ICA_comp_", comp, ".png"), width = 180, height = 140, type = "cairo", units = "mm", scale = 1.25)
}