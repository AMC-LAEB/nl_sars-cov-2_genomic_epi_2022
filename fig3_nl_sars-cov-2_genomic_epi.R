library('seraphim')
library('diagram')
library('lubridate')
library('stringr')

# set local working directory 
#setwd('~/Dropbox/github_repo/covid19_nl_w2/new_analysis_20211010/nl_sars-cov-2_genomic_epi_2022/')
source("mccExtractions.r")

# 1. Extracting the spatio-temporal information embedded in posterior trees
Alpha_localTreesDirectory = "Alpha_Extracted_trees"
Alpha_mostRecentSamplingDatum = 2021.641096 # alpha 
#Alpha_allTrees = scan(file = "nl_sampled_alpha_cauchyRRW.trees", what="", sep='\n', quiet = T)

Delta_localTreesDirectory = "Delta_Extracted_trees"
Delta_mostRecentSamplingDatum = 2021.663014 # delta 
#Delta_allTrees = scan(file = "nl_sampled_delta_cauchyRWW.trees", what="", sep='\n', quiet = T)

nberOfTreesToSample = 1000

treeExtractionsWrapper = function(localTreesDirectory, allTrees, nberOfTreesToSample, mostRecentSamplingDatum){
  burnIn = 1000 # number of trees to discard as burn-in (NOT states)
  randomSampling = FALSE
  coordinateAttributeName = "location"
  nberOfCores = 8
  treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling,
                  nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, 
                  nberOfCores)
}

# launch tree extraction 
#treeExtractionsWrapper(Alpha_localTreesDirectory, Alpha_allTrees, nberOfTreesToSample, Alpha_mostRecentSamplingDatum)
#treeExtractionsWrapper(Delta_localTreesDirectory, Delta_allTrees, nberOfTreesToSample, Delta_mostRecentSamplingDatum)

# 2. Extracting the spatio-temporal information embedded in the MCC tree
Alpha_mcc_tree = readAnnotatedNexus("./trees/nl_sampled_alpha_cauchyRRW.mcc.tre")
Delta_mcc_tree = readAnnotatedNexus("./trees/nl_sampled_delta_cauchyRRW.mcc.tre")

Alpha_mcc_tab = mccExtractions(Alpha_mcc_tree, Alpha_mostRecentSamplingDatum)
Delta_mcc_tab = mccExtractions(Delta_mcc_tree, Delta_mostRecentSamplingDatum)

# 3. Estimating the HPD region for each time slice
# produces a list of distinct spatial polygon data frames, with one data frame for each time slice
nberOfExtractionFiles = nberOfTreesToSample
prob = 0.95 # probability corresponding to the HPD (highest posterior density) regions.
precision = 0.01916 # time interval used to define the successive time slices. (one week)

Alpha_startDatum = min(Alpha_mcc_tab[,"startYear"])
Alpha_polygons = suppressWarnings(spreadGraphic2(Alpha_localTreesDirectory, nberOfExtractionFiles, prob, Alpha_startDatum, precision))

Delta_startDatum = min(Delta_mcc_tab[,"startYear"]) 
Delta_polygons = suppressWarnings(spreadGraphic2(Delta_localTreesDirectory, nberOfExtractionFiles, prob, Delta_startDatum, precision))

##################

#4. Plotting the dispersal history of lineages

fig_label <- function(text, region="figure", pos="topleft", yadjust=1, cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  text(x1, y1 * yadjust, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

##########

png("../manuscript/figure3.png", width=8.3, height=11.8 * 0.65, units='in', res=330)
template_raster = shapefile("../data/NLD_adm_shp/NLD_adm1.shp")

# set date color scale 
colour_scale = colorRampPalette(brewer.pal(11,"Spectral"))(101)#[21:121]
minYear = 2020.7267759562842; maxYear = 2021.66301369863 # minYear and maxYear to study period 

# set plot layout 
layout(mat = matrix(c(2,3,4,5,6,7,
                      1,1,1,8,8,8), nrow = 2, ncol= 6, byrow=TRUE),
       heights = c(1, 2),
       widths = c(1, 1, 1, 1))

# plot main Alpha 
# set Alpha endYears_colours 
endYears_indices = (((Alpha_mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
# set Alpha polygon colours 
polygons_colours = rep(NA, length(Alpha_polygons))
for (i in 1:length(Alpha_polygons))
{
  date = as.numeric(names(Alpha_polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colour_scale[polygon_index],"40")
}

par(mar = c(0,0,0,0), xpd=NA) # set plot margins 
plot(template_raster, col="gray90", lwd=0.1,) # plot Netherlands map 
fig_label(text = 'b', region = 'plot', cex = 2, font = 2)
#mtext("b", adj=0, line=-1)

for (i in 1:length(Alpha_polygons))
{
  plot(Alpha_polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
} 
# plot arrows 
for (i in 1:dim(Alpha_mcc_tab)[1])
{
  curvedarrow(cbind(Alpha_mcc_tab[i,"startLon"],Alpha_mcc_tab[i,"startLat"]), cbind(Alpha_mcc_tab[i,"endLon"],Alpha_mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=.8, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.3, dr=NA, endhead=F)
}
# plot tips 
for (i in dim(Alpha_mcc_tab)[1]:1)
{
  #if (i == 1)
  #{
  #  points(Alpha_mcc_tab[i,"startLon"], Alpha_mcc_tab[i,"startLat"], pch=16, col=colour_scale[1], cex=0.8)
  #  points(Alpha_mcc_tab[i,"startLon"], Alpha_mcc_tab[i,"startLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)
  #}
  points(Alpha_mcc_tab[i,"endLon"], Alpha_mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
  points(Alpha_mcc_tab[i,"endLon"], Alpha_mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)
}

### Alpha - M1

par(mar = c(0,0,0,0)) # set plot margins 
plot(template_raster, col="gray90", lwd=0.1,) # plot Netherlands map 
fig_label(text = 'a', region = 'plot', cex = 2, font = 2)
#mtext("a", adj=0, line=1)

# first month 
cutoff_startdate <- 2020.915
cutoff_enddate <- 2020.973

# plot arrows 
for (i in 1:dim(Alpha_mcc_tab)[1])
{
  if (Alpha_mcc_tab[i,"startYear"] >= cutoff_startdate & Alpha_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    curvedarrow(cbind(Alpha_mcc_tab[i,"startLon"],Alpha_mcc_tab[i,"startLat"]), cbind(Alpha_mcc_tab[i,"endLon"],Alpha_mcc_tab[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=.8, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.3, dr=NA, endhead=F)  
  }
}
# plot tips 
for (i in dim(Alpha_mcc_tab)[1]:1)
{
  if (Alpha_mcc_tab[i,"startYear"] >= cutoff_startdate & Alpha_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    points(Alpha_mcc_tab[i,"endLon"], Alpha_mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.6)
    points(Alpha_mcc_tab[i,"endLon"], Alpha_mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)  
  }
}

fig_label("2-29 Dec 2020\n(W49-W52)", region = "plot", pos = "top", cex = 1.0, yadjust=0.995)

### Alpha - M1 onwards to w10

par(mar = c(0,0,0,0)) # set plot margins 
plot(template_raster, col="gray90", lwd=0.1,) # plot Netherlands map 
#mtext("a", adj=0, line=1)

# first month 
cutoff_startdate <- 2020.973
cutoff_enddate <- 2021.126

# plot arrows 
for (i in 1:dim(Alpha_mcc_tab)[1])
{
  if (Alpha_mcc_tab[i,"startYear"] >= cutoff_startdate & Alpha_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    curvedarrow(cbind(Alpha_mcc_tab[i,"startLon"],Alpha_mcc_tab[i,"startLat"]), cbind(Alpha_mcc_tab[i,"endLon"],Alpha_mcc_tab[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=.8, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.3, dr=NA, endhead=F)  
  }
}
# plot tips 
for (i in dim(Alpha_mcc_tab)[1]:1)
{
  if (Alpha_mcc_tab[i,"startYear"] >= cutoff_startdate & Alpha_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    points(Alpha_mcc_tab[i,"endLon"], Alpha_mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.6)
    points(Alpha_mcc_tab[i,"endLon"], Alpha_mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)  
  }
}

fig_label("30 Dec 2020-9 Mar 2021\n(W53-W9)", region = "plot", pos = "top", cex = 1.0, yadjust=0.995)

### Alpha - week 10 onwards 

par(mar = c(0,0,0,0)) # set plot margins 
plot(template_raster, col="gray90", lwd=0.1,) # plot Netherlands map 
#mtext("a", adj=0, line=1)

# first month 
cutoff_startdate <- 2021.126
cutoff_enddate <- 2021.663

# plot arrows 
for (i in 1:dim(Alpha_mcc_tab)[1])
{
  if (Alpha_mcc_tab[i,"startYear"] >= cutoff_startdate & Alpha_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    curvedarrow(cbind(Alpha_mcc_tab[i,"startLon"],Alpha_mcc_tab[i,"startLat"]), cbind(Alpha_mcc_tab[i,"endLon"],Alpha_mcc_tab[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=.8, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.3, dr=NA, endhead=F)  
  }
}
# plot tips 
for (i in dim(Alpha_mcc_tab)[1]:1)
{
  if (Alpha_mcc_tab[i,"startYear"] >= cutoff_startdate & Alpha_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    points(Alpha_mcc_tab[i,"endLon"], Alpha_mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.6)
    points(Alpha_mcc_tab[i,"endLon"], Alpha_mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)  
  }
}

fig_label("10 Mar-31 Aug 2021\n(W10-W34)", region = "plot", pos = "top", cex = 1.0, yadjust=0.995)

# set Delta endYears_colours 
endYears_indices = (((Delta_mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
polygons_colours = rep(NA, length(Delta_polygons))
# set Delta polygon colours 
for (i in 1:length(Delta_polygons))
{
  date = as.numeric(names(Delta_polygons[[i]]))
  if (date >= 2021.11){
    polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
    polygons_colours[i] = paste0(colour_scale[polygon_index],"40")
    print (polygon_index)
  }
}

# Delta - first month 
par(mar = c(0,0,0,0)) # set plot margins 
plot(template_raster, col="gray90", lwd=0.1,) # plot Netherlands map 
fig_label(text = 'c', region = 'plot', cex = 2, font = 2)
#mtext("c", adj=0, line=1)

j = 1
cutoff_startdate = 2021.0; cutoff_enddate = 2021.375

# plot arrows 
for (i in 1:dim(Delta_mcc_tab)[1])
{
  if (Delta_mcc_tab[i,"startYear"] >= cutoff_startdate & Delta_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    curvedarrow(cbind(Delta_mcc_tab[i,"startLon"],Delta_mcc_tab[i,"startLat"]), cbind(Delta_mcc_tab[i,"endLon"],Delta_mcc_tab[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=.8, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.3, dr=NA, endhead=F)
  }
}

# plot tips 
for (i in dim(Delta_mcc_tab)[1]:1)
{
  if (Delta_mcc_tab[i,"startYear"] >= cutoff_startdate & Delta_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    points(Delta_mcc_tab[i,"endLon"], Delta_mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.6)
    points(Delta_mcc_tab[i,"endLon"], Delta_mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)
  }
}

fig_label("20 Apr-18 May 2021\n(W16-W19)", region = "plot", pos = "top", cex = 1.0, yadjust=0.995)

# Delta - second month 
par(mar = c(0,0,0,0)) # set plot margins 
plot(template_raster, col="gray90", lwd=0.1,) # plot Netherlands map 

cutoff_startdate = 2021.375; cutoff_enddate = 2021.471

# plot arrows 
for (i in 1:dim(Delta_mcc_tab)[1])
{
  if (Delta_mcc_tab[i,"startYear"] >= cutoff_startdate & Delta_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    curvedarrow(cbind(Delta_mcc_tab[i,"startLon"],Delta_mcc_tab[i,"startLat"]), cbind(Delta_mcc_tab[i,"endLon"],Delta_mcc_tab[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=.8, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.3, dr=NA, endhead=F)
  }
}

# plot tips 
for (i in dim(Delta_mcc_tab)[1]:1)
{
  if (Delta_mcc_tab[i,"startYear"] >= cutoff_startdate & Delta_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    points(Delta_mcc_tab[i,"endLon"], Delta_mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.6)
    points(Delta_mcc_tab[i,"endLon"], Delta_mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)
  }
}

fig_label("19 May-21 Jun 2021\n(W20-W24)", region = "plot", pos = "top", cex = 1.0, yadjust=0.995)

# Delta - second month 
par(mar = c(0,0,0,0)) # set plot margins 
plot(template_raster, col="gray90", lwd=0.1,) # plot Netherlands map 

cutoff_startdate = 2021.471; cutoff_enddate = 2021.529

# plot arrows 
for (i in 1:dim(Delta_mcc_tab)[1])
{
  if (Delta_mcc_tab[i,"startYear"] >= cutoff_startdate & Delta_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    curvedarrow(cbind(Delta_mcc_tab[i,"startLon"],Delta_mcc_tab[i,"startLat"]), cbind(Delta_mcc_tab[i,"endLon"],Delta_mcc_tab[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=.8, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.3, dr=NA, endhead=F)
  }
}

# plot tips 
for (i in dim(Delta_mcc_tab)[1]:1)
{
  if (Delta_mcc_tab[i,"startYear"] >= cutoff_startdate & Delta_mcc_tab[i,"endYear"] <= cutoff_enddate) {
    points(Delta_mcc_tab[i,"endLon"], Delta_mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.6)
    points(Delta_mcc_tab[i,"endLon"], Delta_mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)
  }
}

fig_label("22 Jun-13 Jul 2021\n(W25-W27)", region = "plot", pos = "top", cex = 1.0, yadjust=0.995)

# plot main Delta 
par(mar = c(0,0,0,0)) # set plot margins 
plot(template_raster, col="gray90", lwd=0.1,) # plot Netherlands map 
fig_label(text = 'd', region = 'plot', cex = 2, font = 2)
#mtext("d", adj=0, line=1)

for (i in 1:length(Delta_polygons))
{
  date = as.numeric(names(Delta_polygons[[i]]))
  if (date >= 2021.11){
    plot(Delta_polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
  }
} 
# plot arrows 
for (i in 1:dim(Delta_mcc_tab)[1])
{
  if (Delta_mcc_tab[i,"startYear"] >= 2021.11 & Delta_mcc_tab[i,"endYear"] >= 2021.11){
    curvedarrow(cbind(Delta_mcc_tab[i,"startLon"],Delta_mcc_tab[i,"startLat"]), cbind(Delta_mcc_tab[i,"endLon"],Delta_mcc_tab[i,"endLat"]), arr.length=0,
                arr.width=0, lwd=.8, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.3, dr=NA, endhead=F)
  }
}

# plot tips 
for (i in dim(Delta_mcc_tab)[1]:1)
{
  if (Delta_mcc_tab[i,"startYear"] >= 2021.11 & Delta_mcc_tab[i,"endYear"] >= 2021.11){
    points(Delta_mcc_tab[i,"endLon"], Delta_mcc_tab[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.8)
    points(Delta_mcc_tab[i,"endLon"], Delta_mcc_tab[i,"endLat"], pch=1, col="gray10", cex=0.8, lwd=0.2)
  }
}

rast = raster(matrix(nrow=1, ncol=2)); 
rast[1] = 2020.743; rast[2] = 2021.663
plot(rast, legend.only=T, col=colour_scale, smallplot=c(0.02,0.4,0.05,0.07),
     legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.8, at=c(2020.743, 2020.896, 2021.049, 2021.203, 2021.356, 2021.51, 2021.663), 
                    labels=c('W39', 'W47', 'W2', 'W10', 'W18', 'W26', 'W34')))

dev.off()

write.csv(Alpha_mcc_tab, file = "./Alpha_mcc_tab.csv")
write.csv(Delta_mcc_tab, file = "./Delta_mcc_tab.csv")


