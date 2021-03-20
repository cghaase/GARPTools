## ----message = FALSE, warning = FALSE, comment = "R>"--------------------
library(GARPTools)

## ----comment = "R>"------------------------------------------------------
data("wtdeer_df")
head(wtdeer_df)

## ----comment = "R>"------------------------------------------------------
data("wtdeer_locations")
wtdeer_locations

## ----eval = TRUE, echo = FALSE-------------------------------------------
library(raster)
data("env_layers")
writeRaster(env_layers[[1]], "C:/GARP/rasters/alt.asc", overwrite=TRUE)
writeRaster(env_layers[[2]], "C:/GARP/rasters/bio_1.asc", overwrite=TRUE)
writeRaster(env_layers[[3]], "C:/GARP/rasters/bio_12.asc", overwrite=TRUE)


## ----eval = FALSE, echo = TRUE-------------------------------------------
#  files <- list.files(path=paste(system.file(package="GARPTools"),"/data", sep = ""),
#                      pattern = ".asc", full.names = TRUE)
#  
#  env_layers <- stack(files)
#  
#  writeRaster(env_layers[[1]], "C:/GARP/rasters/alt.asc", overwrite=TRUE)
#  writeRaster(env_layers[[2]], "C:/GARP/rasters/bio_1.asc", overwrite=TRUE)
#  writeRaster(env_layers[[3]], "C:/GARP/rasters/bio_12.asc", overwrite=TRUE)

## ----comment = "R>", eval = FALSE, fig.height=2, fig.align = "center", fig.cap = "Three example environmental layers used in DesktopGARP experimental runs."----
#  plot(env_layers, nc = 3)

## ----comment = "R>", echo = -1, fig.width=4, fig.height=3, fig.align = "center", fig.cap = "Altitude of sampling area with sampling points."----
par(mar=c(3,2,2,1))
data("nc_boundary")
data("alt")


## ---- eval = TRUE, echo = FALSE------------------------------------------
library(SDMTools)

## ----comment = "R>", warning = FALSE, message = FALSE, results = "hide"----
rasterPrep(file.path = "C:/GARP/rasters", mask = nc_boundary, 
           output.path = "C:/GARP/resampled/")

## ----comment = "R>"------------------------------------------------------
data("wtdeer_locations")
alt_resample <- raster("C:/GARP/resampled/alt_resample.asc")

wtdeer_centroids <- centroid(x = alt_resample, points = wtdeer_locations, 
                             xy = wtdeer_df[,c("Latitude","Longitude")],
                             species = wtdeer_df$Species)

## ----comment = "R>", echo = -1, fig.width=4, fig.height=3, fig.align = "center", fig.cap = "Resampled altitude map of sampling area with sampling points (black) and centroid points (purple)."----
par(mar=c(2,2,2,1))
plot(alt_resample)
points(wtdeer_locations, pch = 16)
points(wtdeer_centroids, col = "purple")

## ---- eval=TRUE, echo=FALSE----------------------------------------------
library(rJava)
library(xlsx)

## ----comment = "R>"------------------------------------------------------
splitData(points = wtdeer_centroids, p = 0.75, type = "all", iterations = 10, 
          output = TRUE, output.path = "C:/GARP/wtdeer")

## ----eval = TRUE, echo = FALSE-------------------------------------------
library(raster)
data(tasks)

writeRaster(tasks[[1]], filename = "task1.asc", format = "ascii", overwrite=TRUE)
writeRaster(tasks[[2]], "C:/GARP/BestSubsets/task2.asc",overwrite=TRUE)
writeRaster(tasks[[3]], "C:/GARP/BestSubsets/task3.asc", overwrite=TRUE)
writeRaster(tasks[[4]], "C:/GARP/BestSubsets/task4.asc", overwrite=TRUE)
writeRaster(tasks[[5]], "C:/GARP/BestSubsets/task5.asc", overwrite=TRUE)
writeRaster(tasks[[6]], "C:/GARP/BestSubsets/task6.asc", overwrite=TRUE)
writeRaster(tasks[[7]], "C:/GARP/BestSubsets/task7.asc", overwrite=TRUE)
writeRaster(tasks[[8]], "C:/GARP/BestSubsets/task8.asc", overwrite=TRUE)
writeRaster(tasks[[9]], "C:/GARP/BestSubsets/task9.asc", overwrite=TRUE)
writeRaster(tasks[[10]], "C:/GARP/BestSubsets/task10.asc", overwrite=TRUE)

## ----eval = FALSE, echo = TRUE-------------------------------------------
#  files <- list.files(path=paste(system.file(package="GARPTools"),"/data", sep = ""),
#                      pattern = ".asc", full.names = TRUE)
#  
#  tasks <- stack(files)
#  
#  writeRaster(tasks[[1]], "C:/GARP/BestSubsets/task1.asc")
#  writeRaster(tasks[[2]], "C:/GARP/BestSubsets/task2.asc")
#  writeRaster(tasks[[3]], "C:/GARP/BestSubsets/task3.asc")
#  writeRaster(tasks[[4]], "C:/GARP/BestSubsets/task4.asc")
#  writeRaster(tasks[[5]], "C:/GARP/BestSubsets/task5.asc")
#  writeRaster(tasks[[6]], "C:/GARP/BestSubsets/task6.asc")
#  writeRaster(tasks[[7]], "C:/GARP/BestSubsets/task7.asc")
#  writeRaster(tasks[[8]], "C:/GARP/BestSubsets/task8.asc")
#  writeRaster(tasks[[9]], "C:/GARP/BestSubsets/task9.asc")
#  writeRaster(tasks[[10]], "C:/GARP/BestSubsets/task10.asc")

## ----comment = "R>"------------------------------------------------------
r.sum <- sumRasters("C:/GARP/BestSubsets/")

## ----comment = "R>"------------------------------------------------------
summary(r.sum)

## ----comment = "R>"------------------------------------------------------
freq(r.sum)

## ----eval = TRUE, echo = FALSE-------------------------------------------
library(maptools)

## ----echo = -1, message = FALSE, warning = FALSE, fig.width=4, fig.height=3, fig.align = "center", fig.cap = "Best subset of models output by DesktopGARP with testing points."----
par(mar=c(2,2,2,1))
test.pts <- readShapePoints("C:/GARP/wtdeer1_test.shp")
plot(r.sum)
plot(nc_boundary, add = TRUE)
points(test.pts, pch = 16)

## ----echo = -1, comment = "R>", fig.width=4, fig.height=4, fig.align = "center", fig.cap = "Plot of the Receiver Operating Characteristic (ROC) curve."----
par(mar=c(2,2,2,1))
aucGARP(n = 10, x = r.sum, points = test.pts)

## ----comment = "R>"------------------------------------------------------
commissionGARP(n = 10, x = r.sum, points = test.pts, file.path = "C:/GARP/BestSubsets/")
omissionGARP(n = 10, x = r.sum, points = test.pts)

## ----eval = FALSE, comment = "R>"----------------------------------------
#  modelRule <- FindModRule(pathA = "C:/GARP/runs/", pathB = "C:/GARP/BestSubsets/")

## ----eval = FALSE, comment = "R>"----------------------------------------
#  rangelogit <- ExtractRules("C:/GARP/RuleSets.txt", table = modelRule,
#                             pathA = "C:/GARP/resampled/", pathB = "C:/GARP/runs/",
#                             project = FALSE)

## ----eval = FALSE, comment = "R>"----------------------------------------
#  write.csv(rangelogit, "C:/GARP/RangeLogit.csv")

## ----eval = FALSE, comment = "R>"----------------------------------------
#  MedMinMax <- getMedMinMax(pathA = "C:/GARP/RangeLogit.csv",
#                            pathB = "C:/GARP/resampled/")

## ----eval = FALSE, comment = "R>"----------------------------------------
#  recaledM <- Rescale(MedMinMax)

## ----eval = FALSE, comment = "R>"----------------------------------------
#  plotRange(rescaledM, colour = "purple", n=15,
#            xlabel = "Environmental Varibles",
#            ylabel = "Median Range")
#  

## ----eval = FALSE, comment = "R>"----------------------------------------
#  PreNo <- TotPresRules(file = "C:/GARP/RuleSets.txt", table = modelRule,
#                        pathA = "C:/GARP/resampled/", pathB = "C:/GARP/runs/",
#                        project = FALSE)

## ----eval = FALSE, comment = "R>"----------------------------------------
#  prevalence <- Prevalence(pathA = "C:/GARP/RuleSets.txt",
#                           pathB = "C:/GARP/runs/", num = 10)

## ----eval = FALSE, comment = "R>"----------------------------------------
#  write.csv(prevalence, "C:/GARP/Prevalence.csv")

## ----eval = FALSE, comment = "R>"----------------------------------------
#  UnimportIndex <- UnimportIdx(prevalence, rescaled)

## ----eval = FALSE, comment = "R>"----------------------------------------
#  write.csv(UnimportanceIndex,"C:/GARP/UnimportanceIndex.csv")

