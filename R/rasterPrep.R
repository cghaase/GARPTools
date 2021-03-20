#' Prepares input rasters for GARP
#'
#' \code{rasterPrep} The desktop version of GARP requires all input environmental layers to be of the same spatial resolution and extent and saved as an ".asc". This function resamples input rasters to the same cell size and extent of a sample raster.
#' User can define sample raster in \code{res.raster} or define the number of cells in each row of the raster with \code{cells}. If both \code{res.raster} and \code{cell} parameters are \code{NULL},
#' raster with the largest cell size will be automatically defined as the sample raster. The user can also crop rasters to a mask of the study area with \code{mask = "mask filename"}.

#' @param file.path location of rasters to resample and crop
#' @param res.raster optional; a raster object to resample and crop all other rasters to
#' @param cells optional; a numeric value representing the number of cells per row to resample rasters
#' @param mask optional; a spatial object specifying extent to crop rasters
#' @param rescale logical; should each raster be centered on the mean and scaled by the standard deviation?
#' @param output.path a file path specifying where the output rasters should be saved
#'
#' @return For each raster within \code{file.path}, a new ACII file, resampled to the spatial resolution of the sample raster and spatial extent of the sample raster (\code{res.raster}) or \code{mask} file,
#' will be returned to \code{output.path} with the same name and \code{"_resample"} appended to the file name.
#'
#' @import raster
#' @import sp
#'
#' @export


rasterPrep <- function(file.path, res.raster = NULL, cells = NULL, mask = NULL, rescale = FALSE, output.path){

  #List all rasters in folder and create table
  list.rasters <- list.files(path = paste(file.path), full.names=TRUE)

  #Remove info folder
  for(i in 1:length(list.rasters)){
    if("info" %in% strsplit(list.rasters[i],split="/")[[1]][length(strsplit(list.rasters[i], split = "/")[[1]])]){list.rasters <- list.rasters[-i]}
  }

  raster.table <- data.frame(list.rasters)
  colnames(raster.table) <- "RasterName"

  for(i in 1:length(list.rasters)){
    r <- raster(list.rasters[i])
    raster.table$Projection[i]  <- projection(raster(list.rasters[i]))
    raster.table$ResolutionX[i] <- res(raster(list.rasters[i]))[1]
    raster.table$ResolutionY[i] <- res(raster(list.rasters[i]))[2]
  }

  #Print table
  print(raster.table)

  #Pull raster with largest cell size
  if(is.null(res.raster)){
    if(is.null(cells)){
      print("=========== Determining raster with largest cell size ===========")
      max.raster <- raster.table$RasterName[max(raster.table$ResolutionX) == raster.table$ResolutionX]
      max.raster <- raster(paste(as.character(max.raster)[1]))
    }else{
      print("=========== Resampling raster to defined cell size ===========")
      max.raster <- raster(list.rasters[1])
      max.raster <- aggregate(max.raster, fact = cells, fun = mean, expand = FALSE, na.rm = TRUE)
    }
  }else{
    print("=========== Determining raster with largest cell size ===========")
    max.raster <- raster(list.rasters[1])
    for(i in 2:length(list.rasters)){
      max.raster.res <- max(res(max.raster))
      max.raster <- if(max(res(raster(list.rasters[i]))) > max.raster.res)
      {raster(list.rasters[i])} else {max.raster}
    }
  }

  #Check to make sure all projections are the same
  print("=========== Determining raster projections are the same ===========")
  for(i in 1:length(list.rasters)){
    r <- raster(list.rasters[i])
    if(is.na(projection(r))){
      projection(r) <- projection(max.raster)
      warning(paste("Projection of raster ",r@data@names," is NA. Reprojecting to sample raster.", sep = ""))
    }else if(projection(r) != projection(max.raster)){
      r.project <- projectRaster(from = r, to = max.raster, alignOnly = TRUE)
      warning(paste("Projection of raster ",r@data@names," differs. Reprojecting to sample raster.", sep = ""))
    }
  }

  #Read in mask to clip rasters
  if(is.null(mask) == FALSE){
    print("=========== Cropping rasters to the extent of mask ===========")
    if(is.na(projection(max.raster) ==  projection(mask))){
      projection(mask) <- projection(max.raster)
      max.raster <- crop(max.raster, mask, snap = "out")
      warning("Projection of mask shapefile differs from raster. Reprojecting to sample raster.")
    }
    else if(projection(max.raster) !=  projection(mask)){
      projection(mask) <- projection(max.raster)
      max.raster <- crop(max.raster, mask, snap = "out")
      warning("Projection of mask shapefile differs from raster. Reprojecting to sample raster.")
    }else{
      max.raster <- crop(max.raster, mask, snap = "out")
    }
  }

  #Resample rasters
  print("=========== Resampling rasters to sample raster resolution ===========")
  for(i in 1:length(list.rasters)){
    r <- raster(list.rasters[i])

    #Rescale rasters
    if(rescale == TRUE){
      r <- scale(r, center = TRUE, scale = TRUE)
    }

    #Extract name of raster
    name <- r@data@names
    print(paste("=========== Outputting ",name," ===========", sep = ""))

    if(is.null(mask) == FALSE){
      #Resample and output as .asc file
      r.crop <- crop(r, extent(mask), snap = "out")
      r.crop <- mask(r.crop, mask = mask)
      r <- resample(r.crop, max.raster, overwrite=TRUE, filename = paste(output.path, name, "_resample.asc", sep=""))
    } else{r <- resample(r, max.raster, overwrite=TRUE, filename = paste(output.path, name, "_resample.asc", sep=""))
    }

    #Re-read back in for correct format for GARP
    # r <- raster(paste(output.path, name, "_resample.asc", sep=""))
    # write.asciigrid(r, paste(output.path, name, "_resample.asc", sep=""))
  }
}
