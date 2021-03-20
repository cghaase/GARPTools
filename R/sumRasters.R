#' Sums rasters of the best subset models output from GARP
#'
#' \code{sumRasters} GARP outputs binary presence/absence grids for the sample area and selects the best models based on user-defined qualifications. This function
#' sums the best models into a single raster that represents the number of models that agree on a predicted presence location.
#'
#' @param file.path path to the location of output rasters to sum (usually found in "GARPruns" folder output from GARP)
#' @param output logical. Should the summated raster be saved to file?
#' @param output.path a file path specifying where the output raster should be saved, including file name and type (e.g. "C:/Users/sumraster1.asc"). Can be of type ".grd", ".asc", ".tif", etc. See \href{https://cran.r-project.org/web/packages/raster/raster.pdf}{writeFormats} for supported types
#' @param format the format type of the raster to be saved. Can be of type "raster", "ascii", "GTiff", etc. See \href{https://cran.r-project.org/web/packages/raster/raster.pdf}{writeFormats} for supported types
#' @param over.write logical. Should the raster be overwritten?
#'
#' @return A single raster with cell values that represent the number of models that agree on a predicted presence location.
#'
#' @import raster
#' @import sp
#'
#' @export

sumRasters <- function(file.path, output = FALSE, output.path = NULL, format = NULL, over.write = FALSE){

  #List all ASCII files in folder to sum
  list.files <- list.files(path = paste(file.path), pattern = ".asc$",  full.names = TRUE)
  for(i in 1:length(list.files)){
    if("info" %in% strsplit(list.files[i], split= "/")[[1]][length(strsplit(list.files[i], split= "/")[[1]])]){list.files <- list.files[-i]}
  }

  #Start with a 0 raster
  sum.raster = 0

  #Iteratively add the rasters from the file
  for(i in 1:length(list.files)){
    task <- raster(list.files[i])
    sum.raster <- (task + sum.raster)
  }

  #Write out raster
  if(output==TRUE){
    if(is.null(format) == TRUE){
      stop('Output format not correctly specified.')
    }
    else{writeRaster(sum.raster, filename = paste(output.path), format = paste(format), overwrite = over.write)
      return(sum.raster)
    }
  }
  else{return(sum.raster)
  }
}
