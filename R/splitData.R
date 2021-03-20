#' Splits presence locations into training and testing data sets for use by GARP
#'
#' \code{splitData} Uses a user-defined precentage to split presence locations into training and testing datasets. Both are output as .csv files (need to be converted to .xls for use in the desktop version of GARP) and a shapefile object (for testing using \code{\link{aucGARP}}).
#'
#' @param points a data.frame object containing latitude and longitude values of presence locations
#' @param p a numeric value specifying the precentage of locations to keep for training
#' @param type a character value specifying which points should be returned: training ("train"), testing ("test"), or both ("all")
#' @param iterations a numeric value specifying how many resampling events should be run; defaults to 1
#' @param output logical. Should the datasets be saved to file?
#' @param output.path a file path specifying the file name of the objects to be saved. "_train" or "_test" will be appended to file name (e.g. "C:/Users/")
#' @param output.name the file name of the objects to be saved. "_train" or "_test" will be appended to file name (e.g. "Homo_sapiens_test.shp")
#'
#' @return A list object containing data.frames for training locations (when \code{type = "train"}), testing locations (when \code{type = "test"}), or both (when \code{type = "all"}). When output = TRUE, outputs both a .csv and shapefile to specified file path.
#'
#' @details Use the \code{iterations} element when performing multiple experiments in GARP. The function will return a single .csv datasheet with all training data combined for use by GARP, with
#' the "Species" vector containing the species name from the \code{points} object with the iteration number attached (e.g. "Homo_sapiens1_train", "Homo_sapiens2_train", "Homo_sapiens3_train", etc.).
#'
#'
#' @examples
#' \dontrun{
#'   hs <- data.frame("Latitude" = c(-89, 72, 63, 42, 54), "Longitude" = c(-12, 13, 24, 26, 87), "Species" = rep("Homo sapiens", 5))
#'   splitData(points = hs, p = 0.7, type = "all", output.path = "C/GARP", output.name = "Homo_sapiens")
#' }
#' @import raster
#' @import sp
#'
#' @export

splitData <- function(points, p, type = c("train", "test", "all"), iterations = 1, output = TRUE, output.path, output.name){

  if("Latitude" %in% names(points) == FALSE) {stop("Coordinates must be labeled 'Latitude' and 'Longitude'.")}
  if("Longitude" %in% names(points) == FALSE) {stop("Coordinates must be labeled 'Latitude' and 'Longitude'.")}
  if("Species" %in% names(points) == FALSE) {stop("Must contain vector of species names labeled 'Species'.")}
  species <- levels(points$Species)

  df.train <- data.frame()
  for(n in 1:paste(iterations)){
    train.pts.format.all <- data.frame()
    test.pts.format.all  <- data.frame()

    for(i in species){

      #Subset for species
      points$Species <- as.factor(points$Species)
      points.sp <- subset(points, points$Species == i)

      #Add unique id
      points.sp$id <- seq(1, length(points.sp$Species))

      #Pull a random sample dataset
      sample <- sample(points.sp$id,ceiling(length(points.sp) * p))

      #Subset points into training points
      train.pts <- points.sp[(points.sp$id %in% sample), ]

      #Subset other points as testing points
      test.pts <- points.sp[!(points.sp$id %in% sample), ]

      #Make sure in correct format for GARP
      train.pts.format <- data.frame(paste(train.pts$Species,n,sep=""), train.pts$Longitude, train.pts$Latitude)
      test.pts.format  <- data.frame(paste(test.pts$Species,n,sep=""), test.pts$Longitude, test.pts$Latitude)
      colnames(train.pts.format) <- c("Species", "Longitude", "Latitude")
      colnames(test.pts.format)  <- c("Species", "Longitude", "Latitude")

      train.pts.format.all <- rbind(train.pts.format.all, train.pts.format)
      test.pts.format.all  <- rbind(test.pts.format.all, test.pts.format)
    }

    #Convert to spatial points dataframe
    test.pts.all  <- SpatialPointsDataFrame(as.matrix(data.frame(test.pts.format.all$Longitude,test.pts.format.all$Latitude)), data = test.pts.format.all)

    #Output testing files
    if(type == "test"){
      if(output == TRUE){
        #Write csv file
        write.csv(test.pts.format.all, paste(output.path,output.name,n,"_test.csv", sep = ""), row.names=FALSE)
        #Write shapefile
        shapefile(test.pts.all, filename = paste(output.path, output.name,n,"_test.shp", sep=""), overwrite = TRUE)
      } else if(output == FALSE){return(test.pts.format.all)}
    } else if(type == "all"){
      if(output == TRUE){
        #Write csv file
        write.csv(test.pts.format.all, paste(output.path,output.name,n,"_test.csv", sep = ""), row.names=FALSE)
        #Write shapefile
        shapefile(test.pts.all, filename = paste(output.path, output.name,n,"_test.shp", sep=""), overwrite = TRUE)
      } else if(output == FALSE){return(test.pts.format.all)}
    }

    #Append training files
    df.train <- rbind(df.train, train.pts.format.all)
  }

  #Convert to spatial points datafame
  df.train.sp <- SpatialPointsDataFrame(as.matrix(data.frame(df.train$Longitude,df.train$Latitude)), data = df.train)

  #Output files based on type call
  if(type == "train"){
    if(output == TRUE){
      #Write csv file
      write.csv(df.train.sp, paste(output.path,output.name,"_train.csv", sep = ""), row.names=FALSE)
      #Write shapefile
      shapefile(df.train.sp, filename = paste(output.path, output.name,n,"_train.shp", sep=""), overwrite = TRUE)
    } else if(output == FALSE){return(df.train.sp)}
  } else if(type == "all"){
    if(output == TRUE){
      #Write csv file
      write.csv(df.train.sp, paste(output.path,output.name,"_train.csv", sep = ""), row.names=FALSE)
      #Write shapefile
      shapefile(df.train.sp, filename = paste(output.path, output.name,n,"_train.shp", sep=""), overwrite = TRUE)
    } else if(output == FALSE){return(df.train.sp)}
  }
}
