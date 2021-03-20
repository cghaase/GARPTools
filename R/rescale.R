#' Rescales environmental variables from 0 to 1
#'
#' \code{rescale} Rescales the the minimum and maximum median values for each environment variable used in dominant rules across a GARP Best Subsets Procedure to be between 0 to 1. This will allow the user to evalaute the range of values across variables regardless of units of measure
#'
#' @param file a data.frame of the original median of the minimum and maximum values for each environmental variable produced by \code{\link{getMedMinMax}}
#'
#' @return Returns a data.frame containing the rescaled range, and the median of the minimum and maximum values for each environmental variable.
#'
#' @seealso \code{\link{getMedMinMax}}
#'
#' @import sp
#' @import raster
#'
#' @export

rescale<-function(file){
  MedMinMax<-file
  if (colnames(MedMinMax)[1]!='X'){
    MedMaxMin<-MedMinMax[,2:3]
    MedMinimum<-c()
    MedMaximum<-c()
    Range<-c()
    RuleNames<-c()
    for (j in seq(1:length(MedMinMax$MedMin))){
      a<-(MedMinMax$MedMin[j]-MedMinMax$Min[j])/MedMinMax$Range[j]
      b<-(MedMinMax$MedMax[j]-MedMinMax$Min[j])/MedMinMax$Range[j]
      MedMinimum[j]<-a
      MedMaximum[j]<-b
      Range[j]<-b-a
      RuleNames[j]<-as.character(MedMinMax$Names[j])
    }
  } else {
    MedMaxMin<-MedMinMax[,3:4]
    MedMinimum<-c()
    MedMaximum<-c()
    Range<-c()
    RuleNames<-c()
    for (j in seq(1:length(MedMaxMin$MedMin))){
      a<-(MedMinMax$MedMin[j]-MedMinMax$Min[j])/MedMinMax$Range[j]
      b<-(MedMinMax$MedMax[j]-MedMinMax$Min[j])/MedMinMax$Range[j]
      MedMinimum[j]<-a
      MedMaximum[j]<-b
      Range[j]<-b-a
      RuleNames[j]<-as.character(MedMinMax$Names[j])
    }
  }
  MedRescale<-data.frame(RuleNames,MedMinimum,MedMaximum,Range)
  return(MedRescale)
}
