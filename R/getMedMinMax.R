#' Pull rulesets created by GARP
#'
#' \code{getMedMinMax} Returns the minimum and maximum median values for variables in dominant rules from a given rule-set
#'
#' @param pathA a vector object that contains the full file names for the extracted rule .csv files created by \code{\link{extractRules}}
#' @param pathB a file path of the directory containing environmental layers used in GARP (must be in ESRI ASCII raster format)
#'
#' @return Returns a data.frame that contains the range, and median of the minimum and maximum values for each environment variable used in GARP.
#'
#' @seealso \code{\link{extractRules}}
#'
#' @examples
#' \dontrun{
#' MedMinMax <- getMedMinMax(pathA = "C:/GARP/RangeLogit.csv", pathB = "C:/GARP/resampled/")
#' }
#'
#' @import sp
#' @import raster
#'
#' @export

getMedMinMax<-function(pathA,pathB){
  RuleDir<-pathA
  envdir<-pathB
  DataAll<-data.frame()
  for (i in seq(1:length(RuleDir))){
    DataAll<-rbind(DataAll,assign(paste("Data", i, sep = "."), read.csv(RuleDir[i])))
  }
  envdir<-envdir
  envnames<-list.files(envdir,pattern=".asc")
  Names<-c()
  MedMax<-c()
  MedMin<-c()
  Max<-c()
  Min<-c()
  MedRange<-c()
  Range<-c()
  for (i in seq(1:length(envnames))){
    Names=append(Names,strsplit(envnames[i],'.asc')[[1]])
    a<-median(assign(paste('R.',strsplit(envnames[i],'.asc')[[1]]),subset(DataAll,DataAll$RuleName==strsplit(envnames[i],'.asc')[[1]]))$Maximum)
    b<-median(assign(paste('R.',strsplit(envnames[i],'.asc')[[1]]),subset(DataAll,DataAll$RuleName==strsplit(envnames[i],'.asc')[[1]]))$Minimum)
    c<-as.numeric(max(assign(paste('R.',strsplit(envnames[i],'.asc')[[1]]),subset(DataAll,DataAll$RuleName==strsplit(envnames[i],'.asc')[[1]]))$Maximum))
    d<-as.numeric(min(assign(paste('R.',strsplit(envnames[i],'.asc')[[1]]),subset(DataAll,DataAll$RuleName==strsplit(envnames[i],'.asc')[[1]]))$Minimum))
    MedMax<-append(MedMax,a)
    MedMin<-append(MedMin,b)
    Max<-append(Max,c)
    Min<-append(Min,d)
    MedRange<-append(MedRange,a-b)
    Range<-append(Range,c-d)
  }
  Names<-as.character(Names)
  MedMaxMin<-data.frame(Names,MedMin,MedMax,MedRange,Min, Max, Range)
  return(MedMaxMin)
}
