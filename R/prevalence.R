#' Calculate the prevalence of each environmental covariate in the rulesets
#'
#' \code{prevalence} Calculate the prevalence of each environmental covariate in the ruleset file
#'
#' @param pathA a vector object that contains the full file names for the extracted rule .csv files created by \code{\link{extractRules}}
#' @param pathB the primary directory for ruleset grids
#' @param num numeric; the total number of dominant presence rules; can be generated from \code{\link{totPresRules}}
#'
#' @return A table recording the prevalence of each environmental covariate.
#'
#' @details Ruleset grids called by the \code{pathB} path must be in ESRI ASCII format.
#'
#' @import sp
#' @import raster
#'
#' @export

prevalence<-function(pathA,pathB,num){
  RuleDir<-pathA
  envdir<-pathB
  TotPresNum<-num
  DataAll<-data.frame()
  for (i in seq(1:length(RuleDir))){
    DataAll<-rbind(DataAll,assign(paste("Data", i, sep = "."), read.csv(RuleDir[i])))
  }
  envdir<-envdir
  envnames<-list.files(envdir,pattern=".asc")
  Names<-c()
  PresNumEn<-c()
  TotPres<-c()
  Prevalence<-c()
  for (i in seq(1:length(envnames))){
    Names=append(Names,strsplit(envnames[i],'.asc')[[1]])
    a<-nrow(assign(paste('R.',strsplit(envnames[i],'.asc')[[1]],sep=''),subset(DataAll,DataAll$RuleName==strsplit(envnames[i],'.asc')[[1]])))
    b<-a/TotPresNum
    PresNumEn<-append(PresNumEn,a)
    Prevalence<-append(Prevalence,b)
    TotPres<-append(TotPres,TotPresNum)
  }
  Names<-as.character(Names)
  Prevalence<-data.frame(Names,PresNumEn,TotPres,Prevalence)
  return(Prevalence)
}
