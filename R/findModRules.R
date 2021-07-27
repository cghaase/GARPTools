#' Pull rulesets created by GARP
#'
#' \code{findModRules} Gets the dominant rules produced by GARP which cover more than 90% of the landscape for the best subset models
#'
#' @param pathA file path to the directory for ruleset grids produced by GARP
#' @param pathB a vector object that contains the directories of the best subsets for each GARP experiment
#'
#' @return Returns a data.frame that includes the model numbers (ruleset number) and dominant rule numbers
#'
#' @seealso \code{\link{extractRules}}
#'
#' @examples
#' \dontrun{
#' modelRule <- findModRules(pathA = "C:/GARP/runs/", pathB = "C:/GARP/BestSubsets/")
#' }
#'
#' @import sp
#' @import raster
#'
#' @export

findModRules <- function(pathA,pathB){
  rsetPridir<-pathA
  bsdirlist<-pathB
  BestRules<-c()
  bestmodel<-c()
  rsetlist<-list.files(path=rsetPridir, pattern = 'grid')
  #find the rset for best subsets, do the zonal stats, and find the rules covered 90%
  for (i in 1:length(bsdirlist)){
    bslist <- list.files(path = bsdirlist[i],pattern='task', all.files = TRUE, full.names = FALSE)
    for (j in 1: length(bslist)){
      if (grepl('\\.[a-z]',bslist[j])==FALSE){
        best <- regmatches(bslist[j], regexpr("\\d+", bslist[j]))
        best <- as.numeric(best)
        for (l in 1:length(rsetlist)){
          rsetdir <- paste(rsetPridir,rsetlist[l],sep = '')
          rsetFile <- list.files(path = rsetdir, pattern = 'rset_\\d+_0')
          for (k in 1:length(rsetFile)){
            if (grepl('\\.[a-z]',rsetFile[k])==FALSE){
              rset <- regmatches(rsetFile[k], regexpr('\\d+_',rsetFile[k]))
              rset <- regmatches(rset, regexpr('\\d+',rset))
              rset <- as.numeric(rset)
              if (rset==best){
                rsetRa <- paste(rsetdir,'/',sep = '')
                rsetRa <- paste(rsetRa, rsetFile[k], sep = '')
                rsetRa <- raster(rsetRa)
                dbf <- zonal(rsetRa,rsetRa,'mean')
                dbf<-dbf[order(-dbf$length),]
                dbf<-data.frame(dbf)
                time=0
                total=0
                for (m in 1:length(dbf$length)){
                  total=total+dbf$length[m]
                  if (total<=0.9*sum(dbf$length)){
                    if (dbf$zones[m] != 0){
                      #print (dbf$zones[m])
                      bestmodel<-append(bestmodel,best)
                      BestRules<-append(BestRules,as.numeric(paste(dbf$zones[m])))
                    }
                  }else {
                    time=time+1
                    if (time == 1){
                      if (dbf$zones[m] != 0){
                        bestmodel<-append(bestmodel,best)
                        BestRules<-append(BestRules,as.numeric(paste(dbf$zones[m])))
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  modelRule<-cbind(bestmodel,BestRules)
  modelRule<-data.frame(modelRule)
  return(modelRule)
}
