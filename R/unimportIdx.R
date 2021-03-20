#' Calculates the unimportance index of each environmental covariate for variable selection
#'
#' \code{unimportIdx} Calculates the unimportance index of each environmental covariate for variable selection
#'
#' @param table1 a data.frame recording the prevalence of each environmental covariate; can be produced by \code{\link{prevalence}}
#' @param table2 a data.frame with the rescaled median range for each environmental variable; can be produced by \code{\link{rescale}}
#'
#' @return A table recording the normalized unimportance index for each environmental covariate
#'
#'
#'
#' @export

unimportIdx<-function(table1,table2){
  table3<-table1
  table3$Names=as.character(table3$Names)
  table2$RuleNames=as.character(table2$RuleNames)
  table3<-table3[!table3$Names=='mask',]
  table2<-table2[!table2$RuleNames=='mask',]
  for (i in 1:length(table3$Names)){
    for(j in 1:length(table2$RuleNames)){
      if (table3$Names[i]==table2$RuleNames[j]){
        table3$MedRange[i]<-table2$Range[j]
      }
    }
  }

  table3$Unimportance<-(1-table3$Prevalence) * table3$MedRange

  for (k in 1:length(table3$MedRange)){
    table3$RescaledUnimIDx[k]<-(table3$Unimportance[k]-min(table3$Unimportance))/(max(table3$Unimportance)-min(table3$Unimportance))
  }
  return(table3)
}
