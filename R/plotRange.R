#' Plots the rescaled median ranges for environment variables
#'
#' \code{plotRange} Plots the rescaled median ranges for environment variables used across dominant rules from a GARP experiment
#'
#' @param file a data.frame of the rescaled range, and the median of the minimum and maximum values for each environmental variables used in GARP. Produced by \code{\link{rescale}}
#' @param colour a character string designating the color of the bars of the plot
#' @param n numeric; specifies the width of the bars
#' @param xlabel label for the x-axis
#' @param ylabel label for the y-axis
#'
#' @return Returns a plot of the rescaled median ranges for environment variables used in DesktopGARP.
#'
#' @import ggplot2
#'
#' @export


plotRange<-function(file, colour = "purple",n = 1,xlabel = "Environmental Variables",ylabel = "Median Range"){
  rescaledcsv<-file
  rescaledcsv$Mean<-rescaledcsv$Range/2
  MedRange<-ggplot(rescaledcsv, aes(x=RuleNames, y=Mean))+coord_flip()+geom_errorbar(aes(ymin=rescaledcsv$MedMinimum,
                                                                                         ymax=rescaledcsv$MedMaximum), colour=colour,width=0, size=n)+xlab(xlabel)+ylab(ylabel)
  return(MedRange)
}
