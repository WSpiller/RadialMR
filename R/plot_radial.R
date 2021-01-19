#' plot_radial
#'
#' A function for producing radial IVW and MR-Egger plots either individually or simultaneously. The function allows for a variety of aesthetic and scaling options, utilising the output from the \code{IVW_radial} and \code{egger_radial} functions.
#'
#' @param r_object An object of class \code{"IVW"} or \code{"egger"}. For visualising both estimates simultaneously, both objects should be included as a vector c(\code{A},\code{B}), where \code{A} and \code{B} denote the \code{"IVW"} and \code{"egger"} objects respectively.
#' @param radial_scale Indicates whether to produce a plot including a full radial scale (\code{TRUE}), or a scatterplot showing only the effect estimates (\code{FALSE}).
#' @param show_outliers Indicates whether display only the set of variants identified as outliers (\code{TRUE}) or the complete set of variants (\code{FALSE}). Note that when (\code{show_outliers=TRUE}), non-outlying variants further from the origin than the furthest outlier will cause an error message that one or more points have been omitted. These are non-outlying variants beyond the scale. If no outliers are present, a plot will be produced using the full set of variants, with an accompanying message indicating the absence of outliers.
#' @param scale_match Indicates whether x and y axes should have the same range(\code{TRUE}), or different ranges (\code{FALSE}) This improves the interpretation of the radial scale, and is set to \code{FALSE} when the radial scale is omitted from the plot.
#' @return A ggplot object containing a radial plot of either the IVW, MR-Egger, or both estimates simultaneously.
#' @author Wes Spiller; Jack Bowden; Tom Palmer.
#' @references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.
#' @importFrom stats quantile median
#' @export
#' @examples
#' ldl.dat <- data_radial[data_radial[,10]<5*10^-8,]
#' ldl.fdat <- format_radial(ldl.dat[,6], ldl.dat[,9],
#'                           ldl.dat[,15], ldl.dat[,21], ldl.dat[,1])
#' ivw.object <- ivw_radial(ldl.fdat, 0.05, 1, 0.0001, TRUE)
#' plot_radial(ivw.object)

plot_radial<-function(r_object,radial_scale,show_outliers,scale_match){


  if(missing(radial_scale)) {
    radial_scale <- TRUE
  }

  if(missing(show_outliers)) {
    show_outliers <- FALSE
  }

  if(missing(scale_match)) {
    scale_match <- TRUE
  }

  # for R CMD check
  Wj <- BetaWj <- Outliers <- out_types <- layout <- NULL

  #Function for producing scale and confidence interval circle indicators

  circle.radial <- function(center = c(0,0), radius, num.points, START,END){

    R = radius

    tt <- seq(START,END,length.out = num.points)

    #Generates x-axis values for circle
    xx <- center[1] + R * cos(tt)

    #Generates y-axis values for circle
    yy <- center[2] + R * sin(tt)

    return(data.frame(x = xx, y = yy))

  }

  if(length(r_object)>=6 && length(r_object)<=13){

    if(class(r_object)=="IVW"){

      r_object$coef<-c(r_object$coef[2,])
      r_object$coef<-as.numeric(r_object$coef)

      #Define radius for scale and IVW circles
      maximal<-atan(max(abs(r_object$data[,3]))/max(r_object$data[,2]))
      R.All<-max(abs(r_object$data[,3]))/sin(maximal)
      R.IVW<-R.All+(R.All/12)

      #Define scaling parameter
      Label.Scaling<-R.IVW/24

      Scalemin<-min(r_object$data[,3]/r_object$data[,2],r_object$confint[1])

      Scalemax<-max(r_object$data[,3]/r_object$data[,2],r_object$confint[2])

      #Plot circle for overall scale
      circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
      cxAll<-circledat_ALLEST$x
      cyAll<-circledat_ALLEST$y

      #Plot circle for IVW estimate
      circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
      cxIVW<-circledat_IVWEST$x
      cyIVW<-circledat_IVWEST$y

      #Define minimum and maximum values of y axis for scaling
      Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
      Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)

      Y_Range<-Scalemax-Scalemin

      #Function for rounding to nearest 0.05 increment
      mround <- function(x){
        0.5*round(x/0.5)
      }

      #If only outliers are to be shown
      if(show_outliers==TRUE){

        #If the full scale is to be shown
        if(radial_scale==TRUE){

          #If axes scales should be fixed
          if(scale_match==TRUE){

            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){

              #Produce plot showing full scale and all variants
              B<-ggplot2::ggplot(r_object$data,ggplot2::aes(x=Wj,y=BetaWj))+ggplot2::labs(title="IVW Radial")+ ggplot2::geom_point(ggplot2::aes(colour="Variant"))+
                ggplot2::geom_path(ggplot2::aes(x=cxAll,y=cyAll))+ggplot2::geom_path(ggplot2::aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                ggplot2::geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),
                                 axis.line = ggplot2::element_line(colour = "black"))+ggplot2::ylab(expression(hat(beta)[j]~sqrt(W[j])))+ggplot2::xlab(expression(sqrt(W[j])))+
                ggplot2::geom_segment(ggplot2::aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                ggplot2::scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+ggplot2::scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                ggplot2::coord_fixed()+ggplot2::scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))


              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)

                for (i in 1:length(Y_Scale)){
                  B<- B + ggplot2::geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)

                }
              }

              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){

                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)

                for (i in 1:length(Y_Scale)){
                  B<- B + ggplot2::geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)

                }
              }
              for(i in 1:length(r_object$data[,3])){
                Theta<-atan(r_object$data[,3]/r_object$data[,2])

                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All

                B<- B + ggplot2::geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }

              print("No significant Outliers")

            }

            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){

              #Redefine scale and IVW circle radii
              maximal<-atan(max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$BetaWj))/max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$Wj)))
              R.All<-max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$BetaWj))/sin(maximal)
              R.IVW<-R.All+(R.All/12)

              #Plot circle for overall scale
              circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y

              #Plot circle for IVW estimate
              circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y

              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)

              Y_Range<-Scalemax-Scalemin

              #Generate plot showing only outliers using the full scale
              B<-ggplot2::ggplot(r_object$data,ggplot2::aes(x=Wj,y=BetaWj))+ggplot2::labs(title="IVW Radial") + ggplot2::geom_point(ggplot2::aes(colour = Outliers))+ggplot2::geom_path(ggplot2::aes(x=cxAll,y=cyAll))+ggplot2::geom_path(ggplot2::aes(x=cxIVW,y=cyIVW,colour="IVW"))+

                ggplot2::geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[1],digits=3),size=3)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[2],digits=3),size=3)+
                ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"))+ggplot2::ylab(expression(hat(beta)[j]~sqrt(W[j])))+ggplot2::xlab(expression(sqrt(W[j])))+
                ggplot2::geom_segment(ggplot2::aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                ggplot2::scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+ggplot2::scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                ggplot2::coord_fixed()+ggplot2::theme(legend.title=ggplot2::element_blank())+ggplot2::scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))

              #Function for rounding to nearest 0.05 increment
              mround <- function(x){
                0.5*round(x/0.5)
              }


              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)

                for (i in 1:length(Y_Scale)){
                  B<- B + ggplot2::geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)

                }
              }

              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){

                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)

                for (i in 1:length(Y_Scale)){
                  B<- B + ggplot2::geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)

                }
              }


              #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
              Theta<-rep(0,length(r_object$data[r_object$data$Outliers == "Outlier", ]$Outliers))

              for(i in 1:length(r_object$data[r_object$data$Outliers == "Outlier", ][,3])){
                Theta<-atan(r_object$data[r_object$data$Outliers == "Outlier", ][,3]/r_object$data[r_object$data$Outliers == "Outlier", ][,2])
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + ggplot2::geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                B<- B + ggplot2::geom_segment(x = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], xend = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i],
                                     y = r_object$data[r_object$data$Outliers == "Outlier", ][,3][i], yend =r_object$data[r_object$data$Outliers == "Outlier", ][,2][i]*
                                       r_object$coef[1],linetype="solid",colour="#56B4E9")
              }

            }

          }

          #If axes scales should be fixed
          if(scale_match==FALSE){
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){

              #Produce plot showing full scale and all variants
              B<-ggplot2::ggplot(r_object$data,ggplot2::aes(x=Wj,y=BetaWj))+ggplot2::labs(title="IVW Radial")+ ggplot2::geom_point(ggplot2::aes(colour="Variant"))+
                ggplot2::geom_path(ggplot2::aes(x=cxAll,y=cyAll))+ggplot2::geom_path(ggplot2::aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                ggplot2::geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),
                                 axis.line = ggplot2::element_line(colour = "black"))+ggplot2::ylab(expression(hat(beta)[j]~sqrt(W[j])))+ggplot2::xlab(expression(sqrt(W[j])))+
                ggplot2::geom_segment(ggplot2::aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                ggplot2::scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+ggplot2::scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                ggplot2::scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))


              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)

                for (i in 1:length(Y_Scale)){
                  B<- B + ggplot2::geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)

                }
              }

              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){

                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)

                for (i in 1:length(Y_Scale)){
                  B<- B + ggplot2::geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)

                }
              }
              for(i in 1:length(r_object$data[,3])){
                Theta<-atan(r_object$data[,3]/r_object$data[,2])

                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All

                B<- B + ggplot2::geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }

              print("No significant Outliers")

            }

            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){

              #Redefine scale and IVW circle radii
              maximal<-atan(max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$BetaWj))/max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$Wj)))
              R.All<-max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$BetaWj))/sin(maximal)
              R.IVW<-R.All+(R.All/12)

              #Plot circle for overall scale
              circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y

              #Plot circle for IVW estimate
              circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y

              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)

              Y_Range<-Scalemax-Scalemin

              #Generate plot showing only outliers using the full scale
              B<-ggplot2::ggplot(r_object$data,ggplot2::aes(x=Wj,y=BetaWj))+ggplot2::labs(title="IVW Radial") + ggplot2::geom_point(ggplot2::aes(colour = Outliers))+ggplot2::geom_path(ggplot2::aes(x=cxAll,y=cyAll))+ggplot2::geom_path(ggplot2::aes(x=cxIVW,y=cyIVW,colour="IVW"))+

                ggplot2::geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[1],digits=3),size=3)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[2],digits=3),size=3)+
                ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"))+ggplot2::ylab(expression(hat(beta)[j]~sqrt(W[j])))+ggplot2::xlab(expression(sqrt(W[j])))+
                ggplot2::geom_segment(ggplot2::aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                ggplot2::scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+ggplot2::scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                ggplot2::theme(legend.title=ggplot2::element_blank())+ggplot2::scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))




              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)

                for (i in 1:length(Y_Scale)){
                  B<- B + ggplot2::geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)

                }
              }

              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){

                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)

                for (i in 1:length(Y_Scale)){
                  B<- B + ggplot2::geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)

                }
              }
              #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
              Theta<-rep(0,length(r_object$data[r_object$data$Outliers == "Outlier", ]$Outliers))

              for(i in 1:length(r_object$data[r_object$data$Outliers == "Outlier", ][,3])){
                Theta<-atan(r_object$data[r_object$data$Outliers == "Outlier", ][,3]/r_object$data[r_object$data$Outliers == "Outlier", ][,2])
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + ggplot2::geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                B<- B + ggplot2::geom_segment(x = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], xend = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i],
                                     y = r_object$data[r_object$data$Outliers == "Outlier", ][,3][i], yend =r_object$data[r_object$data$Outliers == "Outlier", ][,2][i]*
                                       r_object$coef[1],linetype="solid",colour="#56B4E9")

              }
            }
          }

        }

        #If the full scale is to be shown
        if(radial_scale==FALSE){

          #If axes scales should be fixed
          if(scale_match==TRUE){

            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){

              #Produce plot showing full scale and all variants
              B<-ggplot2::ggplot(r_object$data,ggplot2::aes(x=Wj,y=BetaWj))+ggplot2::labs(title="IVW Radial")+ ggplot2::geom_point(ggplot2::aes(colour="Variant"))+
                ggplot2::geom_path(ggplot2::aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                ggplot2::geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),
                                 axis.line = ggplot2::element_line(colour = "black"))+ggplot2::ylab(expression(hat(beta)[j]~sqrt(W[j])))+ggplot2::xlab(expression(sqrt(W[j])))+
                ggplot2::geom_segment(ggplot2::aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                ggplot2::scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+ggplot2::scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                ggplot2::coord_fixed()+ggplot2::scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))

              print("No significant Outliers")

            }

            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){

              #Redefine scale and IVW circle radii
              maximal<-atan(max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$BetaWj))/max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$Wj)))
              R.All<-max(abs(r_object$data[r_object$data$Outliers == "Outlier", ]$BetaWj))/sin(maximal)
              R.IVW<-R.All+(R.All/12)

              #Plot circle for overall scale
              circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y

              #Plot circle for IVW estimate
              circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y

              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)

              Y_Range<-Scalemax-Scalemin



              #Generate plot showing only outliers using the full scale
              B<-ggplot2::ggplot(r_object$data,ggplot2::aes(x=Wj,y=BetaWj))+ggplot2::labs(title="IVW Radial") + ggplot2::geom_point(ggplot2::aes(colour = Outliers))+ggplot2::geom_path(ggplot2::aes(x=cxIVW,y=cyIVW,colour="IVW"))+

                ggplot2::geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[1],digits=3),size=3)+
                ggplot2::geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[2],digits=3),size=3)+
                ggplot2::theme_bw()+ggplot2::theme(panel.border = ggplot2::element_blank(),panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(),axis.line = ggplot2::element_line(colour = "black"))+ggplot2::ylab(expression(hat(beta)[j]~sqrt(W[j])))+ggplot2::xlab(expression(sqrt(W[j])))+
                ggplot2::geom_segment(ggplot2::aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                ggplot2::scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+ggplot2::scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                ggplot2::coord_fixed()+ggplot2::theme(legend.title=ggplot2::element_blank())+ggplot2::scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))

              #Draw lines indicating individual Q contributions and give text displaying each value
              for(i in 1:length(r_object$data[r_object$data$Outliers == "Outlier", ][,3])){




                B<- B + ggplot2::geom_segment(x = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], xend = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], y = r_object$data[r_object$data$Outliers == "Outlier", ][,3][i], yend =r_object$data[r_object$data$Outliers == "Outlier", ][,2][i]*r_object$coef[1],linetype="solid",colour="#56B4E9")
