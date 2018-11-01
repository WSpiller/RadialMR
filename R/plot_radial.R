#' plot_radial
#'
#' A function for producing radial IVW and MR-Egger plots either individually or simultaneously. The function allows for a variety of aesthetic and scaling options, utilising the output from the \code{IVW_radial} and \code{egger_radial} functions.
#'
#' @param r_object An object of class \code{"IVW"} or \code{"egger"}. For visualising both estimates simultaneously, both objects should be included as a vector c(\code{A},\code{B}), where \code{A} and \code{B} denote the \code{"IVW"} and \code{"egger"} objects respectively.
#' @param radial_scale Indicates whether to produce a plot including a full radial scale (\code{TRUE}), or a scatterplot showing only the effect estimates (\code{FALSE}).
#' @param show_outliers Indicates whether display only the set of variants identified as outliers (\code{TRUE}) or the complete set of variants (\code{FALSE}). Note that when (\code{show_outliers=TRUE}), non-outlying variants further from the origin than the furthest outlier will cause an error message that one or more points have been omitted. These are non-outlying variants beyond the scale. If no outliers are present, a plot will be produced using the full set of variants, with an accompanying message indicating the absence of outliers.
#' @param scale_match Indicates whether x and y axes should have the same range(\code{TRUE}), or different ranges (\code{FALSE}) This improves the interpretation of the radial scale, and is set to \code{FALSE} when the radial scale is omitted from the plot.
#' @return A ggplot object containing a radial plot of either the IVW, MR-Egger, or both estimates simultaneously.
#'@author Wes Spiller; Jack Bowden.
#'@references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.
#'@export
#'@examples
#'
#'plot_radial(r_object,TRUE,TRUE,TRUE) 
#'
#'plot_radial(c(r_object,r_object),TRUE,TRUE,FALSE)

plot_radial<-function(r_object,radial_scale,show_outliers,scale_match){
  
  if(missing(radial_scale)) {
    radial_scale<-T
  }
  
  if(missing(show_outliers)) {
    show_outliers<-F
  }
  
  if(missing(scale_match)) {
    scale_match<-T
  }
  
  #Load ggplot2 library
  library(ggplot2)
  
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
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial")+ geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              for(i in 1:length(r_object$data[,3])){
                Theta<-atan(r_object$data[,3]/r_object$data[,2])
                
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
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
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))
              
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
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              
              #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
              Theta<-rep(0,length(r_object$data[r_object$data$Outliers == "Outlier", ]$Outliers))
              
              for(i in 1:length(r_object$data[r_object$data$Outliers == "Outlier", ][,3])){
                Theta<-atan(r_object$data[r_object$data$Outliers == "Outlier", ][,3]/r_object$data[r_object$data$Outliers == "Outlier", ][,2])
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                B<- B + geom_segment(x = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], xend = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i],
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
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial")+ geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              for(i in 1:length(r_object$data[,3])){
                Theta<-atan(r_object$data[,3]/r_object$data[,2])
                
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
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
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))
              
              
              
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
              Theta<-rep(0,length(r_object$data[r_object$data$Outliers == "Outlier", ]$Outliers))
              
              for(i in 1:length(r_object$data[r_object$data$Outliers == "Outlier", ][,3])){
                Theta<-atan(r_object$data[r_object$data$Outliers == "Outlier", ][,3]/r_object$data[r_object$data$Outliers == "Outlier", ][,2])
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                B<- B + geom_segment(x = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], xend = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i],
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
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial")+ geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
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
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))
              
              #Draw lines indicating individual Q contributions and give text displaying each value
              for(i in 1:length(r_object$data[r_object$data$Outliers == "Outlier", ][,3])){
                
                
                
                
                B<- B + geom_segment(x = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], xend = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], y = r_object$data[r_object$data$Outliers == "Outlier", ][,3][i], yend =r_object$data[r_object$data$Outliers == "Outlier", ][,2][i]*r_object$coef[1],linetype="solid",colour="#56B4E9")
                
                
                
                
                B<- B + geom_text(x=r_object$data[r_object$data$Outliers == "Outlier", ][,2][i],y=(r_object$data[r_object$data$Outliers == "Outlier", ][,3][i]*0.9),label=round(r_object$data[r_object$data$Outliers == "Outlier", ][,4][i],digits=2),size=2.5)
                
              }
              
            }
            
          }
          
          #If axes scales should be fixed
          if(scale_match==FALSE){
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){
              
              #Produce plot showing full scale and all variants
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial")+ geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
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
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))
              
              #Draw lines indicating individual Q contributions and give text displaying each value
              for(i in 1:length(r_object$data[r_object$data$Outliers == "Outlier", ][,3])){
                
                
                
                
                B<- B + geom_segment(x = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], xend = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], y = r_object$data[r_object$data$Outliers == "Outlier", ][,3][i], yend =r_object$data[r_object$data$Outliers == "Outlier", ][,2][i]*r_object$coef[1],linetype="solid",colour="#56B4E9")
                
                
                
                
                B<- B + geom_text(x=r_object$data[r_object$data$Outliers == "Outlier", ][,2][i],y=(r_object$data[r_object$data$Outliers == "Outlier", ][,3][i]*0.9),label=round(r_object$data[r_object$data$Outliers == "Outlier", ][,4][i],digits=2),size=2.5)
                
              }
              
              
            }
            
          }
          
        }
        
      }
      
      #If only outliers are to be shown
      if(show_outliers==FALSE){
        
        #If the full scale is to be shown
        if(radial_scale==TRUE){
          
          #If axes scales should be fixed
          if(scale_match==TRUE){
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){
              
              #Produce plot showing full scale and all variants
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial")+ geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              for(i in 1:length(r_object$data[,3])){
                Theta<-atan(r_object$data[,3]/r_object$data[,2])
                
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
              
              print("No significant Outliers")
              
            }
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){
              
              #Plot circle for overall scale
              circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              
              circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing all variants, and IVW estimate using the full scale
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#56B4E9","#E69F00","#000000"))
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              #Draw lines showing the ratio estimate for each individual variant
              for(i in 1:length(r_object$data[,3])){
                Theta<-atan(r_object$data[,3]/r_object$data[,2])
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
            }
          }
          
          if(scale_match==FALSE){
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){
              
              #Produce plot showing full scale and all variants
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial")+ geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              for(i in 1:length(r_object$data[,3])){
                Theta<-atan(r_object$data[,3]/r_object$data[,2])
                
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
              
              print("No significant Outliers")
              
            }
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){
              
              #Plot circle for overall scale
              circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              
              circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing all variants, and IVW estimate using the full scale
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#56B4E9","#E69F00","#000000"))
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                Y_Scale<-c(Y_Scale,0)
                Y_Scale<-mround(Y_Scale)
                Y_Scale<-unique(Y_Scale)
                
                for (i in 1:length(Y_Scale)){
                  B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                  
                }
              }
              
              #Draw lines showing the ratio estimate for each individual variant
              for(i in 1:length(r_object$data[,3])){
                Theta<-atan(r_object$data[,3]/r_object$data[,2])
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
            }
          }
          
          
          
          
        }#radial scale
        
        #If the full scale should not be shown
        if(radial_scale==FALSE){
          
          #If axes scales should be fixed
          if(scale_match==TRUE){
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){
              
              #Produce plot showing full scale and all variants
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial")+ geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              
              print("No significant Outliers")
              
            }
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){
              
              #Plot circle for overall scale
              circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              
              circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing all variants, and IVW estimate using the full scale
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#56B4E9","#E69F00","#000000"))
              
            }
          }
          
          #If axes scales should be fixed
          if(scale_match==FALSE){
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){
              
              #Produce plot showing full scale and all variants
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj))+labs(title="IVW Radial")+ geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5),
                          label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling)-c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1]),
                          label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                 axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)+c(quantile(r_object$data$Wj,seq(0,0.1,0.1))[1])),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              
              print("No significant Outliers")
              
            }
            
            #If no outliers are present
            if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){
              
              #Plot circle for overall scale
              circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              
              circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing all variants, and IVW estimate using the full scale
              B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$coef[1],digits=3),size=4)+
                geom_text(x=cos(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$confint[1],digits=3),size=3)+
                geom_text(x=cos(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$confint[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#56B4E9","#E69F00","#000000"))
              
            }
          }
          
        }
      }
      
      return(B)
      
    }
    
    if(class(r_object)=="egger"){
      
      
      maximal<-atan(max(abs(r_object$data[,3]))/max(r_object$data[,2]))
      R.All<-max(abs(r_object$data[,3]- r_object$coef[1,1]))/sin(maximal)
      R.Egger<-R.All+(R.All/12)
      
      #Define scaling parameter
      Label.Scaling_Egg<-R.Egger/24
      
      Scalemin<-min(r_object$data[,3]/r_object$data[,2],r_object$confint[1])
      
      Scalemax<-max(r_object$data[,3]/r_object$data[,2],r_object$confint[2])
      
      #Plot circle for overall scale
      circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
      cxAll<-circledat_ALLEST$x
      cyAll<-circledat_ALLEST$y
      
      #Plot circle for Egger estimate
      circledat_EGGEREST <- circle.radial(c(0,r_object$coef[1,1]),R.Egger,length(r_object$data[,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
      cxEgger<-circledat_EGGEREST$x
      cyEgger<-circledat_EGGEREST$y
      
      
      #Define minimum and maximum range of values for y axis
      Y_MINEGGER<-min(  c(min(r_object$data[,3]),(sin(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg))+r_object$coef[1,1],r_object$coef[1,1]-Label.Scaling_Egg))
      Y_MAXEGGER<-max(  c(max(r_object$data[,3]),(sin(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg))+r_object$coef[1,1],r_object$coef[1,1]+Label.Scaling_Egg))
      
      #If all variants are to be shown
      if(show_outliers==TRUE){
        
        #If no outliers are present
        if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){
          
          #Generates a plot showing all variants not using the full scale
          B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
          )+labs(title="MR-Egger Radial") + geom_point(aes(colour="variant"))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(r_object$coef[2,1]))*(R.Egger+Label.Scaling_Egg), y=sin(atan(r_object$coef[2,1]))*(R.Egger+Label.Scaling_Egg) +
                        r_object$coef[1,1],label=round(r_object$coef[2,1],digits=3),size=4)+geom_text(x=cos(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                        r_object$coef[1,1],label=round(r_object$confint[1],digits=3),size=3)+geom_text(x=cos(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                                                                                                         r_object$coef[1,1],label=round(r_object$confint[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[2,1]))*R.Egger, y = r_object$coef[1,1], yend = (sin(atan(r_object$coef[2,1]))*R.Egger)+r_object$coef[1,1],colour="MR-Egger"))+
            scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINEGGER,Y_MAXEGGER))+scale_color_manual(breaks=c("variant","MR-Egger"),values=c("MR-Egger"="#D55E00","variant"="#000000"))+
            theme(legend.title=element_blank())
          
          #Indicates that no outliers are present
          print("No significant Outliers")
          
          
        }
        
        if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){
          
          
          maximal<-atan(max(abs(r_object$data[r_object$data$Outliers == "Outlier", ][,3]))/max(abs(r_object$data[r_object$data$Outliers == "Outlier", ][,2])))
          R.All<-max(abs(r_object$data[r_object$data$Outliers == "Outlier", ][,3]- r_object$coef[1,1]))/sin(maximal)
          R.Egger<-R.All+(R.All/12)
          
          #Define scaling parameter
          Label.Scaling_Egg<-R.Egger/24
          
          #Plot circle for Egger estimate
          circledat_EGGEREST <- circle.radial(c(0,r_object$coef[1,1]),R.Egger,length(r_object$data[r_object$data$Outliers == "Outlier", ][,3]),atan(r_object$confint[1]),atan(r_object$confint[2]))
          cxEgger<-circledat_EGGEREST$x
          cyEgger<-circledat_EGGEREST$y
          
          #Generate a plot showing only outlying variants and MR-Egger estimate without full scale
          B<-ggplot(r_object$data[r_object$data$Outliers == "Outlier", ],aes(x=Wj,y=BetaWj)
          )+labs(title="MR-Egger Radial") + geom_point(aes(colour=Outliers))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(r_object$coef[2,1]))*(R.Egger+Label.Scaling_Egg), y=sin(atan(r_object$coef[2,1]))*(R.Egger+Label.Scaling_Egg)+
                        r_object$coef[1,1],label=round(r_object$coef[2,1],digits=3),size=4)+geom_text(x=cos(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                        r_object$coef[1,1],label=round(r_object$confint[1],digits=3),size=3)+geom_text(x=cos(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                                                                                                         r_object$coef[1,1],label=round(r_object$confint[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[2,1]))*R.Egger, y = r_object$coef[1,1], yend = (sin(atan(r_object$coef[2,1]))*R.Egger)+r_object$coef[1,1],colour="MR-Egger"))+
            scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+
            scale_y_continuous(limits = c(Y_MINEGGER,Y_MAXEGGER))+scale_colour_manual(breaks=c("Variant","Outlier","MR-Egger"),values=c("#D55E00","#E69F00","#000000"))+
            theme(legend.title=element_blank())
          
          #Draw lines indicating individual Q contributions and give text displaying each value
          for(i in 1:length(r_object$data[r_object$data$Outliers == "Outlier", ][,3])){
            B<- B + geom_segment(x = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], xend = r_object$data[r_object$data$Outliers == "Outlier", ][,2][i], y = r_object$data[r_object$data$Outliers == "Outlier", ][,3][i], yend =r_object$data[r_object$data$Outliers == "Outlier", ][,2][i]*r_object$coef[2,1]+r_object$coef[1,1],linetype="solid",colour="#D55E00")
            B<- B + geom_text(x=r_object$data[r_object$data$Outliers == "Outlier", ][,2][i],y=(r_object$data[r_object$data$Outliers == "Outlier", ][,3][i]*0.9),label=round(r_object$data[r_object$data$Outliers == "Outlier", ][,4][i],digits=2),size=2.5)
            
          }
          
        }
        
      }
      
      
      
      #If all variants are to be shown
      if(show_outliers==FALSE){
        
        
        #If no outliers are present
        if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])==0){
          
          
          #Generates a plot showing all variants not using the full scale
          B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
          )+labs(title="MR-Egger Radial") + geom_point(aes(colour="variant"))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(r_object$coef[2,1]))*(R.Egger+Label.Scaling_Egg), y=sin(atan(r_object$coef[2,1]))*(R.Egger+Label.Scaling_Egg) +
                        r_object$coef[1,1],label=round(r_object$coef[2,1],digits=3),size=4)+geom_text(x=cos(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                        r_object$coef[1,1],label=round(r_object$confint[1],digits=3),size=3)+geom_text(x=cos(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                                                                                                         r_object$coef[1,1],label=round(r_object$confint[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[2,1]))*R.Egger, y = r_object$coef[1,1], yend = (sin(atan(r_object$coef[2,1]))*R.Egger)+r_object$coef[1,1],colour="MR-Egger"))+
            scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINEGGER,Y_MAXEGGER))+scale_color_manual(breaks=c("variant","MR-Egger"),values=c("MR-Egger"="#D55E00","variant"="#000000"))+
            theme(legend.title=element_blank())
          
          #Indicates that no outliers are present
          print("No significant Outliers")
        }
        
        #If no outliers are present
        if(nrow(r_object$data[r_object$data$Outliers == "Outlier", ])>0){
          
          #Generate a plot showing all variants not using full scale
          B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
          )+labs(title="MR-Egger Radial") + geom_point(aes(colour=Outliers))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(r_object$coef[2,1]))*(R.Egger+Label.Scaling_Egg), y=sin(atan(r_object$coef[2,1]))*(R.Egger+Label.Scaling_Egg)+
                        r_object$coef[1,1],label=round(r_object$coef[2,1],digits=3),size=4)+geom_text(x=cos(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(r_object$confint[1]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                        r_object$coef[1,1],label=round(r_object$confint[1],digits=3),size=3)+geom_text(x=cos(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(r_object$confint[2]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                                                                                                         r_object$coef[1,1],label=round(r_object$confint[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(r_object$coef[2,1]))*R.Egger, y = r_object$coef[1,1], yend = (sin(atan(r_object$coef[2,1]))*R.Egger)+r_object$coef[1,1],colour="MR-Egger"))+
            scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINEGGER,Y_MAXEGGER))+scale_colour_manual(breaks=c("Variant","Outlier","MR-Egger"),values=c("#D55E00","#E69F00","#000000"))+
            theme(legend.title=element_blank())
          
        }
        
      }
      
      
    }
    
  }
  
  if(length(r_object)==19){
    
    names(r_object)<-c("IVW.coef","IVW.qstatistic","IVW.df","IVW.outliers","data","IVW.confint","it.coef","Fex.coef",
                       "Rex.coef","It.confint","Fex.confint","Rex.confint","meanF","egger.coef","egger.qstatistic","egger.df","egger.outliers","egger.data","egger.confint")
    
    
    #Define scale, IVW and MR-Egger circle radii
    maximal<-atan(max(abs(r_object$data[,3]))/max(abs(r_object$data[,2])))
    R.All<-max(abs(r_object$data[,3]- r_object$egger.coef[1,1]))/sin(maximal)+median(r_object$data[,2])
    R.IVW<-R.All+(R.All/12)
    R.Egger<-R.All+(R.All/4)
    
    #Define general scaling parameter
    Label.Scaling<-R.IVW/24
    
    #Define scaling parameter for MR-Egger
    Label.Scaling_Egg<-R.Egger/8
    
    Scalemin<-min(r_object$data[,3]/r_object$data[,2],r_object$IVW.confint[1],r_object$egger.confint[1])
    
    Scalemax<-max(r_object$data[,3]/r_object$data[,2],r_object$IVW.confint[2],r_object$egger.confint[2])
    
    #Define minimum and maximum values of y axis for scaling
    Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling_Egg)
    Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling_Egg)
    
    Y_Range<-Scalemax-Scalemin
    
    #Function for rounding to nearest 0.05 increment
    mround <- function(x){ 
      0.5*round(x/0.5) 
    }
    
    #Define circle plotting points for scale, IVW and MR-Egger estimates
    circledat_EGGEREST <- circle.radial(c(0, r_object$egger.coef[1,1]),R.Egger,length(r_object$data[,3]),atan(r_object$egger.confint[1]),atan(r_object$egger.confint[2]))
    cxEgger<-circledat_EGGEREST$x
    cyEgger<-circledat_EGGEREST$y
    
    #Plot circle for overall scale
    circledat_ALLEST <- circle.radial(c(0,0),R.All,length(r_object$data[,3]),atan(Scalemin),atan(Scalemax))
    cxAll<-circledat_ALLEST$x
    cyAll<-circledat_ALLEST$y
    
    circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(r_object$data[,3]),atan(r_object$IVW.confint[1]),atan(r_object$IVW.confint[2]))
    cxIVW<-circledat_IVWEST$x
    cyIVW<-circledat_IVWEST$y
    
    which(grepl("Outlier", r_object$egger.data[,6]))
    
    r_object$data$out_types<-rep(0,length(r_object$data[,6]))
    
    for(i in 1:length(r_object$data[,3])){
      
      if(r_object$data[,6][i]=="Variant" & r_object$egger.data[,6][i]=="Variant"){
        
        r_object$data$out_types[i]<-"Variant"
        
      }
      
      if(r_object$data[,6][i]=="Variant" & r_object$egger.data[,6][i]=="Outlier"){
        
        r_object$data$out_types[i]<-"MR-Egger Outlier"
        
      }
      
      if(r_object$data[,6][i]=="Outlier" & r_object$egger.data[,6][i]=="Variant"){
        
        r_object$data$out_types[i]<-"IVW Outlier"
        
      }
      
      if(r_object$data[,6][i]=="Outlier" & r_object$egger.data[,6][i]=="Outlier"){
        
        r_object$data$out_types[i]<-"IVW and MR-Egger Outlier"
        
      }
      
    }
    
    r_object$data$out_types<-factor(r_object$data$out_types)
    
    #If only outliers are to be shown
    if(show_outliers==TRUE){
      
      #If the full scale is to be shown
      if(radial_scale==TRUE){
        
        #If axes scales should be fixed
        if(scale_match==TRUE){
          
          #If no outliers are present
          if(nrow(r_object$data[r_object$data$out_types == "Variant", ])==length(r_object$data[,3])){
            
            #Generate a plot showing both IVW and MR-Egger estimates, using the full scale and all variants
            B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant"))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+coord_fixed()+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+
                                 r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+ theme(legend.title=element_blank())                   
            
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            for(i in 1:length(r_object$data[,3])){
              Theta<-atan(r_object$data[,3]/r_object$data[,2])
              
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
            
            print("No significant Outliers")
            
            
          }
          
          if(nrow(r_object$data[r_object$data$out_types == "Variant", ])<length(r_object$data[,3])){
            
            temp.dat<-r_object$data[unique(c(which(grepl("Outlier", r_object$egger.data[,6])),
                                             which(grepl("Outlier", r_object$data[,6])))), ]
            
            #Redefine IVW, MR-Egger and scale radii
            maximal<-atan(max(abs(temp.dat[,3]))/max(abs(temp.dat[,2])))
            R.All<-max(abs(r_object$data[,3]- r_object$egger.coef[1,1]))/sin(maximal)+median(r_object$data[,2])
            R.IVW<-R.All+(R.All/12)
            R.Egger<-R.All+(R.All/4)
            
            #Redefine general label scaling parameter
            Label.Scaling<-R.IVW/24
            
            #Redefine MR-Egger label scaling parameter
            Label.Scaling_Egg<-R.Egger/8
            
            #Redefine values for plotting IVW, MR-Egger and scale circles
            circledat_EGGEREST <- circle.radial(c(0, r_object$egger.coef[1,1]),R.Egger,length(temp.dat[,3]),atan(r_object$egger.confint[1]),atan(r_object$egger.confint[2]))
            cxEgger<-circledat_EGGEREST$x
            cyEgger<-circledat_EGGEREST$y
            
            #Plot circle for overall scale
            circledat_ALLEST <- circle.radial(c(0,0),R.All,length(temp.dat[,3]),atan(Scalemin),atan(Scalemax))
            cxAll<-circledat_ALLEST$x
            cyAll<-circledat_ALLEST$y
            
            circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(temp.dat[,3]),atan(r_object$IVW.confint[1]),atan(r_object$IVW.confint[2]))
            cxIVW<-circledat_IVWEST$x
            cyIVW<-circledat_IVWEST$y
            
            #Define minimum and maximum values of y axis for scaling
            Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
            Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
            
            Y_Range<-Scalemax-Scalemin
            
            #Generates a plot showing the MR-Egger and IVW estimates, showing only outliers and using the full scale
            B<-ggplot(temp.dat,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+coord_fixed()+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="white","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())
            
            if(max(abs(temp.dat[,3]/temp.dat[,2]))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            if(max(abs(temp.dat[,3]/temp.dat[,2]))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
            Theta<-rep(0,length(temp.dat[,3]))
            
            for(i in 1:length(temp.dat[,3])){
              Theta<-atan(temp.dat[,3]/temp.dat[,2])
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
            
            
            
          }
          
        }
        
        #If axes scales should be fixed
        if(scale_match==FALSE){
          
          #If no outliers are present
          if(nrow(r_object$data[r_object$data$out_types == "Variant", ])==length(r_object$data[,3])){
            
            
            #Generate a plot showing both IVW and MR-Egger estimates, using the full scale and all variants
            B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant"))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+
                                 r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+ theme(legend.title=element_blank())                   
            
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            for(i in 1:length(r_object$data[,3])){
              Theta<-atan(r_object$data[,3]/r_object$data[,2])
              
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
            
            print("No significant Outliers")
            
            
          }
          
          if(nrow(r_object$data[r_object$data$out_types == "Variant", ])<length(r_object$data[,3])){
            
            temp.dat<-r_object$data[unique(c(which(grepl("Outlier", r_object$egger.data[,6])),
                                             which(grepl("Outlier", r_object$data[,6])))), ]
            
            #Redefine IVW, MR-Egger and scale radii
            maximal<-atan(max(abs(temp.dat[,3]))/max(abs(temp.dat[,2])))
            R.All<-max(abs(r_object$data[,3]- r_object$egger.coef[1,1]))/sin(maximal)+median(r_object$data[,2])
            R.IVW<-R.All+(R.All/12)
            R.Egger<-R.All+(R.All/4)
            
            #Redefine general label scaling parameter
            Label.Scaling<-R.IVW/24
            
            #Redefine MR-Egger label scaling parameter
            Label.Scaling_Egg<-R.Egger/8
            
            #Redefine values for plotting IVW, MR-Egger and scale circles
            circledat_EGGEREST <- circle.radial(c(0, r_object$egger.coef[1,1]),R.Egger,length(temp.dat[,3]),atan(r_object$egger.confint[1]),atan(r_object$egger.confint[2]))
            cxEgger<-circledat_EGGEREST$x
            cyEgger<-circledat_EGGEREST$y
            
            #Plot circle for overall scale
            circledat_ALLEST <- circle.radial(c(0,0),R.All,length(temp.dat[,3]),atan(Scalemin),atan(Scalemax))
            cxAll<-circledat_ALLEST$x
            cyAll<-circledat_ALLEST$y
            
            circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(temp.dat[,3]),atan(r_object$IVW.confint[1]),atan(r_object$IVW.confint[2]))
            cxIVW<-circledat_IVWEST$x
            cyIVW<-circledat_IVWEST$y
            
            #Define minimum and maximum values of y axis for scaling
            Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
            Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
            
            Y_Range<-Scalemax-Scalemin
            
            #Generates a plot showing the MR-Egger and IVW estimates, showing only outliers and using the full scale
            B<-ggplot(temp.dat,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="white","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())
            
            if(max(abs(temp.dat[,3]/temp.dat[,2]))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            if(max(abs(temp.dat[,3]/temp.dat[,2]))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
            Theta<-rep(0,length(temp.dat[,3]))
            
            for(i in 1:length(temp.dat[,3])){
              Theta<-atan(temp.dat[,3]/temp.dat[,2])
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
            
            
            
          }
          
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
      }
      
      #If the full scale is to be shown
      if(radial_scale==FALSE){
        
        #If no outliers are present
        if(nrow(r_object$data[r_object$data$out_types == "Variant", ])==length(r_object$data[,3])){
          
          MINCILIM<-min(sin(atan(r_object$IVW.confint[1]))*R.IVW + Label.Scaling,sin(atan(r_object$egger.confint[1]-1))*(R.Egger +Label.Scaling_Egg-r_object$egger.coef[1,1]),r_object$egger.coef[1,1])
          
          MAXCILIM<-max(sin(atan(r_object$IVW.confint[2]))*R.IVW + Label.Scaling,sin(atan(r_object$egger.confint[2]+1))*(R.Egger +Label.Scaling_Egg-r_object$egger.coef[1,1]),r_object$egger.coef[1,1])
          
          if(min(r_object$data[,3]) < MINCILIM){
            Y_MIN<-min(r_object$data[,3])
          } 
          
          if(min(r_object$data[,3]) > MINCILIM){
            Y_MIN<-MINCILIM
          } 
          
          if(max(r_object$data[,3]) >MAXCILIM){
            Y_MAX<-max(r_object$data[,3])
          } 
          
          if(max(r_object$data[,3]) <MAXCILIM){
            Y_MAX<-MAXCILIM
          } 
          
          
          Y_Range<-Y_MAX-Y_MIN
          
          #Generate a plot showing the MR-Egger and IVW estimates, showing all variants and not using a full scale
          B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
          )+labs(title="Radial Estimates")+ geom_point(aes(colour=Outliers)) +geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
            geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
            geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
            geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
            geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
            geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
            geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
            theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
            scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+
                                                                                                                r_object$egger.coef[1,1],colour="MR-Egger"))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+scale_color_manual(breaks=c("Variant","Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000","Outlier"="#E69F00"))+
            theme(legend.title=element_blank())
          
          #Indicate that no outliers are present
          print("No significant Outliers")
          
        }
        
        if(nrow(r_object$data[r_object$data$out_types == "Variant", ])<length(r_object$data[,3])){
          
          temp.dat<-r_object$data[unique(c(which(grepl("Outlier", r_object$egger.data[,6])),
                                           which(grepl("Outlier", r_object$data[,6])))), ]
          
          #Redefine IVW, MR-Egger and scale radii
          maximal<-atan(max(abs(temp.dat[,3]))/max(abs(temp.dat[,2])))
          R.All<-max(abs(r_object$data[,3]- r_object$egger.coef[1,1]))/sin(maximal)+median(r_object$data[,2])
          R.IVW<-R.All+(R.All/12)
          R.Egger<-R.All+(R.All/4)
          
          #Redefine general label scaling parameter
          Label.Scaling<-R.IVW/24
          
          #Redefine MR-Egger label scaling parameter
          Label.Scaling_Egg<-R.Egger/8
          
          #Redefine values for plotting IVW, MR-Egger and scale circles
          circledat_EGGEREST <- circle.radial(c(0, r_object$egger.coef[1,1]),R.Egger,length(temp.dat[,3]),atan(r_object$egger.confint[1]),atan(r_object$egger.confint[2]))
          cxEgger<-circledat_EGGEREST$x
          cyEgger<-circledat_EGGEREST$y
          
          #Plot circle for overall scale
          circledat_ALLEST <- circle.radial(c(0,0),R.All,length(temp.dat[,3]),atan(Scalemin),atan(Scalemax))
          cxAll<-circledat_ALLEST$x
          cyAll<-circledat_ALLEST$y
          
          circledat_IVWEST <- circle.radial(c(0,0),R.IVW,length(temp.dat[,3]),atan(r_object$IVW.confint[1]),atan(r_object$IVW.confint[2]))
          cxIVW<-circledat_IVWEST$x
          cyIVW<-circledat_IVWEST$y
          
          MINCILIM<-min(sin(atan(r_object$IVW.confint[1]))*R.IVW + Label.Scaling,sin(atan(r_object$egger.confint[1]-1))*(R.Egger +Label.Scaling_Egg-r_object$egger.coef[1,1]),r_object$egger.coef[1,1])
          
          MAXCILIM<-max(sin(atan(r_object$IVW.confint[2]))*R.IVW + Label.Scaling,sin(atan(r_object$egger.confint[2]+1))*(R.Egger +Label.Scaling_Egg-r_object$egger.coef[1,1]),r_object$egger.coef[1,1])
          
          if(min(temp.dat[,3]) < MINCILIM){
            Y_MIN<-min(temp.dat[,3])
          } 
          
          if(min(temp.dat[,3]) > MINCILIM){
            Y_MIN<-MINCILIM
          } 
          
          if(max(temp.dat[,3]) >MAXCILIM){
            Y_MAX<-max(temp.dat[,3])
          } 
          
          if(max(temp.dat[,3]) <MAXCILIM){
            Y_MAX<-MAXCILIM
          } 
          
          Y_Range<-Y_MAX-Y_MIN
          
          #Generate a plot showing both IVW and MR-Egger estimates, only showing outlying variants and not using the full scale
          B<-ggplot(temp.dat,aes(x=Wj,y=BetaWj)
          )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types)) +geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
            geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
            geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
            geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
            geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
            geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
            geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
            scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+r_object$egger.coef[1,1],colour="MR-Egger"))+
            scale_y_continuous(limits = c(Y_MIN,Y_MAX))+scale_color_manual(breaks=c("IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="white","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())                   
          
          #Draw lines showing Q statistic contributions for each variant using respective colours
          for(i in 1:length(temp.dat[,3])){
            if(abs((temp.dat[,3][i]-temp.dat[,2][i]*r_object$IVW.coef[1]))>abs((temp.dat[,3][i]-(temp.dat[,2][i]*r_object$egger.coef[2,1]+r_object$egger.coef[1,1])))){
              
              B<- B + geom_segment(x = temp.dat[,2][i], xend = temp.dat[,2][i], y = temp.dat[,3][i], yend =temp.dat[,2][i]*r_object$IVW.coef[1],linetype="solid",colour="#56B4E9")
              
              B<- B + geom_segment(x = temp.dat[,2][i], xend = temp.dat[,2][i], y = temp.dat[,3][i], yend =temp.dat[,2][i]*r_object$egger.coef[2,1]+r_object$egger.coef[1,1],linetype="solid",colour="#D55E00")
              
            }
            
            if(abs((temp.dat[,3][i]-temp.dat[,2][i]*r_object$IVW.coef[1]))<abs((temp.dat[,3][i]-(temp.dat[,2][i]*r_object$egger.coef[2,1]+r_object$egger.coef[1,1])))){
              
              B<- B + geom_segment(x = temp.dat[,2][i], xend = temp.dat[,2][i], y = temp.dat[,3][i], yend =temp.dat[,2][i]*r_object$egger.coef[2,1]+r_object$egger.coef[1,1],linetype="solid",colour="#D55E00")
              
              B<- B + geom_segment(x = temp.dat[,2][i], xend = temp.dat[,2][i], y = temp.dat[,3][i], yend =temp.dat[,2][i]*r_object$IVW.coef[1],linetype="solid",colour="#56B4E9")
              
            }
            
          }
          
        }
        
      }
      
    }
    
    #If only outliers are to be shown
    if(show_outliers==FALSE){
      
      #If the full scale is to be shown
      if(radial_scale==TRUE){
        
        #If axes scales should be fixed
        if(scale_match==TRUE){
          
          #If no outliers are present
          if(nrow(r_object$data[r_object$data$out_types == "Variant", ])==length(r_object$data[,3])){
            
            
            #Generate a plot showing both IVW and MR-Egger estimates, using the full scale and all variants
            B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant"))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+coord_fixed()+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+
                                 r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+ theme(legend.title=element_blank())                   
            
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            for(i in 1:length(r_object$data[,3])){
              Theta<-atan(r_object$data[,3]/r_object$data[,2])
              
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
            
            print("No significant Outliers")
            
            
          }
          
          if(nrow(r_object$data[r_object$data$out_types == "Variant", ])<length(r_object$data[,3])){
            
            #Generates a plot showing the MR-Egger and IVW estimates, showing only outliers and using the full scale
            B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+coord_fixed()+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="black","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())
            
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
            Theta<-rep(0,length(r_object$data[,3]))
            
            for(i in 1:length(r_object$data[,3])){
              Theta<-atan(r_object$data[,3]/r_object$data[,2])
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
            
            
            
          }
          
        }
        
        #If axes scales should be fixed
        if(scale_match==FALSE){
          
          #If no outliers are present
          if(nrow(r_object$data[r_object$data$out_types == "Variant", ])==length(r_object$data[,3])){
            
            #Generate a plot showing both IVW and MR-Egger estimates, using the full scale and all variants
            B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant"))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+
                                 r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+ theme(legend.title=element_blank())                   
            
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            for(i in 1:length(r_object$data[,3])){
              Theta<-atan(r_object$data[,3]/r_object$data[,2])
              
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
            
            print("No significant Outliers")
            
            
          }
          
          if(nrow(r_object$data[r_object$data$out_types == "Variant", ])<length(r_object$data[,3])){
            
            #Generates a plot showing the MR-Egger and IVW estimates, showing only outliers and using the full scale
            B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
              geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
              geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
              geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="black","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())
            
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            if(max(abs(r_object$data[,3]/r_object$data[,2]))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              Y_Scale<-c(Y_Scale,0)
              Y_Scale<-mround(Y_Scale)
              Y_Scale<-unique(Y_Scale)
              
              for (i in 1:length(Y_Scale)){
                B<- B + geom_text(x=(cos(atan(Y_Scale[i]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[i]))*(R.All+Label.Scaling)), label=Y_Scale[i],size=2.5)
                
              }
            }
            
            #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
            Theta<-rep(0,length(r_object$data[,3]))
            
            for(i in 1:length(r_object$data[,3])){
              Theta<-atan(r_object$data[,3]/r_object$data[,2])
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
            
            
            
          }
          
        }
        
      }
      
      #If the full scale is to be shown
      if(radial_scale==FALSE){
        
        #If no outliers are present
        if(nrow(r_object$data[r_object$data$out_types == "Variant", ])==length(r_object$data[,3])){
          
          MINCILIM<-min(sin(atan(r_object$IVW.confint[1]))*R.IVW + Label.Scaling,sin(atan(r_object$egger.confint[1]-1))*(R.Egger +Label.Scaling_Egg-r_object$egger.coef[1,1]),r_object$egger.coef[1,1])
          
          MAXCILIM<-max(sin(atan(r_object$IVW.confint[2]))*R.IVW + Label.Scaling,sin(atan(r_object$egger.confint[2]+1))*(R.Egger +Label.Scaling_Egg-r_object$egger.coef[1,1]),r_object$egger.coef[1,1])
          
          if(min(r_object$data[,3]) < MINCILIM){
            Y_MIN<-min(r_object$data[,3])
          } 
          
          if(min(r_object$data[,3]) > MINCILIM){
            Y_MIN<-MINCILIM
          } 
          
          if(max(r_object$data[,3]) >MAXCILIM){
            Y_MAX<-max(r_object$data[,3])
          } 
          
          if(max(r_object$data[,3]) <MAXCILIM){
            Y_MAX<-MAXCILIM
          } 
          
          Y_Range<-Y_MAX-Y_MIN
          
          #Generate a plot showing the MR-Egger and IVW estimates, showing all variants and not using a full scale
          B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
          )+labs(title="Radial Estimates")+ geom_point(aes(colour=Outliers)) +geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
            geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
            geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
            geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
            geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
            geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
            geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
            theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
            scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+
                                                                                                                r_object$egger.coef[1,1],colour="MR-Egger"))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+scale_color_manual(breaks=c("Variant","Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000","Outlier"="#E69F00"))+
            theme(legend.title=element_blank())
          
          #Indicate that no outliers are present
          print("No significant Outliers")
          
        }
        
        if(nrow(r_object$data[r_object$data$out_types == "Variant", ])<length(r_object$data[,3])){
          
          #Generates a plot showing the MR-Egger and IVW estimates, showing only outliers and using the full scale
          B<-ggplot(r_object$data,aes(x=Wj,y=BetaWj)
          )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
            geom_text(x=cos(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), y=sin(atan(r_object$IVW.coef[1]))*(R.IVW+Label.Scaling*1.5), label=round(r_object$IVW.coef[1],digits=3),size=4)+
            geom_text(x=cos(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[1])-0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[1],digits=3),size=3)+
            geom_text(x=cos(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(r_object$IVW.confint[2])+0.02)*(R.IVW+Label.Scaling),label=round(r_object$IVW.confint[2],digits=3),size=3)+
            theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                              axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(r_object$IVW.coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$IVW.coef[1]))*R.IVW,colour="IVW"))+
            scale_x_continuous(limits = c(0,max(cxEgger)+Label.Scaling_Egg*2),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
            geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5), y=sin(atan(r_object$egger.coef[2,1]))*(R.Egger+Label.Scaling*1.5) +r_object$egger.coef[1,1],label=round(r_object$egger.coef[2,1],digits=3),size=4)+
            geom_text(x=cos(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[1]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[1],digits=3),size=3)+
            geom_text(x=cos(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(r_object$egger.confint[2]))*(R.Egger+Label.Scaling*1.5))+r_object$egger.coef[1,1],label=round(r_object$egger.confint[2],digits=3),size=3)+
            geom_segment(aes(x = 0, xend = cos(atan(r_object$egger.coef[2,1]))*R.Egger, y = r_object$egger.coef[1,1], yend = (sin(atan(r_object$egger.coef[2,1]))*R.Egger)+r_object$egger.coef[1,1],colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="black","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())
          
          
        }
        
      }
      
    }
    
  }
  
  return(B)
  
}