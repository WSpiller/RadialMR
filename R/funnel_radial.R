#' funnel_radial
#'
#' A function for producing generalized radial IVW and MR-Egger funnel plots either individually or simultaneously. The function also allows for lines indicating the magnitude for the MR Egger transformation for each variant, though it should be noted that this distance is a function of the weight attributed to the variant, and is therefore not indicative of outliers.
#'
#' @param r_object An object of class \code{"IVW"} or \code{"egger"}. For visualising both estimates simultaneously, both objects should be included as a vector c(\code{A},\code{B}), where \code{A} and \code{B} denote the \code{"IVW"} and \code{"egger"} objects respectively.
#' @param show_transform Indicates whether to produce a plot including lines showing the magnitude of the MR Egger transformation for each variant (\code{TRUE}), or a scatterplot showing only the variants and corresponding estimates (\code{FALSE}).
#' @return A ggplot object containing a radial funnel plot of either the IVW, MR-Egger, or both estimates simultaneously.
#'@author Wes Spiller; Jack Bowden.
#'@references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.
#'@export
#'@examples
#'
#'funnel_radial(r_object,TRUE) 

funnel_radial<-function(r_object,show_transform){
  
  if(missing(show_transform)) {
    show_transform<-FALSE
  }
  
  
  #Load ggplot2 library
  library(ggplot2)
  
  if(length(r_object)>=6 && length(r_object)<=13){
    
    if(class(r_object)=="IVW"){
      
      #Produce plot showing full scale and all variants
      B<-ggplot(r_object$data,aes(x=r_object$data$BetaWj/r_object$data$Wj,y=Wj))+labs(title="Radial Funnel Plot")+ geom_point(aes(colour="Variant"))+
        theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black"))+xlab(expression(hat(beta)[j]))+ylab(expression(sqrt(W[j])))+
        geom_segment(aes(x = r_object$coef[1], xend = r_object$coef[1], y = 0, yend = max(r_object$data$Wj),colour="IVW Estimate"))+
        geom_segment(aes(x = r_object$confint[1], xend = r_object$confint[2], y = 0, yend = 0,colour="IVW Estimate"))+
        scale_y_continuous(limits = c(-0.2,max(r_object$data$Wj)))+
        scale_color_manual(name="",breaks=c("Variant","IVW Estimate"),values=c("Variant"="black","IVW Estimate"="#56B4E9"))+
        geom_text(x=r_object$coef[1], y=-0.2,
                  label=paste(round(r_object$coef[1],digits=3),paste("(",round(r_object$confint[1],digits=2),",",round(r_object$confint[2],digits=2),")"),collapse=""),size=3)
      
      #     geom_text(x=r_object$coef[1], y=-0.2,
      #              label=round(r_object$coef[1],digits=3),size=3)
    }
    
    if(class(r_object)=="egger"){
      
      transeg<-(r_object$data$BetaWj/r_object$data$Wj)-(r_object$coef[1,1]/r_object$data$Wj)
      Wj<-r_object$data$Wj
      temp<-data.frame(transeg,Wj)
      
      #Produce plot showing full scale and all variants
      B<-ggplot(temp,aes(x=transeg,y=Wj))+labs(title="Radial Funnel Plot")+ geom_point(aes(colour="Variant"))+
        theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                         axis.line = element_line(colour = "black"))+xlab(expression(hat(beta)[j]))+ylab(expression(sqrt(W[j])))+
        geom_segment(aes(x = r_object$coef[2,1], xend = r_object$coef[2,1], y = 0, yend = max(r_object$data$Wj),colour="MR Egger Estimate"))+
        geom_segment(aes(x = r_object$confint[1], xend = r_object$confint[2], y = 0, yend = 0,colour="MR Egger Estimate"))+
        scale_y_continuous(limits = c(-0.2,max(r_object$data$Wj)))+
        scale_color_manual(name="",breaks=c("Variant","MR Egger Estimate"),values=c("Variant"="black","MR Egger Estimate"="#D55E00"))+
        geom_text(x=r_object$coef[2,1], y=-0.2,
                  label=paste(round(r_object$coef[2,1],digits=3),paste("(",round(r_object$confint[1],digits=2),",",round(r_object$confint[2],digits=2),")"),collapse=""),size=3)
    }
    
  }
  
  if(length(r_object)==19){
    
    names(r_object)<-c("IVW.coef","IVW.qstatistic","IVW.df","IVW.outliers","data","IVW.confint","it.coef","Fex.coef",
                       "Rex.coef","It.confint","Fex.confint","Rex.confint","meanF","egger.coef","egger.qstatistic","egger.df","egger.outliers","egger.data","egger.confint")
    
    transeg<-(r_object$data$BetaWj/r_object$data$Wj)-(r_object$egger.coef[1,1]/r_object$data$Wj)
    
    betavec<-c(r_object$data$BetaWj/r_object$data$Wj,transeg)
    Wj<-c(r_object$data$Wj,r_object$data$Wj)
    tran.indicator<-c(rep("Variant",length(betavec)/2),rep("Egger transform",length(betavec)/2))
    
    temp<-data.frame(betavec,Wj,tran.indicator)
    
    
    #Produce plot showing full scale and all variants
    B <-ggplot(temp,aes(x=betavec,y=Wj))+labs(title="Radial Funnel Plot")+ geom_point(aes(colour=tran.indicator))+
      theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"))+xlab(expression(hat(beta)[j]))+ylab(expression(sqrt(W[j])))+
      geom_segment(aes(x = r_object$IVW.coef[1], xend = r_object$IVW.coef[1], y = -.5, yend = max(r_object$data$Wj),colour="IVW Estimate"))+
      geom_segment(aes(x = r_object$egger.coef[2,1], xend = r_object$egger.coef[2,1], y = 0, yend = max(r_object$data$Wj),colour="Egger Estimate"))+
      geom_segment(aes(x = r_object$egger.confint[1], xend = r_object$egger.confint[2], y = 0, yend = 0,colour="Egger Estimate"))+
      geom_segment(aes(x = r_object$IVW.confint[1], xend = r_object$IVW.confint[2], y = -.5, yend = -.5,colour="IVW Estimate"))+
      scale_y_continuous(limits = c(-1,max(r_object$data$Wj)))+
      scale_color_manual(name="",breaks=c("Variant","Egger transform","IVW Estimate","Egger Estimate"),values=c("Variant"="black","Egger transform"="#D55E00","IVW Estimate"="#56B4E9","Egger Estimate"="#D55E00"))+
      geom_text(x=r_object$IVW.coef[1], y=-0.7,label=paste(round(r_object$IVW.coef[1],digits=3),paste("(",round(r_object$IVW.confint[1],digits=2),",",round(r_object$IVW.confint[2],digits=2),")"),collapse=""),size=3)+
      geom_text(x=r_object$egger.coef[2,1], y=-0.2,label=paste(round(r_object$egger.coef[2,1],digits=3),paste("(",round(r_object$egger.confint[1],digits=2),",",round(r_object$egger.confint[2],digits=2),")"),collapse=""),size=3)
    
    
    if(show_transform==TRUE){
      
      for(i in 1:length(transeg)){
        
        B <- B + geom_segment(x = temp$betavec[i], xend = temp$betavec[i+length(transeg)], y = temp$Wj[i], yend = temp$Wj[i],linetype="dotted")
        
      }
      
    }
    
  }
  
  
  return(B)
  
}












