#' RadialMR
#'
#' Fits radial forms of Inverse Variance Weighted (IVW) and MR-Egger regression. It provides output summaries for each method, and plotting options for visualising goodness-of-fit using Cochran's Q statistic. This allows for identification of outlying observations and is coding invariant unlike conventional MR-Egger regression.
#'
#' @param BXG A numeric vector of beta-coefficient values for genetic associations with the first variable (exposure).
#' @param BYG A numeric vector of beta-coefficient values for genetic associations with the second variable (outcome).
#' @param seBXG The standard errors corresponding to the beta-coefficients BXG.
#' @param seBYG The standard errors corresponding to the beta-coefficients BYG.
#' @param RSID A vector of names for the genetic variants included in the analysis. If unavailable, \code{RSID="NO"} will generate placeholder ID numbers for each variant.
#' @param METHOD Indicates whether to perform radial IVW (\code{"IVW"}),radial MR-Egger (\code{"EGGER"}), or both radial methods (\code{"BOTH"}) simultaneously.
#' @param SCALE Indicates whether to produce a plot including a full radial scale (\code{"YES"}), or a scatterplot showing only the effect estimates (\code{"NO"}).
#' @param OUTLIERS Indicates whether display only the set of variants identified as outliers (\code{"YES"}) or the complete set of variants (\code{"NO"}). Note that when (\code{"OUTLIERS=YES"}), non-outlying variants further from the origin than the furthest outlier will cause an error message that one or more points have been omitted. These are non-outlying variants beyond the scale.
#' @param ALPHA A value specifying the statistical significance threshold for identifying outliers (\code{"0.05"} specifies a p-value threshold of 0.05)
#' @param SCALEMATCH Indicates whether x and y axes should have the same range(\code{"YES"}), or different ranges (\code{"NO"}).
#' @return An object containing at least the following components:\describe{
#' \item{\code{plot}}{A radial plot showing one or more radial MR estimates.}
#' \item{\code{Q}}{A data frame giving the overall Q statistic for the model method and respective chi-squared test. If \code{METHOD="BOTH"}, separate Q statistics and chi-square values are given for each model fit.}
#' \item{\code{outliers}}{A data frame containing variants identified as outliers, with respective Q statistics, chi-squared tests and SNP identification. In cases where \code{METHOD="BOTH"}, the dataframe will contain variants which are identified as an outlier for either model.}
#' \item{\code{IVW}}{A summary output for the radial IVW regression fit. If \code{METHOD="EGGER"}, will return \code{"NULL"}.}
#' \item{\code{egger}}{A summary output for the radial MR-Egger regression fit If \code{METHOD="IVW"}, will return \code{NULL}.}
#'}
#'@author Wes Spiller; Jack Bowden.
#'@references Bowden, J. et al., 2017. "Improving the visualisation, interpretation and analysis of two-sample summary data Mendelian randomization via the radial plot and radial regression." BioRxiv. doi: https://doi.org/10.1101/200378. Available at: \url{https://www.biorxiv.org/content/early/2017/10/11/200378}
#'@export
#'@examples
#'
#' RadialMR(Data$beta.exposure,Data$beta.outcome,Data$se.exposure,Data$se.outcome,Data$SNP,"BOTH","YES","NO",0.05,"NO")$plot
#' 

RadialMR<-function(BXG,BYG,seBXG,seBYG,RSID,METHOD,SCALE,OUTLIERS,ALPHA,SCALEMATCH){
  
  library(ggplot2)

  #Assign temporary SNP IDs if not provided
  suppressWarnings(if(RSID=="NO"){
    RSID<-seq(from=1,to=length(BYG),by=1)
  })
  
  #Define function for plotting circles
  circleFun <- function(center = c(0,0),R,npoints = length(BYG),START,END){
    r = R
    tt <- seq(START,END,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  
  #Define radial IVW procedure
  if(METHOD=="IVW"){
    #Define ratio estimates
    Ratios<-BYG/BXG
    
    #Define inverse variance weights
    W<-((BXG^2)/(seBYG^2))
    
    #Define vector of squared weights
    Wj<-sqrt(W)
    
    #Define vector of weights * ratio estimates
    BetaWj<-Ratios*Wj
    
    #Define IVW Model
    IVW.Model<-lm(BetaWj~-1+Wj)
    
    #Define output of IVW.Model fit
    EstimatesIVW<-summary(lm(BetaWj~-1+Wj))
    
    #Define Egger summary as NULL
    EstimatesEGGER<-NULL
    
    #Define dataframe with weights, ratios and betaweights
    Dat<-data.frame(Wj,BetaWj,Ratios,RSID)
    
    #Define vector of fitted values using IVW.Model
    Dat$Pred<-predict(IVW.Model)
    
    #Define slope parameter for IVW.Model
    IVW.Slope<-EstimatesIVW$coefficients[1]
    
    #Define standard error for slope parameter of IVW.Model
    IVW.SE<-EstimatesIVW$coefficients[2]
    
    #Define confidence interval for IVW.Model slope parameter
    IVW_CI<-confint(IVW.Model)
    
    #Define Q statistic for each individual variant
    Qj<-W*(Ratios-IVW.Slope)^2
    
    #Add Qj values into dataframe
    Dat$Qj<-Qj
    
    #Define total Q statistic
    Total_Q<-sum(Qj)
    
    #Perform chi square test for overall Q statistic
    Total_Q_chi<-1-pchisq(Total_Q,length(BXG-1))
    
    #Generate a dataframe giving overall Q statistic and chi square value for model fit
    QDisplay<-data.frame(Total_Q,Total_Q_chi)
    
    #Redefine column names of QDisplay dataframme
    colnames(QDisplay) = c("Overall_Q","Chi_Sq")
    
    #Define a placeholder vector of 0 values for chi square tests
    Qj_Chi<-rep(0,length(BXG))
    
    #Perform chi square tests for each Qj value
    for(i in 1:length(Qj)){
      Qj_Chi[i]<-1-pchisq(Qj[i],1)
    }
    
    #Add chi squared values to Dat dataframe
    Dat$Qj_Chi<-Qj_Chi
    
    #Define a placeholder vector of 0 values fo identifying outliers
    Out_Indicator<-rep(0,length(BXG))
    
    #For loop defining outlier indicator variable.
    for( i in 1:length(BXG)){
      if(Qj_Chi[i]<ALPHA){
        Out_Indicator[i]<-1
      }
    }
    
    #Include the outlier identifier cariable in the dataframe and define as a factor
    Dat$Outliers<-factor(Out_Indicator)
    levels(Dat$Outliers)[levels(Dat$Outliers)=="0"] <- "Variant"
    levels(Dat$Outliers)[levels(Dat$Outliers)=="1"] <- "Outlier"
    
    #If no outliers are present indicate this is the case
    if(sum(Out_Indicator==0)){
      Q_DATA<-"No significant outliers"
    }
    
    #If outliers are present produce dataframe containing individual Q statistics and chi square values for each outlier
    if(sum(Out_Indicator>0)){
      #Generate a subset of data containing only outliers
      Out_Dat<-subset(Dat, Outliers == "Outlier")
      
      #Construct a datafrae containing SNP IDs, Q statistics and chi-square values for each outlying variant
      Q_DATA<-data.frame(Out_Dat$RSID,Out_Dat$Qj,Out_Dat$Qj_Chi)
      
      #Redefine column names
      colnames(Q_DATA) = c("SNP","Q_statistic","Q_Chisq")
      
      #Order Q_DATA dataframe by smallest chi square value
      Q_DATA<-Q_DATA[with(Q_DATA, order(Q_Chisq)), ]
    }
    
    #Define Radial IVW function
    R.IVW<-function(BXG,BYG,seBXG,seBYG,FULL,OUT_ONLY){
      
      #Define radius for scale and IVW circles
      maximal<-atan(max(abs(Dat$BetaWj))/max(Dat$Wj))
      R.All<-max(abs(Dat$BetaWj))/sin(maximal)
      R.IVW<-R.All+(R.All/12)
      
      #Define scaling parameter
      Label.Scaling<-R.IVW/24
      
      Scalemin<-min(c((max(abs(Ratios))*-1),confint(IVW.Model)[1]))
      
      Scalemax<-max(c(max(abs(Ratios)),confint(IVW.Model)[2]))
      
      #Plot circle for overall scale
      circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),atan(Scalemin),atan(Scalemax))
      cxAll<-circledat_ALLEST$x
      cyAll<-circledat_ALLEST$y
      
      #Plot circle for IVW estimate
      circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
      cxIVW<-circledat_IVWEST$x
      cyIVW<-circledat_IVWEST$y
      
      #Define minimum and maximum values of y axis for scaling
      Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
      Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
      
      Y_Range<-Scalemax-Scalemin
      
      #If only outliers are to be shown
      if(OUT_ONLY=="YES"){
        
        #If the full scale is to be shown
        if(FULL=="YES"){
          
          #If axes scales should be fixed
          if(SCALEMATCH=="YES"){
            
            #If no outliers are present
            if(sum(Out_Indicator)==0){
              
              #Produce plot showing full scale and all variants
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
              print("No significant Outliers")
            }
            
            
            #If outliers are present
            if(sum(Out_Indicator)>0){
              
              #Create a subset of data containing only outliers
              Out_Dat<-subset(Dat, Outliers == "Outlier")
              
              #Redefine scale and IVW circle radii
              maximal<-atan(max(abs(Out_Dat$BetaWj))/max(abs(Out_Dat$Wj)))
              R.All<-max(abs(Out_Dat$BetaWj))/sin(maximal)
              R.IVW<-R.All+(R.All/12)
              
              #Plot circle for overall scale
              circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate plot showing only outliers using the full scale
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                coord_fixed()+theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
              Theta<-rep(0,length(Out_Dat$Outliers))
              a<-rep(0,length(Out_Dat$Outliers))
              b<-rep(0,length(Out_Dat$Outliers))
              for(i in 1:length(Out_Dat$Outliers)){
                Theta[i]<-atan(Out_Dat$Ratios[i])
                b[i]<-sin(Theta[i])*R.All
                a[i]<-cos(Theta[i])*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$Pred[i],linetype="solid",colour="#56B4E9")
                
              }
            }
            
          }
          
          #If axes scales should not be fixed
          if(SCALEMATCH=="NO"){
            
            #If no outliers are present
            if(sum(Out_Indicator)==0){
              
              #Produce plot showing full scale and all variants
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
              print("No significant Outliers")
            }
            
            
            #If outliers are present
            if(sum(Out_Indicator)>0){
              
              #Create a subset of data containing only outliers
              Out_Dat<-subset(Dat, Outliers == "Outlier")
              
              #Redefine scale and IVW circle radii
              maximal<-atan(max(abs(Out_Dat$BetaWj))/max(abs(Out_Dat$Wj)))
              R.All<-max(abs(Out_Dat$BetaWj))/sin(maximal)
              R.IVW<-R.All+(R.All/12)
              
              #Plot circle for overall scale
              circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate plot showing only outliers using the full scale
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              #Draw individual ratio estimate lines for each outlier and Q distance indicator lines
              Theta<-rep(0,length(Out_Dat$Outliers))
              a<-rep(0,length(Out_Dat$Outliers))
              b<-rep(0,length(Out_Dat$Outliers))
              for(i in 1:length(Out_Dat$Outliers)){
                Theta[i]<-atan(Out_Dat$Ratios[i])
                b[i]<-sin(Theta[i])*R.All
                a[i]<-cos(Theta[i])*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$Pred[i],linetype="solid",colour="#56B4E9")
                
              }
            }
          }
        }
        
        #If the full scale is not to be shown
        if(FULL=="NO") {
          
          #If no outliers are present
          if(sum(Out_Indicator)==0){
            
            if(min(BetaWj) < (sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)){
              Y_MINIVW<-min(BetaWj)
            } 
            
            if(min(BetaWj) >(sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)){
              Y_MINIVW<-(sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)
            } 
            
            if(max(BetaWj) >(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)){
              Y_MAXIVW<-max(BetaWj)
            } 
            
            if(max(BetaWj) <(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)){
              Y_MAXIVW<-(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)
            } 
            
            Y_Range<-Y_MAXIVW-Y_MINIVW
            
            #Generate plot showing all variants without using the full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="IVW Radial") + geom_point(aes(colour="Variant"))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1]))*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1]))*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2]))*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2]))*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,cos(atan(IVW.Slope))*(R.IVW+(Label.Scaling*3))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
              scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
            
            #Show "No significant Outliers" to explain why the full set of variants is displayed
            print("No significant Outliers")
            
          }
          
          #If outliers are present
          if(sum(Out_Indicator)>0){
            
            #Create a subset of data containing only outliers
            Out_Dat<-subset(Dat, Outliers == "Outlier")
            
            #Redefine IVW and scale circle radii
            maximal<-atan(max(abs(Out_Dat$BetaWj))/max(abs(Out_Dat$Wj)))
            R.All<-max(abs(Out_Dat$BetaWj))/sin(maximal)
            R.IVW<-R.All+(R.All/12)
            
            circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),-pi/2,pi/2)
            cxAll<-circledat_ALLEST$x
            cyAll<-circledat_ALLEST$y
            
            circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
            cxIVW<-circledat_IVWEST$x
            cyIVW<-circledat_IVWEST$y
            
            if(min(Out_Dat$BetaWj) < (sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)){
              Y_MINIVW<-min(Out_Dat$BetaWj)
            } 
            
            if(min(Out_Dat$BetaWj) >(sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)){
              Y_MINIVW<-(sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)
            } 
            
            if(max(Out_Dat$BetaWj) >(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)){
              Y_MAXIVW<-max(Out_Dat$BetaWj)
            } 
            
            if(max(Out_Dat$BetaWj) <(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)){
              Y_MAXIVW<-(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)
            } 
            
            Y_Range<-Y_MAXIVW-Y_MINIVW
            
            #Generate plot showing only outliers and not using the full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1]))*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1]))*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2]))*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2]))*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,cos(atan(IVW.Slope))*(R.IVW+(Label.Scaling*3))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
              scale_colour_manual(breaks=c("Outlier","IVW"),values=c("#56B4E9","#E69F00","white"))
            
            #Draw lines indicating individual Q contributions and give text displaying each value
            for(i in 1:length(Out_Dat$Outliers)){
              B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$Pred[i],linetype="solid",colour="#56B4E9")
              B<- B + geom_text(x=Out_Dat$Wj[i],y=(Out_Dat$BetaWj[i]*0.9),label=round(Out_Dat$Qj[i],digits=2),size=2.5)
            }
          }
        }
      }
      
      #If all variants are to be shown
      if(OUT_ONLY=="NO"){
        
        #If the full scale is to be used
        if(FULL=="YES"){
          
          #If axes scales should not be fixed
          if(SCALEMATCH=="NO"){
            if(sum(Out_Indicator)==0){
              
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour="Variant"))+
                geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                
              }
              
              #Indicate that no outliers are present
              print("No significant Outliers")
            }
            
            #If outliers are present
            if(sum(Out_Indicator)>0){
              
              #Create a subset consiting of only outlying variants
              Out_Dat<-subset(Dat, Outliers == "Outlier")
              
              #Plot circle for overall scale
              circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              
              circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing all variants, and IVW estimate using the full scale
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
                geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
                theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#56B4E9","#E69F00","#000000"))
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Draw lines showing the ratio estimate for each individual variant
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
            }
            
            
            
            
          }
          
          ##If axes scales should be fixed
          if(SCALEMATCH=="YES"){
          #If no outliers are present
          if(sum(Out_Indicator)==0){
            
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="IVW Radial") + geom_point(aes(colour="Variant"))+
              geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
              coord_fixed()+scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
            
            if(max(abs(Ratios))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
              
            }
            
            if(max(abs(Ratios))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
              
            }
            for(i in 1:length(BYG)){
              Theta<-atan(Ratios)
              
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              
            }
            
            #Indicate that no outliers are present
            print("No significant Outliers")
          }
          
          #If outliers are present
          if(sum(Out_Indicator)>0){
            
            #Create a subset consiting of only outlying variants
            Out_Dat<-subset(Dat, Outliers == "Outlier")
            
            #Plot circle for overall scale
            circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),atan(Scalemin),atan(Scalemax))
            cxAll<-circledat_ALLEST$x
            cyAll<-circledat_ALLEST$y
            
            circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
            cxIVW<-circledat_IVWEST$x
            cyIVW<-circledat_IVWEST$y
            
            #Define minimum and maximum values of y axis for scaling
            Y_MINIVW<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
            Y_MAXIVW<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
            
            Y_Range<-Scalemax-Scalemin
            
            #Generate a plot showing all variants, and IVW estimate using the full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,R.IVW+(Label.Scaling*5)),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
              coord_fixed()+theme(legend.title=element_blank())+scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#56B4E9","#E69F00","#000000"))
            
            if(max(abs(Ratios))>2){
              Y_Scale<-Y_Range/6
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
              
            }
            
            if(max(abs(Ratios))<2){
              
              Y_Scale<-Y_Range/4
              Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
              B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
              B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
              
            }
            
            #Draw lines showing the ratio estimate for each individual variant
            for(i in 1:length(BYG)){
              Theta<-atan(Ratios)
              b<-sin(Theta)*R.All
              a<-cos(Theta)*R.All
              B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
            }
          }
          }
        }
        
        #If the full scale is not to be used
        if(FULL=="NO") {
          
          #If no outliers are present
          if(sum(Out_Indicator)==0){
            
            if(min(BetaWj) < (sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)){
              Y_MINIVW<-min(BetaWj)
            } 
            
            if(min(BetaWj) >(sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)){
              Y_MINIVW<-(sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)
            } 
            
            if(max(BetaWj) >(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)){
              Y_MAXIVW<-max(BetaWj)
            } 
            
            if(max(BetaWj) <(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)){
              Y_MAXIVW<-(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)
            } 
            
            Y_Range<-Y_MAXIVW-Y_MINIVW
            
            
            #Generate plot showing all variants and the IVW estimate not using the full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="IVW Radial") + geom_point(aes(colour="Variant"))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1]))*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1]))*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2]))*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2]))*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,cos(atan(IVW.Slope))*(R.IVW+(Label.Scaling*3))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
              scale_color_manual(name="Estimate",breaks=c("Variant","IVW"),values=c("IVW"="#56B4E9","Variant"="#000000"))
            
            #Provide indicator that no outliers are present
            print("No significant Outliers")
            
            
            
          }
          
          #If outliers are present
          if(sum(Out_Indicator)>0){
            
            #Create a subset of data containing only outlying variants
            Out_Dat<-subset(Dat, Outliers == "Outlier")
            
            #Redefine IVW and scale circle radii
            circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),-pi/2,pi/2)
            cxAll<-circledat_ALLEST$x
            cyAll<-circledat_ALLEST$y
            
            circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
            cxIVW<-circledat_IVWEST$x
            cyIVW<-circledat_IVWEST$y
            
            if(min(BetaWj) < (sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)){
              Y_MINIVW<-min(BetaWj)
            } 
            
            if(min(BetaWj) >(sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)){
              Y_MINIVW<-(sin(atan(confint(IVW.Model)[1]))*R.IVW+Label.Scaling)
            } 
            
            if(max(BetaWj) >(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)){
              Y_MAXIVW<-max(BetaWj)
            } 
            
            if(max(BetaWj) <(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)){
              Y_MAXIVW<-(sin(atan(confint(IVW.Model)[2]))*R.IVW+Label.Scaling)
            } 
            
            Y_Range<-Y_MAXIVW-Y_MINIVW
            
            #Generate a plot showing all variants and the IVW estimate, not using the full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="IVW Radial") + geom_point(aes(colour = Outliers))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1]))*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1]))*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2]))*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2]))*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,cos(atan(IVW.Slope))*(R.IVW+(Label.Scaling*3))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINIVW,Y_MAXIVW))+
              scale_colour_manual(breaks=c("Variant","Outlier","IVW"),values=c("#56B4E9","#E69F00","#000000"))
          }
        }
      }
      
      #Return the respective plot "B"
      return(B)
    }
    
    #Run the R.IVW function and save the plot output as "A"
    A<-R.IVW(BXG,BYG,seBXG,seBYG,SCALE,OUTLIERS)
  }
  
  #Define radial MR-Egger procedure
  if(METHOD=="EGGER"){
    
    #Define ratio estimates using MR-Egger beta coefficients
    Ratios<-BYG /BXG
    
    #Define a vector of weights W
    W<-((BXG ^2)/(seBYG^2))
    
    #Define a vector of root weights Wj for each variant
    Wj<-sqrt(W)
    
    #Define a vector of root weights multiplied by the ratio estimate for each variant
    BetaWj<-Ratios*Wj
    
    #Fit radial MR-Egger model
    Egger.Model<-lm(BetaWj~Wj)
    
    #Define an object containing summary estimates for Egger.Model
    EstimatesEGGER<-summary(Egger.Model)
    
    #Define indicator that IVW estimates are not present
    EstimatesIVW<-NULL
    
    # Construct a data frame containing weights, beta weights and SNP IDs
    Dat<-data.frame(Wj,BetaWj,RSID)
    
    #Define a vector of fitted values using Egger.Model and input data. 
    Dat$Pred<-predict(Egger.Model)
    
    #Define intercept of Egger.Model
    Egger.Intercept<-EstimatesEGGER$coefficients[1]
    
    #Define slope of Egger.Model
    Egger.Slope<-EstimatesEGGER$coefficients[2]
    
    #Define standard error for intercept of Egger.Model
    Eggerintercept.SE<-EstimatesEGGER$coefficients[3]
    
    #Define standard error for slope of Egger.Model
    Eggerslope.SE<-EstimatesEGGER$coefficients[4]
    
    #Define 95% confidence interval for Egger.Model intercept
    Eggerint_CI<-as.numeric(confint(Egger.Model)[1,])
    
    #Define 95% confidence interval for Egger.Model slope
    Eggerslope_CI<-as.numeric(confint(Egger.Model)[2,])
    
    #Define Q statistic for each individual variant
    Qj<-W*(Ratios-(Egger.Intercept/Wj)-Egger.Slope)^2
    
    #Add Qj values into dataframe
    Dat$Qj<-Qj
    
    #Define total Q statistic
    Total_Q<-sum(Qj)
    
    #Perform chi square test for overall Q statistic
    Total_Q_chi<-1-pchisq(Total_Q,length(BXG-1))
    
    #Define a dataframe to show overall Q statistic and respective chi-square value
    QDisplay<-data.frame(Total_Q,Total_Q_chi)
    
    #Rename columns of QDisplay dataframe
    colnames(QDisplay) = c("Overall_Q","Chi_Sq")
    
    #Define a placeholder vector of 0 values for chi square tests
    Qj_Chi<-rep(0,length(BXG))
    
    #Perform chi square tests for each Qj value
    for(i in 1:length(Qj)){
      Qj_Chi[i]<-1-pchisq(Qj[i],1)
    }
    
    #Add individual Q chi square results for each variant to Dat data frame
    Dat$Qj_Chi<-Qj_Chi
    
    #Define a placeholder vector of 0 values fo identifying outliers
    Out_Indicator<-rep(0,length(BXG))
    
    #For loop defining outlier indicator variable.
    for( i in 1:length(BXG)){
      if(Qj_Chi[i]<ALPHA){
        
        Out_Indicator[i]<-1
        
      }
    }
    
    #Include the outlier identifier cariable in the dataframe and define as a factor
    Dat$Outliers<-factor(Out_Indicator)
    levels(Dat$Outliers)[levels(Dat$Outliers)=="0"] <- "Variant"
    levels(Dat$Outliers)[levels(Dat$Outliers)=="1"] <- "Outlier"
    
    #If no outliers are present
    if(sum(Out_Indicator==0)){
      #Provide an indicator no outliers are present
      Q_DATA<-"No significant outliers"
    }
    
    #If outliers are present
    if(sum(Out_Indicator>0)){
      
      #Define a subset containing only outlying variants
      Out_Dat<-subset(Dat, Outliers == "Outlier")
      
      #Create a dataframe containing individual Q statistics, chi-square values and SNP IDs for each outlier
      Q_DATA<-data.frame(Out_Dat$RSID,Out_Dat$Qj,Out_Dat$Qj_Chi)
      
      #Rename the columns of Q_DATA
      colnames(Q_DATA) = c("SNP","Q_statistic","Q_Chisq")
      
      #Order the Q_DATA dataset by smallest chi square value
      Q_DATA<-Q_DATA[with(Q_DATA, order(Q_Chisq)), ]
    }
    
    #Define R.EGGER function
    R.EGGER<-function(BXG,BYG,seBXG,seBYG,OUT_ONLY){
      
      #Define MR-Egger circle radius
      maximal<-atan(max(abs(Dat$BetaWj))/max(abs(Dat$Wj)))
      R.Egger<-max(abs(Dat$BetaWj))/sin(maximal)
      
      #Define scaling parameter for labels
      Label.Scaling_Egg<-R.Egger/24
      
      #Define minimum and maximum range of values for y axis
      Y_MINEGGER<-min(  c(min(BetaWj),(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg))+Egger.Intercept  ))
      Y_MAXEGGER<-max(  c(max(BetaWj),(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg))+Egger.Intercept  ))
      
      #Define values for MR-Egger circle estimate
      circledat_EGGEREST <- circleFun(c(0, Egger.Intercept),R.Egger,npoints = length(BYG),atan(Eggerslope_CI[1]),atan(Eggerslope_CI[2]))
      cxEgger<-circledat_EGGEREST$x
      cyEgger<-circledat_EGGEREST$y
      
      #If all variants are to be shown
      if(OUT_ONLY=="NO"){
        
        #If no outliers are present
        if(sum(Out_Indicator)==0){
          
          #Generates a plot showing all variants not using the full scale
          B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
          )+labs(title="EGGER Radial") + geom_point(aes(colour="variant"))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling_Egg), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling_Egg) +
                        Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg))+
                                                                                              Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                                                                                         Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+
            scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling_Egg*3))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINEGGER,Y_MAXEGGER))+scale_color_manual(breaks=c("variant","MR-Egger"),values=c("MR-Egger"="#D55E00","variant"="#000000"))+
            theme(legend.title=element_blank())
          
          #Indicates that no outliers are present
          print("No significant Outliers")
        }
        
        #If outliers are present
        if(sum(Out_Indicator)>0){
          
          #Generate a plot showing all variants not using full scale
          B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
          )+labs(title="EGGER Radial") + geom_point(aes(colour=Outliers))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling_Egg), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling_Egg)+
                        Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg))+
                                                                                              Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                                                                                         Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+
            scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling_Egg*3))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINEGGER,Y_MAXEGGER))+scale_colour_manual(breaks=c("Variant","Outlier","MR-Egger"),values=c("#D55E00","#E69F00","#000000"))+
            theme(legend.title=element_blank())
        }
      }
      
      #If only outliers are to be shown
      if(OUT_ONLY=="YES"){
        
        #If no outliers are present
        if(sum(Out_Indicator)==0){
          
          #Generates a plot showing all variants not using the full scale
          B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
          )+labs(title="EGGER Radial") + geom_point(aes(colour="variant"))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling_Egg), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling_Egg) +
                        Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg))+
                                                                                              Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                                                                                         Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+
            scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling_Egg*3))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINEGGER,Y_MAXEGGER))+scale_color_manual(breaks=c("variant","MR-Egger"),values=c("MR-Egger"="#D55E00","variant"="#000000"))+
            theme(legend.title=element_blank())
          
          #Indicates that no outliers are present
          print("No significant Outliers")
        }
        
        #If outliers are present
        if(sum(Out_Indicator)>0){
          
          #Create a subset containing only outliers
          Out_Dat<-subset(Dat, Outliers == "Outlier")
          
          #Redefine MR-Egger circle radius
          maximal<-atan(max(abs(Out_Dat$BetaWj))/max(abs(Out_Dat$Wj)))
          R.Egger<-max(abs(Out_Dat$BetaWj))/sin(maximal)
          
          #Define scaling parameter
          Label.Scaling_Egg<-R.Egger/24
          
          circledat_EGGEREST <- circleFun(c(0, Egger.Intercept),R.Egger,npoints = length(BYG),atan(Eggerslope_CI[1]),atan(Eggerslope_CI[2]))
          cxEgger<-circledat_EGGEREST$x
          cyEgger<-circledat_EGGEREST$y
          
          #Generate a plot showing only outlying variants and MR-Egger estimate without full scale
          B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
          )+labs(title="EGGER Radial") + geom_point(aes(colour=Outliers))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
            geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling_Egg), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling_Egg) +
                        Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling_Egg))+
                                                                                              Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling_Egg))+
                                                                                                                                                                         Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+theme_bw() +
            theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
            geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+
            scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling_Egg*3))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MINEGGER,Y_MAXEGGER))+scale_colour_manual(breaks=c("Outlier","MR-Egger"),values=c("#D55E00","#E69F00","white"))+
            theme(legend.title=element_blank())
          
          #Draw lines indicating Q statistic for each outlier, and present a label for each statistic
          for(i in 1:length(Out_Dat$Outliers)){
            B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$Pred[i],linetype="solid",colour="#D55E00")
            B<- B + geom_text(x=Out_Dat$Wj[i],y=(Out_Dat$BetaWj[i]*0.9),label=round(Out_Dat$Qj[i],digits=2),size=2.5)
          }
        }
      }
      
      #Return the specified plot "B"
      return(B)
    }
    
    #Define the plot output from the R.EGGER function as "A"
    A<-R.EGGER(BXG,BYG,seBXG,seBYG,OUTLIERS)
  }
  
  #Define procedure for implementing IVW and MR-Egger simulatenously
  if(METHOD=="BOTH"){
    
    #Define a vector of ratio estimates for each variant
    Ratios<-BYG /BXG
    
    #Define a vector of weights for each variant
    W<-((BXG ^2)/(seBYG^2))
    
    #Define a vector of root weights for each variant
    Wj<-sqrt(W)
    
    #Define a vector of Ratio multiplied by root weights for each variant
    BetaWj<-Ratios*Wj
    
    #Fit an MR-Egger model
    Egger.Model<-lm(BetaWj~Wj)
    
    #Define an object containing summary estimates from Egger.Model
    EstimatesEGGER<-summary(Egger.Model)
    
    #Construct a data frame containing root weights, betaweights, ratio estimates and SNP IDs
    Dat<-data.frame(Wj,BetaWj,Ratios,RSID)
    
    #Defines Egger.Model intercept
    Egger.Intercept<-EstimatesEGGER$coefficients[1]
    
    #Defines Egger.Model slope
    Egger.Slope<-EstimatesEGGER$coefficients[2]
    
    #Defines standard error for Egger.Model intercept
    Eggerintercept.SE<-EstimatesEGGER$coefficients[3]
    
    #Defines standard error for Egger.Model slope
    Eggerslope.SE<-EstimatesEGGER$coefficients[4]
    
    #Define 95% confidence interval for Egger.Model intercept
    Eggerint_CI<-as.numeric(confint(Egger.Model)[1,])
    
    #Define 95% confidence interval for Egger.Model slope
    Eggerslope_CI<-as.numeric(confint(Egger.Model)[2,])
    
    #Fit IVW model
    IVW.Model<-lm(BetaWj~-1+Wj)
    
    #Define summary estimates from IVW.Model
    EstimatesIVW<-summary(lm(BetaWj~-1+Wj))
    
    #Define slope estimate for IVW.Model
    IVW.Slope<-EstimatesIVW$coefficients[1]
    
    #Define standard error for IVW slope estimate
    IVW.SE<-EstimatesIVW$coefficients[2]
    
    #Define 95% confidence interval for IVW.Model slope estimate
    IVW_CI<-confint(IVW.Model)
    
    #Define IVW Q statistic for each individual variant
    QjIVW<-W*(Ratios-IVW.Slope)^2
    
    #Add IVW Qj values into dataframe
    Dat$QjIVW<-QjIVW
    
    #Define total IVW Q statistic
    Total_QIVW<-sum(QjIVW)
    
    #Perform chi square test for overall Q IVW statistic
    Total_Q_chiIVW<-1-pchisq(Total_QIVW,length(BXG-1))
    
    #Define a placeholder vector of 0 values for Qj IVW chi square tests
    Qj_ChiIVW<-rep(0,length(BXG))
    
    #Perform chi square tests for each  Qj IVW value
    for(i in 1:length(QjIVW)){
      Qj_ChiIVW[i]<-1-pchisq(QjIVW[i],1)
    }
    
    #Add individual chi square tests to Dat dataframe (IVW)
    Dat$Qj_ChiIVW<-Qj_ChiIVW
    
    #Define a placeholder vector of 0 values fo identifying IVW outliers
    Out_IndicatorIVW<-rep(0,length(BXG))
    
    #For loop defining IVW outlier indicator variable.
    for( i in 1:length(BXG)){
      if(Qj_ChiIVW[i]<ALPHA){
        Out_IndicatorIVW[i]<-1
      }
    }
    
    #Include the IVW outlier identifier variable in the dataframe and define as a factor
    Dat$OutliersIVW<-factor(Out_IndicatorIVW)
    levels(Dat$OutliersIVW)[levels(Dat$OutliersIVW)=="0"] <- "Variant"
    levels(Dat$OutliersIVW)[levels(Dat$OutliersIVW)=="1"] <- "Outlier"
    
    #Define MR-Egger Q statistic for each individual variant
    QjEGGER<-W*(Ratios-(Egger.Intercept/Wj)-Egger.Slope)^2
    
    #Add MR-Egger Qj values into dataframe
    Dat$QjEGGER<-QjEGGER
    
    #Define total MR-Egger Q statistic
    Total_QEGGER<-sum(QjEGGER)
    
    #Perform chi square test for overall MR-Egger Q statistic
    Total_Q_chiEGGER<-1-pchisq(Total_QEGGER,length(BXG-1))
    
    #Define a placeholder vector of 0 values for MR-Egger Q chi square tests
    Qj_ChiEGGER<-rep(0,length(BXG))
    
    #Perform chi square tests for each Qj value
    for(i in 1:length(QjEGGER)){
      Qj_ChiEGGER[i]<-1-pchisq(QjEGGER[i],1)
    }
    
    #Add individual MR-Egger Q chi square values to Dat dataframe
    Dat$Qj_ChiEGGER<-Qj_ChiEGGER
    
    #Define a placeholder vector of 0 values fo identifying MR-Egger outliers
    Out_IndicatorEGGER<-rep(0,length(BXG))
    
    #For loop defining MR-Egger outlier indicator variable.
    for( i in 1:length(BXG)){
      if(Qj_ChiEGGER[i]<ALPHA){
        Out_IndicatorEGGER[i]<-1
      }
    }
    
    #Include the MR-Egger outlier identifier variable in the dataframe and define as a factor
    Dat$OutliersEGGER<-factor(Out_IndicatorEGGER)
    levels(Dat$OutliersEGGER)[levels(Dat$OutliersEGGER)=="0"] <- "Variant"
    levels(Dat$OutliersEGGER)[levels(Dat$OutliersEGGER)=="1"] <- "Outlier"
    
    ###Combined###
    
    #Create a dataframe showing overall Q for each model and respective chi square value
    QDisplay<-data.frame(Total_QIVW,Total_Q_chiIVW,Total_QEGGER,Total_Q_chiEGGER)
    
    #Redefine column names
    colnames(QDisplay) = c("Q_IVW","Chi_Sq_IVW","Q_EGGER","Chi_Sq_EGGER")
    
    #If there are no outliers for either model
    if(sum(Out_IndicatorIVW==0) & sum(Out_IndicatorEGGER==0)){
      
      #Indicate there are no outliers detected
      Q_DATA<-"No significant outliers" 
      
    }
    
    #If there are IVW outliers but not MR-Egger outliers
    if(sum(Out_IndicatorIVW>0) & sum(Out_IndicatorEGGER==0)){
      
      #Define a subset of outlying variants
      Out_Dat<-subset(Dat, OutliersIVW == "Outlier")
      
      #Generate a dataframe containing SNP IDs, Q statistics and chi square values for IVW outliers
      Q_DATA<-data.frame(Out_Dat$RSID,Out_Dat$QjIVW,Out_Dat$Qj_ChiIVW)
      
      #Redefine column names for Q_DATA
      colnames(Q_DATA) = c("SNP","Q_statisticIVW","Q_ChisqIVW")
      
      #Order Q_Data by smallest chi-square value
      Q_DATA<-Q_DATA[with (Q_DATA, order(Q_ChisqIVW)), ]
    }
    
    #Order Q_Data by smallest chi-square value
    if(sum(Out_IndicatorIVW==0) & sum(Out_IndicatorEGGER>0)){
      
      #Define a subset of outlying variants
      Out_Dat<-subset(Dat, OutliersEGGER == "Outlier")
      
      #Generate a dataframe containing SNP IDs, Q statistics and chi square values for MR-Egger outliers
      Q_DATA<-data.frame(Out_Dat$RSID,Out_Dat$QjEGGER,Out_Dat$Qj_ChiEGGER)
      
      #Redefine column names for Q_DATA
      colnames(Q_DATA) = c("SNP","Q_statisticEGGER","Q_ChisqEGGER")
      
      #Order Q_Data by smallest chi-square value
      Q_DATA<-Q_DATA[with (Q_DATA, order(Q_ChisqEGGER)), ]
    }
    
    #If there are IVW outliers and MR-Egger outliers
    if(sum(Out_IndicatorIVW>0) & sum(Out_IndicatorEGGER>0)){
      
      #Define a subset of outlying variants
      Out_Dat<-subset(Dat, OutliersIVW == "Outlier"|OutliersEGGER== "Outlier")
      
      #Generate a dataframe containing SNP IDs, Q statistics and chi square values for all outliers
      Q_DATA<-data.frame(Out_Dat$RSID,Out_Dat$QjIVW,Out_Dat$Qj_ChiIVW,Out_Dat$QjEGGER,Out_Dat$Qj_ChiEGGER)
      
      #Redefine column names for Q_DATA
      colnames(Q_DATA) = c("SNP","Q_statisticIVW","Q_ChisqIVW","Q_statisticEgger","Q_ChisqEgger")
      
      #Order Q_Data by smallest chi-square value
      Q_DATA<-Q_DATA[with (Q_DATA, order(Q_ChisqIVW,Q_ChisqEgger)), ]
      
    }
    
    #Define a function for plotting both IVW and MR-Egger estimates simulateously
    BOTH.Function<-function(BXG,BYG,seBXG,seBYG,FULL,OUT_ONLY){
      
      #Define scale, IVW and MR-Egger circle radii
      maximal<-atan(max(abs(Dat$BetaWj))/max(abs(Dat$Wj)))
      R.All<-max(abs(Dat$BetaWj))/sin(maximal)
      R.IVW<-R.All+(R.All/12)
      R.Egger<-R.All+(R.All/4)
      
      #Define general scaling parameter
      Label.Scaling<-R.IVW/24
      
      #Define scaling parameter for MR-Egger
      Label.Scaling_Egg<-R.Egger/8
      
      #Define vector of fitted values for Egger.Model
      Dat$PredEgger<-predict(Egger.Model)
      
      #Define vector of fitted values for IVW.Model
      Dat$PredIVW<-predict(IVW.Model)
      
      Scalemin<-min(c(max(abs(Ratios))*-1,confint(IVW.Model)[1],Eggerslope_CI[1]))
      
      Scalemax<-max(c(max(abs(Ratios)),confint(IVW.Model)[2],Eggerslope_CI[2]))
      
      #Define minimum and maximum values of y axis for scaling
      Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
      Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
      
      Y_Range<-Scalemax-Scalemin
      
      #Define circle plotting points for scale, IVW and MR-Egger estimates
      circledat_EGGEREST <- circleFun(c(0, Egger.Intercept),R.Egger,npoints = length(BYG),atan(Eggerslope_CI[1]),atan(Eggerslope_CI[2]))
      cxEgger<-circledat_EGGEREST$x
      cyEgger<-circledat_EGGEREST$y
      
      #Plot circle for overall scale
      circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),atan(Scalemin),atan(Scalemax))
      cxAll<-circledat_ALLEST$x
      cyAll<-circledat_ALLEST$y
      
      circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
      cxIVW<-circledat_IVWEST$x
      cyIVW<-circledat_IVWEST$y
      
      #Define an outlier indicator for both IVW and MR-Egger outliers
      Out_Combine<-Out_IndicatorIVW+Out_IndicatorEGGER
      
      #Redefine outliers in both models as 1 to define a binary indicator
      for(i in 1:length(BXG)){
        if(Out_Combine[i]==2){
          Out_Combine[i]<-1
        }
      }
      
      out_types<-rep(0,length(BYG))
      
      for(i in 1:length(out_types)){
        
        if(Out_IndicatorIVW[i]==0 & Out_IndicatorEGGER[i]==0){
          out_types[i]<-"Variant"
          
        }
        
        if(Out_IndicatorIVW[i]>0 & Out_IndicatorEGGER[i]>0){
          out_types[i]<-"IVW and MR-Egger Outlier"
          
        }
        
        if(Out_IndicatorIVW[i]==0 & Out_IndicatorEGGER[i]>0){
          out_types[i]<-"MR-Egger Outlier"
          
        }
        
        if(Out_IndicatorIVW[i]>0 & Out_IndicatorEGGER[i]==0){
          out_types[i]<-"IVW Outlier"
        }
        
      }
      
      Dat$out_types<-factor(out_types)
      
      #Include outlier indicator as a factor variable in Dat dataframe
      Dat$Outliers<-factor(Out_Combine)
      levels(Dat$Outliers)[levels(Dat$Outliers)=="0"] <- "Variant"
      levels(Dat$Outliers)[levels(Dat$Outliers)=="1"] <- "Outlier"
      
      #If all variants are to be shown
      if(OUT_ONLY=="NO"){
        
        #If the full scale is to be used
        if(FULL=="YES"){
          
          #If axes scales should not be fixed
          if(SCALEMATCH=="NO"){
            #If no outliers are present
            if(sum(Out_Combine)==0){
              
              #Define minimum and maximum values of y axis for scaling
              Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing both IVW and MR-Egger estimates, using the full scale and all variants
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant"))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*5))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
                geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
                geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
                geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
                geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+
                                   Egger.Intercept,colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+ theme(legend.title=element_blank())                   
              
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Add lines showing ratio estimates for each individual variant
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
              
              #Indicate that no outliers are detected
              print("No significant Outliers")
            }
            
            #If outliers are present
            if(sum(Out_Combine)>0){
              
              #Define minimum and maximum values of y axis for scaling
              Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot which shows the IVW and MR-Egger estimates, using the full scale and showing all variants
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*5))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
                geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
                geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
                geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())                   
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Draw lines showing ratio estimate for each individual variant
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
            }
          }

          #If axes scales should be fixed
          if(SCALEMATCH=="YES"){
            
            #If no outliers are present
            if(sum(Out_Combine)==0){
              
              #Define minimum and maximum values of y axis for scaling
              Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing both IVW and MR-Egger estimates, using the full scale and all variants
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant"))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*5))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+coord_fixed()+
                geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
                geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
                geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
                geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+
                                   Egger.Intercept,colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+ theme(legend.title=element_blank())                   
              
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Add lines showing ratio estimates for each individual variant
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
              
              #Indicate that no outliers are detected
              print("No significant Outliers")
            }
            
            #If outliers are present
            if(sum(Out_Combine)>0){
              
              #Define minimum and maximum values of y axis for scaling
              Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot which shows the IVW and MR-Egger estimates, using the full scale and showing all variants
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                   axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*5))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+coord_fixed()+
                geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
                geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
                geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())                   
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Draw lines showing ratio estimate for each individual variant
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
            }
          }
        }
        
        #If the full scale is not to be used
        if(FULL=="NO"){
          
          #If outliers are present
          if(sum(Out_Combine)==0){
            
            MINCILIM<-min(sin(atan(confint(IVW.Model)[1]))*R.IVW + Label.Scaling,sin(atan(Eggerslope_CI[1]))*(R.Egger +Label.Scaling_Egg-Egger.Intercept))
            
            MAXCILIM<-max(sin(atan(confint(IVW.Model)[2]))*R.IVW + Label.Scaling,sin(atan(Eggerslope_CI[2]))*(R.Egger +Label.Scaling_Egg-Egger.Intercept))
            
            if(min(BetaWj) < MINCILIM){
              Y_MIN<-min(BetaWj)
            } 
            
            if(min(BetaWj) > MINCILIM){
              Y_MIN<-MINCILIM
            } 
            
            if(max(BetaWj) >MAXCILIM){
              Y_MAX<-max(BetaWj)
            } 
            
            if(max(BetaWj) <MAXCILIM){
              Y_MAX<-MAXCILIM
            } 
            
            Y_Range<-Y_MAX-Y_MIN
            
            #Generate a plot showing all variants and both IVW and MR-Egger estimates, without the full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant")) +geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
              geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*3))),expand=c(0,0))+geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+
                                                                                                                                    Egger.Intercept,colour="MR-Egger"))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+
              theme(legend.title=element_blank())
            
            #Indicate that no outliers are detected
            print("No significant Outliers")
          }
          
          #If outliers are present
          if(sum(Out_Combine)>0){
            
            MINCILIM<-min(sin(atan(confint(IVW.Model)[1]))*R.IVW + Label.Scaling,sin(atan(Eggerslope_CI[1]))*(R.Egger +Label.Scaling_Egg-Egger.Intercept))
            
            MAXCILIM<-max(sin(atan(confint(IVW.Model)[2]))*R.IVW + Label.Scaling,sin(atan(Eggerslope_CI[2]))*(R.Egger +Label.Scaling_Egg-Egger.Intercept))
            
            if(min(BetaWj) < MINCILIM){
              Y_MIN<-min(BetaWj)
            } 
            
            if(min(BetaWj) > MINCILIM){
              Y_MIN<-MINCILIM
            } 
            
            if(max(BetaWj) >MAXCILIM){
              Y_MAX<-max(BetaWj)
            } 
            
            if(max(BetaWj) <MAXCILIM){
              Y_MAX<-MAXCILIM
            } 
            
            Y_Range<-Y_MAX-Y_MIN
            
            #Generate a plot showing all variants and both IVW and MR-Egger estimates, without the full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types)) +geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
              geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*3))),expand=c(0,0))+geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+
                                                                                                                                    Egger.Intercept,colour="MR-Egger"))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+scale_color_manual(breaks=c("Variant","IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())
          }
        }
      }
      
      
      #If only outlying variants are to be shown
      if(OUT_ONLY=="YES"){
        
        #If the full scale is to be used
        if(FULL=="YES"){
          
          #If axes scales should not be fixed
          if(SCALEMATCH=="NO"){
            #If no outliers are present
            if(sum(Out_Combine)==0){
              
              
              #Define minimum and maximum values of y axis for scaling
              Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing both IVW and MR-Egger estimates, using the full scale and all variants
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant"))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*5))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
                geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
                geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
                geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
                geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+
                                   Egger.Intercept,colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+ theme(legend.title=element_blank())                   
              
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Add lines showing ratio estimates for each individual variant
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
              
              #Indicate that no outliers are detected
              print("No significant Outliers")
            }
            
            #If outliers are present
            if(sum(Out_Combine)>0){
              
              #Create a subset containing only outlying variants
              Out_Dat<-subset(Dat, Outliers == "Outlier")
              
              #Redefine IVW, MR-Egger and scale radii
              maximal<-atan(max(abs(Out_Dat$BetaWj))/max(abs(Out_Dat$Wj)))
              R.All<-max(abs(Out_Dat$BetaWj))/sin(maximal)
              R.IVW<-R.All+(R.All/12)
              R.Egger<-R.All+(R.All/4)
              
              #Redefine general label scaling parameter
              Label.Scaling<-R.IVW/24
              
              #Redefine MR-Egger label scaling parameter
              Label.Scaling_Egg<-R.Egger/8
              
              #Redefine values for plotting IVW, MR-Egger and scale circles
              circledat_EGGEREST <- circleFun(c(0, Egger.Intercept),R.Egger,npoints = length(BYG),atan(Eggerslope_CI[1]),atan(Eggerslope_CI[2]))
              cxEgger<-circledat_EGGEREST$x
              cyEgger<-circledat_EGGEREST$y
              
              #Plot circle for overall scale
              circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              
              circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generates a plot showing the MR-Egger and IVW estimates, showing only outliers and using the full scale
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*5))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+
                geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
                geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
                geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
                geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+scale_color_manual(breaks=c("IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="white","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())
              
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Draw individual ratio estimates for each outlying variant, and Q statistic indicators of respective colours      
              Theta<-rep(0,length(Out_Dat$Outliers))
              a<-rep(0,length(Out_Dat$Outliers))
              b<-rep(0,length(Out_Dat$Outliers))
              for(i in 1:length(Out_Dat$Outliers)){
                Theta[i]<-atan(Out_Dat$Ratios[i])
                b[i]<-sin(Theta[i])*R.All
                a[i]<-cos(Theta[i])*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                if(abs((Out_Dat$BetaWj[i]-Out_Dat$PredIVW[i]))>abs((Out_Dat$BetaWj[i]-Out_Dat$PredEgger[i]))){
                  
                  B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredIVW[i],linetype="solid",colour="#56B4E9")
                  B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredEgger[i],linetype="solid",colour="#D55E00")
                  
                }
                if(abs((Out_Dat$BetaWj[i]-Out_Dat$PredIVW[i]))<abs((Out_Dat$BetaWj[i]-Out_Dat$PredEgger[i]))){
                  B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredEgger[i],linetype="solid",colour="#D55E00")
                  B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredIVW[i],linetype="solid",colour="#56B4E9")
                  
                }
              }
            }
            
            
          }
          
          #If axes scales should be fixed
          if(SCALEMATCH=="YES"){
            
            #If no outliers are present
            if(sum(Out_Combine)==0){
              
              
              #Define minimum and maximum values of y axis for scaling
              Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generate a plot showing both IVW and MR-Egger estimates, using the full scale and all variants
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="Radial Estimates")+ geom_point(aes(colour="Variant"))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*5))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+coord_fixed()+
                geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
                geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
                geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
                geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+
                                   Egger.Intercept,colour="MR-Egger"))+scale_color_manual(breaks=c("Variant","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000"))+ theme(legend.title=element_blank())                   
              
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Add lines showing ratio estimates for each individual variant
              for(i in 1:length(BYG)){
                Theta<-atan(Ratios)
                b<-sin(Theta)*R.All
                a<-cos(Theta)*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
              }
              
              #Indicate that no outliers are detected
              print("No significant Outliers")
            }
            
            #If outliers are present
            if(sum(Out_Combine)>0){
              
              #Create a subset containing only outlying variants
              Out_Dat<-subset(Dat, Outliers == "Outlier")
              
              #Redefine IVW, MR-Egger and scale radii
              maximal<-atan(max(abs(Out_Dat$BetaWj))/max(abs(Out_Dat$Wj)))
              R.All<-max(abs(Out_Dat$BetaWj))/sin(maximal)
              R.IVW<-R.All+(R.All/12)
              R.Egger<-R.All+(R.All/4)
              
              #Redefine general label scaling parameter
              Label.Scaling<-R.IVW/24
              
              #Redefine MR-Egger label scaling parameter
              Label.Scaling_Egg<-R.Egger/8
              
              #Redefine values for plotting IVW, MR-Egger and scale circles
              circledat_EGGEREST <- circleFun(c(0, Egger.Intercept),R.Egger,npoints = length(BYG),atan(Eggerslope_CI[1]),atan(Eggerslope_CI[2]))
              cxEgger<-circledat_EGGEREST$x
              cyEgger<-circledat_EGGEREST$y
              
              #Plot circle for overall scale
              circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),atan(Scalemin),atan(Scalemax))
              cxAll<-circledat_ALLEST$x
              cyAll<-circledat_ALLEST$y
              
              circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
              cxIVW<-circledat_IVWEST$x
              cyIVW<-circledat_IVWEST$y
              
              #Define minimum and maximum values of y axis for scaling
              Y_MIN<-sin(atan(Scalemin))*(R.IVW+Label.Scaling)
              Y_MAX<-sin(atan(Scalemax))*(R.IVW+Label.Scaling)
              
              Y_Range<-Scalemax-Scalemin
              
              #Generates a plot showing the MR-Egger and IVW estimates, showing only outliers and using the full scale
              B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
              )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types))+geom_path(aes(x=cxAll,y=cyAll))+geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
                geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
                geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
                theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                  axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
                scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*5))),expand=c(0,0))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+coord_fixed()+
                geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
                geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
                geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
                geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
                geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+scale_color_manual(breaks=c("IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="white","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())
              
              
              if(max(abs(Ratios))>2){
                Y_Scale<-Y_Range/6
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[6]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[6]))*(R.All+Label.Scaling)), label=round(Y_Scale[6]),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[7]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[7]))*(R.All+Label.Scaling)), label=round(Y_Scale[7]),size=2.5)
                
              }
              
              if(max(abs(Ratios))<2){
                
                Y_Scale<-Y_Range/4
                Y_Scale<-seq(from=Scalemin,to=Scalemax,by=Y_Scale)
                B<-B + geom_text(x=(cos(atan(Y_Scale[1]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[1]))*(R.All+Label.Scaling)), label=round(Y_Scale[1],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[2]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[2]))*(R.All+Label.Scaling)), label=round(Y_Scale[2],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[3]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[3]))*(R.All+Label.Scaling)), label=round(Y_Scale[3],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[4]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[4]))*(R.All+Label.Scaling)), label=round(Y_Scale[4],digits=1),size=2.5)
                B<-B + geom_text(x=(cos(atan(Y_Scale[5]))*(R.All+Label.Scaling)), y=(sin(atan(Y_Scale[5]))*(R.All+Label.Scaling)), label=round(Y_Scale[5],digits=1),size=2.5)
                
              }
              
              #Draw individual ratio estimates for each outlying variant, and Q statistic indicators of respective colours      
              Theta<-rep(0,length(Out_Dat$Outliers))
              a<-rep(0,length(Out_Dat$Outliers))
              b<-rep(0,length(Out_Dat$Outliers))
              for(i in 1:length(Out_Dat$Outliers)){
                Theta[i]<-atan(Out_Dat$Ratios[i])
                b[i]<-sin(Theta[i])*R.All
                a[i]<-cos(Theta[i])*R.All
                B<- B + geom_segment(x = 0, xend = a[i], y = 0, yend =b[i],linetype="dotted",colour="grey75")
                if(abs((Out_Dat$BetaWj[i]-Out_Dat$PredIVW[i]))>abs((Out_Dat$BetaWj[i]-Out_Dat$PredEgger[i]))){
                  
                  B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredIVW[i],linetype="solid",colour="#56B4E9")
                  B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredEgger[i],linetype="solid",colour="#D55E00")
                  
                }
                if(abs((Out_Dat$BetaWj[i]-Out_Dat$PredIVW[i]))<abs((Out_Dat$BetaWj[i]-Out_Dat$PredEgger[i]))){
                  B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredEgger[i],linetype="solid",colour="#D55E00")
                  B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredIVW[i],linetype="solid",colour="#56B4E9")
                  
                }
              }
            }
          }
        }
        
        #If the full scale is not to be used
        if(FULL=="NO") {
          
          #If no outliers are present
          if(sum(Out_Combine)==0){
            
            MINCILIM<-min(sin(atan(confint(IVW.Model)[1]))*R.IVW + Label.Scaling,sin(atan(Eggerslope_CI[1]))*(R.Egger +Label.Scaling_Egg-Egger.Intercept))
            
            MAXCILIM<-max(sin(atan(confint(IVW.Model)[2]))*R.IVW + Label.Scaling,sin(atan(Eggerslope_CI[2]))*(R.Egger +Label.Scaling_Egg-Egger.Intercept))
            
            if(min(BetaWj) < MINCILIM){
              Y_MIN<-min(BetaWj)
            } 
            
            if(min(BetaWj) > MINCILIM){
              Y_MIN<-MINCILIM
            } 
            
            if(max(BetaWj) >MAXCILIM){
              Y_MAX<-max(BetaWj)
            } 
            
            if(max(BetaWj) <MAXCILIM){
              Y_MAX<-MAXCILIM
            } 
            
            Y_Range<-Y_MAX-Y_MIN
            
            #Generate a plot showing the MR-Egger and IVW estimates, showing all variants and not using a full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour=Outliers)) +geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
              geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+
              theme_bw() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+
              geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*3))),expand=c(0,0))+geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+
                                                                                                                                    Egger.Intercept,colour="MR-Egger"))+scale_y_continuous(limits = c(Y_MIN,Y_MAX))+scale_color_manual(breaks=c("Variant","Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="#000000","Outlier"="#E69F00"))+
              theme(legend.title=element_blank())
            
            #Indicate that no outliers are present
            print("No significant Outliers")
          }
          
          #If outliers are present
          if(sum(Out_Combine)>0){
            
            #Construct a subset of data containing only outliers
            Out_Dat<-subset(Dat, Outliers == "Outlier")
            
            #Redefine IVW, MR-Egger and scale radii
            maximal<-atan(max(abs(Out_Dat$BetaWj))/max(abs(Out_Dat$Wj)))
            R.All<-max(abs(Out_Dat$BetaWj))/sin(maximal)
            R.IVW<-R.All+(R.All/12)
            R.Egger<-R.All+(R.All/4)
            
            #Redefine general scaling parameter
            Label.Scaling<-R.IVW/24
            
            #Redefine MR-Egger scaling parameter
            Label.Scaling_Egg<-R.Egger/8
            
            #Redefine plotting values for scale, MR-Egger and IVW circles
            
            circledat_EGGEREST <- circleFun(c(0, Egger.Intercept),R.Egger,npoints = length(BYG),atan(Eggerslope_CI[1]),atan(Eggerslope_CI[2]))
            cxEgger<-circledat_EGGEREST$x
            cyEgger<-circledat_EGGEREST$y
            
            circledat_ALLEST <- circleFun(c(0,0),R.All,npoints = length(BYG),-pi/2,pi/2)
            cxAll<-circledat_ALLEST$x
            cyAll<-circledat_ALLEST$y
            
            circledat_IVWEST <- circleFun(c(0,0),R.IVW,npoints = length(BYG),atan(IVW_CI[1]),atan(IVW_CI[2]))
            cxIVW<-circledat_IVWEST$x
            cyIVW<-circledat_IVWEST$y
            
            
            
            
            MINCILIM<-min(sin(atan(confint(IVW.Model)[1]))*R.IVW + Label.Scaling,sin(atan(Eggerslope_CI[1]))*(R.Egger +Label.Scaling_Egg-Egger.Intercept))
            
            MAXCILIM<-max(sin(atan(confint(IVW.Model)[2]))*R.IVW + Label.Scaling,sin(atan(Eggerslope_CI[2]))*(R.Egger +Label.Scaling_Egg-Egger.Intercept))
            
            if(min(Out_Dat$BetaWj) < MINCILIM){
              Y_MIN<-min(Out_Dat$BetaWj)
            } 
            
            if(min(Out_Dat$BetaWj) > MINCILIM){
              Y_MIN<-MINCILIM
            } 
            
            if(max(Out_Dat$BetaWj) >MAXCILIM){
              Y_MAX<-max(Out_Dat$BetaWj)
            } 
            
            if(max(Out_Dat$BetaWj) <MAXCILIM){
              Y_MAX<-MAXCILIM
            } 
            
            Y_Range<-Y_MAX-Y_MIN
            
            #Generate a plot showing both IVW and MR-Egger estimates, only showing outlying variants and not using the full scale
            B<-ggplot(Dat,aes(x=Wj,y=BetaWj)
            )+labs(title="Radial Estimates")+ geom_point(aes(colour=out_types)) +geom_path(aes(x=cxIVW,y=cyIVW,colour="IVW"))+
              geom_text(x=cos(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), y=sin(atan(IVW.Slope))*(R.IVW+Label.Scaling*1.5), label=round(IVW.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[1])-0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[1],digits=3),size=3)+
              geom_text(x=cos(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),y=sin(atan(confint(IVW.Model)[2])+0.02)*(R.IVW+Label.Scaling),label=round(confint(IVW.Model)[2],digits=3),size=3)+
              geom_path(aes(x=cxEgger,y=cyEgger,colour="MR-Egger"))+
              geom_text(x=cos(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5), y=sin(atan(Egger.Slope))*(R.Egger+Label.Scaling*1.5) +Egger.Intercept,label=round(Egger.Slope,digits=3),size=4)+
              geom_text(x=cos(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[1]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[1],digits=3),size=3)+
              geom_text(x=cos(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5),y=(sin(atan(Eggerslope_CI[2]))*(R.Egger+Label.Scaling*1.5))+Egger.Intercept,label=round(Eggerslope_CI[2],digits=3),size=3)+theme_bw() +
              theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"))+ylab(expression(hat(beta)[j]~sqrt(W[j])))+xlab(expression(sqrt(W[j])))+geom_segment(aes(x = 0, xend = cos(atan(IVW.Slope))*R.IVW, y = 0, yend = sin(atan(IVW.Slope))*R.IVW,colour="IVW"))+
              scale_x_continuous(limits = c(0,cos(atan(Egger.Slope))*(R.Egger+(Label.Scaling*3))),expand=c(0,0))+geom_segment(aes(x = 0, xend = cos(atan(Egger.Slope))*R.Egger, y = Egger.Intercept, yend = (sin(atan(Egger.Slope))*R.Egger)+Egger.Intercept,colour="MR-Egger"))+
              scale_y_continuous(limits = c(Y_MIN,Y_MAX))+scale_color_manual(breaks=c("IVW and MR-Egger Outlier","MR-Egger Outlier","IVW Outlier","IVW","MR-Egger"),values=c("IVW"="#56B4E9","MR-Egger"="#D55E00","Variant"="white","IVW and MR-Egger Outlier"="#CC79A7","MR-Egger Outlier"="#E69F00","IVW Outlier"="#009E73"))+ theme(legend.title=element_blank())                   
            
            #Draw lines showing Q statistic contributions for each variant using respectvie colours
            for(i in 1:length(Out_Dat$Outliers)){
              if(abs((Out_Dat$BetaWj[i]-Out_Dat$PredIVW[i]))>abs((Out_Dat$BetaWj[i]-Out_Dat$PredEgger[i]))){
                B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredIVW[i],linetype="solid",colour="#56B4E9")
                B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredEgger[i],linetype="solid",colour="#D55E00")
              }
              if(abs((Out_Dat$BetaWj[i]-Out_Dat$PredIVW[i]))<abs((Out_Dat$BetaWj[i]-Out_Dat$PredEgger[i]))){
                B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredEgger[i],linetype="solid",colour="#D55E00")
                B<- B + geom_segment(x = Out_Dat$Wj[i], xend = Out_Dat$Wj[i], y = Out_Dat$BetaWj[i], yend =Out_Dat$PredIVW[i],linetype="solid",colour="#56B4E9")
              }
            }
          }
        }
      }
      
      #Return defined plot "B"
      return(B)
    }
    
    #Return specified plot as "A" using BOTH.Function
    A<-BOTH.Function(BXG,BYG,seBXG,seBYG,SCALE,OUTLIERS)
  }
  
  #Specify multiple return options for RadialMR function: Plot, Overall Q statistics, Outlier information andIVW/MR-Egger model estimates
  multi_return <- function() {
    Out_list <- list("plot" = A,"Q"= QDisplay, "outliers" = Q_DATA,"egger"= EstimatesEGGER,"IVW"=EstimatesIVW)
    return(Out_list) 
  }
  OUT<-multi_return()
}

