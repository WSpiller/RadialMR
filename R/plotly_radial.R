#' plotly_radial
#'
#' A function for producing interactive radial IVW and MR-Egger plots individually.The function utilises the output from the \code{IVW_radial} and \code{egger_radial} functions.
#'
#' @param r_object An object of class \code{"IVW"} or \code{"egger"}.
#' @return A plotly object containing a radial plot of either the IVW or MR-Egger estimates. Hovering the mouse over individual datapoints will highlight the corresponding SNP identification number for that observation.
#'@author Wes Spiller; Jack Bowden.
#'@references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.
#'@export
#'@examples
#'
#'plotly_radial(r_object) 

plotly_radial<-function(r_object,TEST){

  if(missing(TEST)) {
    TEST<-FALSE
  }
  
library(plotly)
  
  if(length(r_object)==13){
    
    r_object$coef<-c(r_object$coef[2,])
    r_object$coef<-as.numeric(r_object$coef)
    
    circle.radial <- function(center = c(0,0), radius, num.points, START,END){
      R = radius
      tt <- seq(START,END,length.out = num.points)
      #Generates x-axis values for circle
      xx <- center[1] + R * cos(tt)
      #Generates y-axis values for circle
      yy <- center[2] + R * sin(tt)
      return(data.frame(x = xx, y = yy))
    }
    
    maxWj<-max(r_object$data$Wj)
    Wjcor<-r_object$data[r_object$data$Wj==maxWj,]
    R.IVW<-Wjcor$Wj^2+Wjcor$BetaWj^2
    R.IVW<-sqrt(R.IVW)*1.05
    IVW_circlemin<-r_object$confint[1]
    IVW_circlemax<-r_object$confint[2]
    
    circledat_IVWEST <- circle.radial(c(0,0),R.IVW,100,atan(IVW_circlemin),atan(IVW_circlemax))
    cxIVW<-circledat_IVWEST$x
    cyIVW<-circledat_IVWEST$y
    
    ay <- list(
      zeroline = FALSE,
      showline = FALSE,
      mirror = "ticks",
      gridcolor = toRGB("white"),
      gridwidth = 0,
      zerolinecolor = toRGB("black"),
      zerolinewidth = 1,
      linecolor = toRGB("white"),
      linewidth = 0,
      title = "Beta*sqrt(Wj)"
    )
    
    ax <- list(
      zeroline = TRUE,
      showline = FALSE,
      mirror = "ticks",
      gridcolor = toRGB("white"),
      gridwidth = 0,
      zerolinecolor = toRGB("black"),
      zerolinewidth = 1,
      linecolor = toRGB("white"),
      linewidth = 0,
      title= "sqrt(Wj)"
    )
    
    if(nrow(r_object$data[r_object$data$Outliers=="Outlier",]) == 0){
      
      T_PLOT<- plot_ly(r_object$data, mode="marker", type = 'scatter') %>%
        add_trace(x = r_object$data[r_object$data$Outliers=="Variant",]$Wj , y = r_object$data[r_object$data$Outliers=="Variant",]$BetaWj, name = 'Variant', mode = 'markers', marker= list(color="black"),text= r_object$data[r_object$data$Outliers=="Variant",]$SNP, hoverinfo = 'text') %>%
        add_trace(x = cxIVW , y = cyIVW, mode = 'line',showlegend=FALSE,name="IVW_CI", line= list(color = "#56B4E9",size=1)) %>%
        layout(xaxis = ax,yaxis = ay) %>%
        add_segments(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,name = "IVW", line= list(color = "#56B4E9",size=1))%>%
        add_annotations(x = c(cxIVW[1],cos(atan(r_object$coef[1]))*R.IVW,cxIVW[100]),
                        y = c(cyIVW[1],sin(atan(r_object$coef[1]))*R.IVW,cyIVW[100]),
                        text = c(round(r_object$confint[1],digits=3),round(r_object$coef[1],digits=3),round(r_object$confint[2],digits=3)),
                        xref = "x",
                        yref = "y",
                        showarrow = TRUE,
                        arrowhead = 0,
                        arrowsize = 0,
                        ax = 25,
                        ay = 1)
      
    }
    
    if(nrow(r_object$data[r_object$data$Outliers=="Outlier",]) > 0){
      
      T_PLOT<- plot_ly(r_object$data, mode="marker", type = 'scatter') %>%
        add_trace(x = r_object$data[r_object$data$Outliers=="Variant",]$Wj , y = r_object$data[r_object$data$Outliers=="Variant",]$BetaWj, name = 'Variant', mode = 'markers', marker= list(color="black"),text= r_object$data[r_object$data$Outliers=="Variant",]$SNP, hoverinfo = 'text') %>%
        add_trace(x = r_object$data[r_object$data$Outliers=="Outlier",]$Wj , y = r_object$data[r_object$data$Outliers=="Outlier",]$BetaWj, name = 'Outlier', mode = 'markers', marker= list(color="#E69F00"),text= r_object$data[r_object$data$Outliers=="Outlier",]$SNP, hoverinfo = 'text') %>%
        add_trace(x = cxIVW , y = cyIVW, mode = 'line',showlegend=FALSE,name="IVW_CI", line= list(color = "#56B4E9",size=1)) %>%
        layout(xaxis = ax,yaxis = ay) %>%
        add_segments(x = 0, xend = cos(atan(r_object$coef[1]))*R.IVW, y = 0, yend = sin(atan(r_object$coef[1]))*R.IVW,name = "IVW", line= list(color = "#56B4E9",size=1))%>%
        add_annotations(x = c(cxIVW[1],cos(atan(r_object$coef[1]))*R.IVW,cxIVW[100]),
                        y = c(cyIVW[1],sin(atan(r_object$coef[1]))*R.IVW,cyIVW[100]),
                        text = c(round(r_object$confint[1],digits=3),round(r_object$coef[1],digits=3),round(r_object$confint[2],digits=3)),
                        xref = "x",
                        yref = "y",
                        showarrow = TRUE,
                        arrowhead = 0,
                        arrowsize = 0,
                        ax = 25,
                        ay = 1)
      
    }
    
  }
  
  if(length(r_object)==6){
    
    if(TEST==T){
    
    circle.radial <- function(center = c(0,r_object$coef[1,1]), radius, num.points, START,END){
      R = radius
      tt <- seq(START,END,length.out = num.points)
      #Generates x-axis values for circle
      xx <- center[1] + R * cos(tt)
      #Generates y-axis values for circle
      yy <- center[2] + R * sin(tt)
      return(data.frame(x = xx, y = yy))
    }
    
    maxWj<-max(r_object$data$Wj)
    Wjcor<-r_object$data[r_object$data$Wj==maxWj,]
    R.Egger<-Wjcor$Wj^2+Wjcor$BetaWj^2
    R.Egger<-sqrt(R.Egger)*1.05
    Egger_circlemin<-r_object$confint[1]
    Egger_circlemax<-r_object$confint[2]
    
    circledat_EggerEST <- circle.radial(c(0,r_object$coef[1,1]),R.Egger,100,atan(Egger_circlemin),atan(Egger_circlemax))
    cxEgger<-circledat_EggerEST$x
    cyEgger<-circledat_EggerEST$y
    
    ay <- list(
      zeroline = FALSE,
      showline = FALSE,
      mirror = "ticks",
      gridcolor = toRGB("white"),
      gridwidth = 0,
      zerolinecolor = toRGB("black"),
      zerolinewidth = 1,
      linecolor = toRGB("white"),
      linewidth = 0,
      title = "Beta*sqrt(Wj)"
    )
    
    ax <- list(
      zeroline = TRUE,
      showline = FALSE,
      mirror = "ticks",
      gridcolor = toRGB("white"),
      gridwidth = 0,
      zerolinecolor = toRGB("black"),
      zerolinewidth = 1,
      linecolor = toRGB("white"),
      linewidth = 0,
      title= "sqrt(Wj)"
    )
    
    if(nrow(r_object$data[r_object$data$Outliers=="Outlier",]) == 0){
      
      T_PLOT<- plot_ly(r_object$data, mode="marker", type = 'scatter') %>%
        add_trace(x = r_object$data[r_object$data$Outliers=="Variant",]$Wj , y = r_object$data[r_object$data$Outliers=="Variant",]$BetaWj, name = 'Variant', mode = 'markers', marker= list(color="black"),text= r_object$data[r_object$data$Outliers=="Variant",]$SNP, hoverinfo = 'text') %>%
        add_trace(x = cxEgger , y = cyEgger, mode = 'line',showlegend=FALSE,name="Egger_CI", line= list(color = "#D55E00",size=1)) %>%
        layout(xaxis = ax,yaxis = ay) %>%
        add_segments(x = 0, xend = cos(atan(r_object$coef[2,1]))*R.Egger, y = r_object$coef[1,1], yend = (sin(atan(r_object$coef[2,1]))*R.Egger)+r_object$coef[1,1],name = "Egger", line= list(color = "#D55E00",size=1))%>%
        add_annotations(x = c(cxEgger[1],cos(atan(r_object$coef[2,1]))*R.Egger,cxEgger[100]),
                        y = c(cyEgger[1],(sin(atan(r_object$coef[2,1]))*R.Egger)+r_object$coef[1,1],cyEgger[100]),
                        text = c(round(r_object$confint[1],digits=3),round(r_object$coef[2,1],digits=3),round(r_object$confint[2],digits=3)),
                        xref = "x",
                        yref = "y",
                        showarrow = TRUE,
                        arrowhead = 0,
                        arrowsize = 0,
                        ax = 25,
                        ay = 1)
      
    }
    
    if(nrow(r_object$data[r_object$data$Outliers=="Outlier",]) > 0){
      
      T_PLOT<- plot_ly(r_object$data, mode="marker", type = 'scatter') %>%
        add_trace(x = r_object$data[r_object$data$Outliers=="Variant",]$Wj , y = r_object$data[r_object$data$Outliers=="Variant",]$BetaWj, name = 'Variant', mode = 'markers', marker= list(color="black"),text= r_object$data[r_object$data$Outliers=="Variant",]$SNP, hoverinfo = 'text') %>%
        add_trace(x = r_object$data[r_object$data$Outliers=="Outlier",]$Wj , y = r_object$data[r_object$data$Outliers=="Outlier",]$BetaWj, name = 'Outlier', mode = 'markers', marker= list(color="#E69F00"),text= r_object$data[r_object$data$Outliers=="Outlier",]$SNP, hoverinfo = 'text') %>%
        add_trace(x = cxEgger , y = cyEgger, mode = 'line',showlegend=FALSE,name="Egger_CI", line= list(color = "#D55E00",size=1)) %>%
        layout(xaxis = ax,yaxis = ay) %>%
        add_segments(x = 0, xend = cos(atan(r_object$coef[2,1]))*R.Egger, y = r_object$coef[1,1], yend = (sin(atan(r_object$coef[2,1]))*R.Egger)+r_object$coef[1,1],name = "Egger", line= list(color = "#D55E00",size=1))%>%
        add_annotations(x = c(cxEgger[1],cos(atan(r_object$coef[2,1]))*R.Egger,cxEgger[100]),
                        y = c(cyEgger[1],(sin(atan(r_object$coef[2,1]))*R.Egger)+r_object$coef[1,1],cyEgger[100]),
                        text = c(round(r_object$confint[1],digits=3),round(r_object$coef[2,1],digits=3),round(r_object$confint[2],digits=3)),
                        xref = "x",
                        yref = "y",
                        showarrow = TRUE,
                        arrowhead = 0,
                        arrowsize = 0,
                        ax = 25,
                        ay = 1)
      
    }
    
    }
    
    if(TEST==FALSE){
      T_PLOT<- as.character("Coming soon")
      
    }
    
  }
  
  return(T_PLOT)
  
}
