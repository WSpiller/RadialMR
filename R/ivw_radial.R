#' ivw_radial
#'
#' Fits a radial inverse variance weighted (IVW) model using first order, second order, or modified second order weights. Outliers are identified using a significance threshold specified by the user. The function returns an object of class \code{"IVW"}, containing regression estimates, a measure of total heterogeneity using Cochran's Q statistic, the individual contribution to overall heterogeneity of each variant, and a data frame for use in constructing the radial plot. 
#' 
#' @param r_input A formatted data frame using the \code{format_radial} function.
#' @param alpha A value specifying the statistical significance threshold for identifying outliers (\code{0.05} specifies a p-value threshold of 0.05).
#' @param weights A value specifying the inverse variance weights used to calculate IVW estimate and Cochran's Q statistic. By default modified second order weights are used, but one can choose to select first order (\code{1}), second order (\code{2}) or modified second order weights (\code{3}).
#' @param summary Indicates whether a summary of the radial IVW model should be printed (default= \code{TRUE}) or withheld (\code{FALSE}).
#' @return An object of class \code{"IVW"} containing the following components:\describe{
#' \item{\code{coef}}{The estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{qstatistic}}{Cochran's Q statistic for overall heterogeneity.}
#' \item{\code{df}}{Degrees of freedom. This is equal to the number of variants -1 when fitting the radial IVW model.}
#' \item{\code{outliers}}{A data frame containing variants identified as outliers, with respective Q statistics, chi-squared tests and SNP identification.}
#' \item{\code{data}}{A data frame containing containing SNP IDs, inverse variance weights, the product of the inverse variance weight and ratio estimate for each variant, contribution to overall heterogeneity with corresponding p-value, and a factor indicator showing outlier status.}
#' \item{\code{confint}}{A vector giving lower and upper confidence limits for the radial IVW effect estimate.}
#'}
#'@author Wes Spiller; Jack Bowden.
#'@references Bowden, J. et al., 2017. "Improving the visualisation, interpretation and analysis of two-sample summary data Mendelian randomization via the radial plot and radial regression." BioRxiv. doi: https://doi.org/10.1101/200378. Available at: \url{https://www.biorxiv.org/content/early/2017/10/11/200378}
#'@export
#'@examples
#'
#' ivw_radial(r_input,0.05,1)
#' 


ivw_radial<-function(r_input,alpha,weights,summary){
  
  #Define ratio estimates
  Ratios<-r_input[,3]/r_input[,2]
  
  cat()
  
  if(missing(alpha)) {
    alpha<-0.05
    warning("Significance threshold for outlier detection not specified: Adopting a 95% threshold")
  }
  
  if(missing(weights)) {
    weights<-3
    warning("Weights not specified: Adopting modified second-order weights")
  }
  
  if(missing(summary)) {
    summary<-TRUE
  }
  
  if(weights==1){
    
    #Define inverse variance weights
    W<-((r_input[,2]^2)/(r_input[,5]^2))
    
  }
  
  if(weights==2){
    
    #Define inverse variance weights
    W<-((r_input[,5]^2/r_input[,2]^2)+((r_input[,3]^2*r_input[,4]^2)/r_input[,2]^4))^-1
    
  }
  
  if(weights==3){
    
    #Define inverse variance weights
    W<-((r_input[,2]^2)/(r_input[,5]^2))
    
  }
  
  #Define vector of squared weights
  Wj<-sqrt(W)
  
  #Define vector of weights * ratio estimates
  BetaWj<-Ratios*Wj
  
  #Define IVW Model
  IVW.Model<-lm(BetaWj~-1+Wj)
  
  #Define output of IVW.Model fit
  EstimatesIVW<-summary(lm(IVW.Model))
  
  #Define slope parameter for IVW.Model
  IVW.Slope<-EstimatesIVW$coefficients[1]
  
  #Define standard error for slope parameter of IVW.Model
  IVW.SE<-EstimatesIVW$coefficients[2]
  
  #Define confidence interval for IVW.Model slope parameter
  IVW_CI<-confint(IVW.Model)
  
  #Define Q statistic for each individual variant
  Qj<-W*(Ratios-IVW.Slope)^2
  
  #Define total Q statistic
  Total_Q<-sum(Qj)
  
  #Perform chi square test for overall Q statistic
  Total_Q_chi<-pchisq(Total_Q,length(r_input[,2])-1,lower.tail = FALSE)
  
  if(weights==3){
    #Define inverse variance weights
    W<- ((r_input[,5]^2+(IVW.Slope^2*r_input[,4]^2))/r_input[,2]^2)^-1
    
    #Define vector of squared weights
    Wj<-sqrt(W)
    
    #Define vector of weights * ratio estimates
    BetaWj<-Ratios*Wj
    
    #Define IVW Model
    IVW.Model<-lm(BetaWj~-1+Wj)
    
    #Define output of IVW.Model fit
    EstimatesIVW<-summary(lm(BetaWj~-1+Wj))
    
    #Define slope parameter for IVW.Model
    IVW.Slope<-EstimatesIVW$coefficients[1]
    
    #Define standard error for slope parameter of IVW.Model
    IVW.SE<-EstimatesIVW$coefficients[2]
    
    #Define confidence interval for IVW.Model slope parameter
    IVW_CI<-confint(IVW.Model)
    
    #Define Q statistic for each individual variant
    Qj<-W*(Ratios-IVW.Slope)^2
    
    #Define total Q statistic
    Total_Q<-sum(Qj)
    
    #Perform chi square test for overall Q statistic
    Total_Q_chi<-pchisq(Total_Q,length(r_input[,2])-1,lower.tail = FALSE)
    
  }
  
  #Define a placeholder vector of 0 values for chi square tests
  Qj_Chi<-0
  
  #Perform chi square tests for each Qj value
  for(i in 1:length(Qj)){
    Qj_Chi[i]<-pchisq(Qj[i],1,lower.tail = FALSE)
  }
  
  #Define dataframe with SNP IDs and outlier statistics
  r_input$Qj<-Qj
  r_input$Qj_Chi<-Qj_Chi
  
  #Define a placeholder vector of 0 values for identifying outliers
  Out_Indicator<-rep(0,length(r_input[,2]))
  
  #For loop defining outlier indicator variable.
  for( i in 1:length(r_input[,2])){
    if(Qj_Chi[i]<alpha){
      Out_Indicator[i]<-1
    }
  }
  
  #Include the outlier identifier cariable in the dataframe and define as a factor
  r_input$Outliers<-factor(Out_Indicator)
  levels(r_input$Outliers)[levels(r_input$Outliers)=="0"] <- "Variant"
  levels(r_input$Outliers)[levels(r_input$Outliers)=="1"] <- "Outlier"
  
  #If no outliers are present indicate this is the case
  if(sum(Out_Indicator==0)){
    outlier_status<-"No significant outliers"
    outtab<-"No significant outliers"
    
  }
  
  #If outliers are present produce dataframe containing individual Q statistics and chi square values for each outlier
  if(sum(Out_Indicator>0)){
    outlier_status<-"Outliers detected"
    
    #Generate a subset of data containing only outliers
    Out_Dat<-subset(r_input, Outliers == "Outlier")
    
    #Construct a datafrae containing SNP IDs, Q statistics and chi-square values for each outlying variant
    outtab<-data.frame(Out_Dat[,1],Out_Dat$Qj,Out_Dat$Qj_Chi)

    #Redefine column names
    colnames(outtab) = c("SNP","Q_statistic","p.value")
    
  }
  
  if(summary==TRUE){
    
  # Print a few summary elements that are common to both lm and plm model summary objects
  cat("\n")
  
  cat("Radial IVW\n")
  
  cat("\n")
  
  print(coef(EstimatesIVW))
  
  cat("\nResidual standard error:", round(EstimatesIVW$sigma,3), "on", EstimatesIVW$df[2], "degrees of freedom")
  
  cat("\n")
  
  cat(paste(c("\nF-statistic:", " on"," and"), round(EstimatesIVW$fstatistic,2), collapse=""),
      "DF, p-value:",
      format.pval(pf(EstimatesIVW$fstatistic[1L], EstimatesIVW$fstatistic[2L], EstimatesIVW$fstatistic[3L], 
                     lower.tail = FALSE), digits=3))
  
  cat("\n")
  
  cat("Q-Statistic for heterogeneity:",Total_Q, "on", length(r_input[,2])-1, "DF",",", "p-value:" , Total_Q_chi)
  
  cat("\n")
  
  cat("\n",outlier_status,"\n")
  
  cat("\n")
  
  }
  
  out_data<-data.frame(r_input[,1],r_input[,6],r_input[,7],r_input[,8])
  out_data$Wj<-Wj
  out_data$BetaWj<-BetaWj
  out_data<-out_data[c(1,5,6,2,3,4)]
  names(out_data)<-c("SNP","Wj","BetaWj","Qj","Qj_Chi","Outliers")
  
  multi_return <- function() {
    Out_list <- list("coef" = EstimatesIVW$coef,"qstatistic"= Total_Q, "df" = length(r_input[,2])-1, "outliers" = outtab, "data" = out_data, "confint" = confint(IVW.Model))
    class(Out_list)<-"IVW"
    
    return(Out_list) 
  }
  OUT<-multi_return()
 
  
}