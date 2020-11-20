#' egger_radial
#'
#' Fits a radial MR-Egger model using first order, second order, or modified second order weights. Outliers are identified using a significance threshold specified by the user. The function returns an object of class \code{"egger"}, containing regression estimates, a measure of total heterogeneity using Rucker's Q statistic, the individual contribution to overall heterogeneity of each variant, and a data frame for use in constructing the radial plot.
#'
#' @param r_input A formatted data frame using the \code{format_radial} function.
#' @param alpha A value specifying the statistical significance threshold for identifying outliers (\code{0.05} specifies a p-value threshold of 0.05).
#' @param weights A value specifying the inverse variance weights used to calculate the MR-Egger estimate and Rucker's Q statistic. By default modified second order weights are used, but one can choose to select first order (\code{1}), second order (\code{2}) or modified second order weights (\code{3}).
#' @param summary A logical argument (\code{T} or \code{F}) indicating whether a summary of results and heterogeneity should be presented (default= \code{TRUE}).
#' @return An object of class \code{"egger"} containing the following components:\describe{
#' \item{\code{coef}}{A matrix giving the intercept and slope coefficient, corresponding standard errors, t-statistics, and (two-sided) p-values.}
#' \item{\code{qstatistic}}{Rucker's Q statistic for overall heterogeneity.}
#' \item{\code{df}}{Degrees of freedom. This is equal to the number of variants -2 when fitting the radial MR-Egger model.}
#' \item{\code{outliers}}{A data frame containing variants identified as outliers, with respective Q statistics, chi-squared tests and SNP identification.}
#' \item{\code{data}}{A data frame containing containing SNP IDs, inverse variance weights, the product of the inverse variance weight and ratio estimate for each variant, contribution to overall heterogeneity with corresponding p-value, and a factor indicator showing outlier status.}
#' \item{\code{confint}}{A vector giving lower and upper confidence limits for the radial MR-Egger effect estimate.}
#'
#'}
#' @author Wes Spiller; Jack Bowden; Tom Palmer.
#' @references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.
#' @importFrom stats lm coef confint optimize pchisq pf pnorm pt qchisq qt sd
#' @export
#' @examples
#'
#'ldl.dat <- data_radial[data_radial[,10]<5*10^-8,]
#'
#'ldl.fdat<-format_radial(ldl.dat[,6], ldl.dat[,9],
#'               ldl.dat[,15], ldl.dat[,21],
#'               ldl.dat[,1])
#'              
#' egger_radial(ldl.fdat, 0.05, 1,T)


egger_radial<-function(r_input,alpha,weights,summary){
  
  # Perform check that r_input has been formatted using format_radial function
  
  if(!("rmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "rmr_format", please resave the object with the output of format_radial().')
    }
  
  # Provide default value for outlier significance threshold if not provided
  if(missing(alpha)) {
    alpha<-0.05
    warning("Significance threshold for outlier detection not specified: Adopting a 95% threshold")
  }
  
  # Set default weighting to modified second order weights if not provided
  if(missing(weights)) {
    weights<-3
    warning("Weights not specified: Adopting modified second-order weights")
  }
  
  if(missing(summary)) {
    summary<-T
  }
  
  # Calculate ratio estimates
  Ratios<-r_input[,3]/r_input[,2]
  
  # Calculate approximate F-statistic for each variant
  F<- r_input[,2]^2/r_input[,4]^2
  
  # Define mean F-statistic across all variants
  mf<- mean(F)
  
  # Calculate first order weights
  if(weights==1){
    W<-((r_input[,2]^2)/(r_input[,5]^2))
  }
  
  # Calculate second order weights
  if(weights==2){
    W<-((r_input[,5]^2/r_input[,2]^2)+((r_input[,3]^2*r_input[,4]^2)/r_input[,2]^4))^-1
  }
  
  # Initially calculate first order weights for downstream use of modified second order weights
  if(weights==3){
    W<-((r_input[,2]^2)/(r_input[,5]^2))
  }
  
  # Define vector of square root weights
  Wj<-sqrt(W)
  
  # Create vector of ratio estimates multiplied by given square root weightings
  BetaWj<-Ratios*Wj
  
  # Define Egger Model
  Egger.Model<-lm(BetaWj~Wj)
  
  # Save summary of Radial Egger model
  EstimatesEGGER<-summary(Egger.Model)
  
  # Define intercept of Radial Egger model
  Egger.Intercept<-EstimatesEGGER$coefficients[1]
  
  # Define causal effect estimate from Radial Egger model
  Egger.Slope<-EstimatesEGGER$coefficients[2]
  
  # Define standard error for intercept of Radial Egger model
  Eggerintercept.SE<-EstimatesEGGER$coefficients[3]
  
  # Define standard error for causal effect from Radial Egger model
  Eggerslope.SE<-EstimatesEGGER$coefficients[4]
  
  # Define 95% confidence interval for intercept from Radial Egger model
  Eggerint_CI<-as.numeric(confint(Egger.Model)[1,])
  
  # Define 95% confidence interval causal effect estimate from Radial Egger model
  Eggerslope_CI<-as.numeric(confint(Egger.Model)[2,])
  
  # Calculate Q statistic for each individual variant
  Qj<-W*(Ratios-(Egger.Intercept/Wj)-Egger.Slope)^2
  
  # Calculate total Q statistic as the sum of individual Q contributions Qj
  Total_Q<-sum(Qj)
  
  # Perform chi square test with respect to global Q statistic Total_Q
  Total_Q_chi<-pchisq(Total_Q,length(r_input[,2])-2,lower.tail = FALSE)
  
  # Perform additional analyses to calculate modified second order weights
  if(weights==3){
    W<- ((r_input[,5]^2+(Egger.Slope^2*r_input[,4]^2))/r_input[,2]^2)^-1
    
    # Create vector of squareroot weights
    Wj<-sqrt(W)
    
    # Create vector of ratio estimates multiplied by given squareroot weightings
    BetaWj<-Ratios*Wj
    
    # Define Radial Egger model
    Egger.Model<-lm(BetaWj~Wj)
    
    # Save summary of Radial Egger model
    EstimatesEGGER<-summary(Egger.Model)
    
    # Define intercept of Radial Egger model
    Egger.Intercept<-EstimatesEGGER$coefficients[1]
    
    # Define causal effect estimate from Radial Egger model
    Egger.Slope<-EstimatesEGGER$coefficients[2]
    
    #Define standard error for intercept of from Radial Egger model
    Eggerintercept.SE<-EstimatesEGGER$coefficients[3]
    
    #Define standard error for causal effect estimate from Radial Egger model
    Eggerslope.SE<-EstimatesEGGER$coefficients[4]
    
    # Define 95% confidence interval for intercept from Radial Egger model
    Eggerint_CI<-as.numeric(confint(Egger.Model)[1,])
    
    # Define 95% confidence interval forcausal effect estimate from Radial Egger model
    Eggerslope_CI<-as.numeric(confint(Egger.Model)[2,])
    
    # Calculate Q statistic for each individual variant
    Qj<-W*(Ratios-(Egger.Intercept/Wj)-Egger.Slope)^2
    
    # Calculate total Q statistic as the sum of individual Q contributions Qj
    Total_Q<-sum(Qj)
    
    # Perform chi square test with respect to global Q statistic Total_Q
    Total_Q_chi<-pchisq(Total_Q,length(r_input[,2])-1,lower.tail = FALSE)
    
  }
  
  # Define a placeholder vector of 0 values for chi square tests
  Qj_Chi<-0
  
  # Perform chi square tests for each Q contribution Qj
  for(i in 1:length(Qj)){
    Qj_Chi[i]<-pchisq(Qj[i],1,lower.tail = FALSE)
  }
  
  # Create data frame with SNP IDs and outlier information
  r_input$Qj<-Qj
  r_input$Qj_Chi<-Qj_Chi
  
  # Define a placeholder vector of 0 values for outlier status variable
  Out_Indicator<-rep(0,length(r_input[,2]))
  
  # Include value of 1 indicating positive outlier status for given sig.threshold
  for( i in 1:length(r_input[,2])){
    if(Qj_Chi[i]<alpha){
      Out_Indicator[i]<-1
    }
  }
  
  # Include the outlier status variable in the data frame as a factor
  r_input$Outliers<-factor(Out_Indicator)
  levels(r_input$Outliers)[levels(r_input$Outliers)=="0"] <- "Variant"
  levels(r_input$Outliers)[levels(r_input$Outliers)=="1"] <- "Outlier"
  
  # Provide indication if no outliers are present
  if(sum(Out_Indicator==0)){
    outlier_status<-"No significant outliers"
    outtab<-"No significant outliers"
  }
  
  # If outliers are present produce data frame with information on outliers
  if(sum(Out_Indicator>0)){
    outlier_status<-"Outliers detected"
    
    # Generate a subset of data containing only outliers
    Out_Dat<-subset(r_input, Outliers == "Outlier")
    
    # Construct a data frame containing SNP IDs, Q statistics and chi-square values for each outlying variant
    outtab<-data.frame(Out_Dat[,1],Out_Dat$Qj,Out_Dat$Qj_Chi)
    
    # Redefine column names
    colnames(outtab) = c("SNP","Q_statistic","p.value")
    
  }
  
  if(summary==TRUE){
    
    # Print a few summary elements that are common to both lm and plm model summary objects
    cat("\n")
    
    cat("Radial MR-Egger\n")
    
    cat("\n")
    
    print(coef(EstimatesEGGER))
    
    cat("\nResidual standard error:", round(EstimatesEGGER$sigma,3), "on", EstimatesEGGER$df[2], "degrees of freedom")
    
    cat("\n")
    
    cat(paste(c("\nF-statistic:", " on"," and"), round(EstimatesEGGER$fstatistic,2), collapse=""),
        "DF, p-value:",
        format.pval(pf(EstimatesEGGER$fstatistic[1L], EstimatesEGGER$fstatistic[2L], EstimatesEGGER$fstatistic[3L],
                       lower.tail = FALSE), digits=3))
    
    cat("\n")
    
    cat("Q-Statistic for heterogeneity:",Total_Q, "on", length(r_input[,2])-2, "DF",",", "p-value:" , Total_Q_chi)
    
    cat("\n")
    
    cat("\n",outlier_status,"\n")
    
    cat("\n")
  }
  
  # Create data frame containing information used to calculate radial estimates and determine outlier status
  out_data<-data.frame(r_input[,1],r_input[,6],r_input[,7],r_input[,8])
  out_data$Wj<-Wj
  out_data$BetaWj<-BetaWj
  out_data<-out_data[c(1,5,6,2,3,4)]
  names(out_data)<-c("SNP","Wj","BetaWj","Qj","Qj_Chi","Outliers")
  
  multi_return <- function() {
    Out_list <- list("coef" = EstimatesEGGER$coef,"qstatistic"= Total_Q, "df" = length(r_input[,2])-1, "outliers" = outtab, "data" = out_data, "confint" = Eggerslope_CI)
    class(Out_list)<-"egger"
    
    return(Out_list)
  }
  OUT<-multi_return()
  
}