#' ivw_radial
#'
#' Fits a radial inverse variance weighted (IVW) model using first order, second order, or modified second order weights. Outliers are identified using a significance threshold specified by the user. The function returns an object of class \code{"IVW"}, containing regression estimates, a measure of total heterogeneity using Cochran's Q statistic, the individual contribution to overall heterogeneity of each variant, and a data frame for use in constructing the radial plot. 
#' 
#' @param r_input A formatted data frame using the \code{format_radial} function.
#' @param alpha A value specifying the statistical significance threshold for identifying outliers (\code{0.05} specifies a p-value threshold of 0.05).
#' @param weights A value specifying the inverse variance weights used to calculate IVW estimate and Cochran's Q statistic. By default modified second order weights are used, but one can choose to select first order (\code{1}), second order (\code{2}) or modified second order weights (\code{3}).
#' @param tol A value indicating the tolerance threshold for performing the iterative IVW approach. The value represents the minimum difference between the coefficients of the previous and current iterations required for a further iteration to be performed (default= \code{0.0001}).
#' @return An object of class \code{"IVW"} containing the following components:\describe{
#' \item{\code{coef}}{The estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{qstatistic}}{Cochran's Q statistic for overall heterogeneity.}
#' \item{\code{df}}{Degrees of freedom. This is equal to the number of variants -1 when fitting the radial IVW model.}
#' \item{\code{outliers}}{A data frame containing variants identified as outliers, with respective Q statistics, chi-squared tests and SNP identification.}
#' \item{\code{data}}{A data frame containing containing SNP IDs, inverse variance weights, the product of the inverse variance weight and ratio estimate for each variant, contribution to overall heterogeneity with corresponding p-value, and a factor indicator showing outlier status.}
#' \item{\code{confint}}{A vector giving lower and upper confidence limits for the radial IVW effect estimate.}
#' \item{\code{it.coef}}{The estimated iterative coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{it.confint}}{A vector giving lower and upper confidence limits for the iterative radial IVW effect estimate.} 
#' \item{\code{fe.coef}}{The estimated fixed effect exact coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{fe.confint}}{A vector giving lower and upper confidence limits for the fixed effect exact radial IVW effect estimate.}
#' \item{\code{re.coef}}{The estimated random effect exact coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
#' \item{\code{re.confint}}{A vector giving lower and upper confidence limits for the random effect exact radial IVW effect estimate.}
#' \item{\code{mf}}{The mean F statistic for the set of genetic variants, indicative of instrument strength.}
#' 
#'}
#'@author Wes Spiller; Jack Bowden.
#'@references Bowden, J., et al., Improving the visualization, interpretation and analysis of two-sample summary data Mendelian randomization via the Radial plot and Radial regression. International Journal of Epidemiology, 2018. 47(4): p. 1264-1278.
#'@export
#'@examples
#'
#' ivw_radial(r_input,0.05,1,0.0001,T)
#' 


ivw_radial<-function(r_input,alpha,weights,tol){
  
  #Define ratio estimates
  Ratios<-r_input[,3]/r_input[,2]
  
  F<- r_input[,2]^2/r_input[,4]^2
  mf<- mean(F)
  
  cat()
  
  if(missing(alpha)) {
    alpha<-0.05
    warning("Significance threshold for outlier detection not specified: Adopting a 95% threshold")
  }
  
  if(missing(weights)) {
    weights<-3
    warning("Weights not specified: Adopting modified second-order weights")
  }
  
  if(missing(tol)) {
    tol<-0.00001
  }
  
    summary<-TRUE

  
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
  
  DF<-length(r_input[,1])-1
  
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
  
  Iterative_ivw<-function(int.tol){
    Diff  <- 1
    Bhat1.Iterative <- 0 
    count <- 0
    while(Diff >= tol){
      
      W    <- 1/(r_input[,5]^2/r_input[,2]^2 + (Bhat1.Iterative^2)*r_input[,4]^2/r_input[,2]^2)
      
      #Define vector of squared weights
      Wj<-sqrt(W)
      
      #Define vector of weights * ratio estimates
      BetaWj<-Ratios*Wj
      
      #Define IVW Model
      new.IVW.Model<-lm(BetaWj~-1+Wj)
      
      #Define output of IVW.Model fit
      new.EstimatesIVW<-summary(lm(BetaWj~-1+Wj))
      
      #Define slope parameter for IVW.Model
      new.IVW.Slope<-new.EstimatesIVW$coefficients[1]
      
      #Define standard error for slope parameter of IVW.Model
      new.IVW.SE<-new.EstimatesIVW$coefficients[2]
      
      #Define confidence interval for IVW.Model slope parameter
      new.IVW_CI<-confint(new.IVW.Model)
      
      #Define Q statistic for each individual variant
      new.Qj<-W*(Ratios-new.IVW.Slope)^2
      
      #Define total Q statistic
      new.Total_Q<-sum(new.Qj)
      
      #Perform chi square test for overall Q statistic
      new.Total_Q_chi<-pchisq(new.Total_Q,length(r_input[,2])-1,lower.tail = FALSE)
      
      Diff  <- abs(Bhat1.Iterative - new.IVW.Slope) 
      Bhat1.Iterative <- new.IVW.Slope
      Bhat1.SE<-new.IVW.SE
      #Bhat1.CI<-new.IVW_CI
      Bhat1.t<-summary(new.IVW.Model)$coefficients[1,3]
      Bhat1.p<-summary(new.IVW.Model)$coefficients[1,4]
      count <- count+1
    }
    
    It.Dat<-data.frame(Bhat1.Iterative,Bhat1.SE,Bhat1.t,Bhat1.p)
    
    multi_return2 <- function() {
      Out_list2 <- list("It.Res" = It.Dat,"count"= count,"It.CI"=new.IVW_CI)
      
      return(Out_list2) 
    }
    OUT2<-multi_return2()
    
    
  }
  
  Bhat1.Iterative<-Iterative_ivw(tol)
  
  ###################EXACT WEIGHTS######################
  
  
  
  
  
  
  
  
  
  
  PL2 = function(a){
    b = a[1]
    w = 1/((phi)*r_input[,5]^2/r_input[,2]^2 + (b^2)*r_input[,4]^2/r_input[,2]^2)
    q =  sum(w*(Ratios - b)^2)
  }
  
  PLfunc = function(a){
    phi    = a[1]
    PL2    = function(a){
      beta   = a[1]
      w      = 1/(phi*r_input[,5]^2/r_input[,2]^2 + (beta^2)*r_input[,4]^2/r_input[,2]^2)
      q      =  (sum(w*(Ratios - beta)^2))
    }#
    b  = optimize(PL2,interval=c(lb,ub))$minimum 
    w    = 1/(phi*r_input[,5]^2/r_input[,2]^2 + (b^2)*r_input[,4]^2/r_input[,2]^2)
    q    =  (sum(w*(Ratios - b)^2) - DF)^2
  }
  
  BootVar = function(sims=1000){
    B = NULL ; pp=NULL
    for(hh in 1:sims){
      L      = length(r_input[,2])
      choice = sample(seq(1,L),L,replace=TRUE)
      bxg    = r_input[,2][choice] ; seX = r_input[,4][choice]
      byg    = r_input[,3][choice] ; seY = r_input[,5][choice]
      Ratios    = byg/bxg
      
      W1        = 1/(seY^2/bxg^2)
      BIVw1     = Ratios*sqrt(W1)
      sW1       = sqrt(W1)
      IVWfitR1  = summary(lm(BIVw1 ~ -1+sW1))
      phi_IVW1  = IVWfitR1$sigma^2
      W2        = 1/(seY^2/bxg^2 + (byg^2)*seX^2/bxg^4)
      BIVw2     = Ratios*sqrt(W2)
      sW2       = sqrt(W2)
      IVWfitR2  = summary(lm(BIVw2 ~ -1+sW2))
      phi_IVW2  = IVWfitR2$sigma^2
      
      phi_IVW2 = max(1,phi_IVW2)
      phi_IVW1 = max(1,phi_IVW1) 
      lb       = IVWfitR1$coef[1] - 10*IVWfitR1$coef[2]
      ub       = IVWfitR1$coef[1] + 10*IVWfitR1$coef[2]
      
      PL2 = function(a){
        b = a[1]
        w = 1/((phi)*seY^2/bxg^2 + (b^2)*seX^2/bxg^2)
        q =  sum(w*(Ratios - b)^2)
      }
      
      PLfunc = function(a){
        phi    = a[1]
        PL2    = function(a){
          beta   = a[1]
          w      = 1/(phi*seY^2/bxg^2 + (beta^2)*seX^2/bxg^2)
          q      =  (sum(w*(Ratios - beta)^2))
        }
        b  = optimize(PL2,interval=c(-lb,ub))$minimum 
        w    = 1/(phi*seY^2/bxg^2 + (b^2)*seX^2/bxg^2)
        q    =  (sum(w*(Ratios - b)^2) - DF)^2
      }
      phi    = optimize(PLfunc,interval=c(phi_IVW2,phi_IVW1+0.001))$minimum 
      B[hh]  = optimize(PL2,interval=c(lb,ub))$minimum 
    }
    se   = sd(B)
    mB   = mean(B)
    return(list(mB=mB,se=se))
  }
  
  CIfunc = function(){
    z = qt(df=DF, 0.975)
    z2 = 2*(1-pnorm(z))
    
    PL3 = function(a){
      b = a[1]
      w = 1/(r_input[,5]^2/r_input[,2]^2 + (b^2)*r_input[,4]^2/r_input[,2]^2)
      q    =  (sum(w*(Ratios - b)^2) - qchisq(1-z2,DF))^2
    }
    
    lb = Bhat - 10*SE
    ub = Bhat + 10*SE
    
    low   = optimize(PL3,interval=c(lb,Bhat))$minimum
    high  = optimize(PL3,interval=c(Bhat,ub))$minimum
    CI    = c(low,high)
    return(list(CI=CI))
  }
  
  
  #######################
  # Fixed effect model  #
  # and Exact Q test    #
  #######################
  
  phi       = 1
  Bhat      = optimize(PL2,interval=c(-2,2))$minimum 
  W         = 1/(r_input[,5]^2/r_input[,2]^2 + (Bhat^2)*r_input[,4]^2/r_input[,2]^2)
  SE        = sqrt(1/sum(W))
  FCI        = CIfunc()
  
  # Qtest
  
  QIVW      = sum(W*(Ratios - Bhat)^2)
  Qp        = 1-pchisq(QIVW,DF)
  Qind      = W*(Ratios - Bhat)^2
  ExactQ    = c(QIVW,Qp)
  ExactQind = Qind
  # Estimation (fixed effects)
  # point estimate, se, t-stat, p-value)
  
  FE_EXACT    = t(c(Bhat,SE,Bhat/SE,2*(1-pt(abs(Bhat/SE),DF))))
  
  FE_EXACT<-data.frame(FE_EXACT)
  
  names(FE_EXACT)<-c("Estimate","Std.Error","t value","Pr(>|t|)")
  
  ########################
  # Random effect model  #
  # and Exact Q test     #
  ########################
  
  BIVW1      = Ratios*sqrt(1/(r_input[,5]^2/r_input[,2]^2))
  IVWfit1    = summary(lm(BIVW1 ~ -1+sqrt(1/(r_input[,5]^2/r_input[,2]^2))))
  phi_IVW1   = IVWfit1$sigma^2
  
  
  BIVW2 <- Ratios*sqrt(1/(r_input[,5]^2/r_input[,2]^2 + (r_input[,3]^2)*r_input[,4]^2/r_input[,2]^4))
  IVWfit2    = summary(lm(BIVW2 ~ -1+sqrt(1/(r_input[,5]^2/r_input[,2]^2 + (r_input[,3]^2)*r_input[,4]^2/r_input[,2]^4))))
  
  phi_IVW2   = IVWfit2$sigma^2
  
  phi_IVW2 = max(1,phi_IVW2)
  phi_IVW1 = max(1,phi_IVW1)+0.001
  
  lb = Bhat - 10*SE
  ub = Bhat + 10*SE
  
  phi       = optimize(PLfunc,interval=c(phi_IVW2,phi_IVW1))$minimum 
  Bhat      = optimize(PL2,interval=c(lb,ub))$minimum 
  Boot      = BootVar()
  SE        = Boot$se 
  
  # point estimate, se, t-stat, p-value)
  
  RCI    = Bhat + c(-1,1)*qt(df=DF, 0.975)*SE
  RE_EXACT  = t(c(Bhat,SE,Bhat/SE,2*(1-pt(abs(Bhat/SE),DF))))
  
  RE_EXACT<-data.frame(RE_EXACT)
  names(RE_EXACT)<-c("Estimate","Std.Error","t value","Pr(>|t|)")
  
  
  
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
    
    Sum.Dat<-data.frame(coef(EstimatesIVW))
    names(Sum.Dat)<-c("Estimate","Std.Error","t value","Pr(>|t|)")
    names(Bhat1.Iterative$It.Res)<-names(Sum.Dat)
    combined.dat<-(rbind(Sum.Dat,Bhat1.Iterative$It.Res))
    combined.dat<-rbind(combined.dat,FE_EXACT)
    combined.dat<-rbind(combined.dat,RE_EXACT)
    
    row.names(combined.dat)<-c("Effect","Iterative","Exact (FE)","Exact (RE)")
    
    if(weights == 1){
      
      row.names(combined.dat)[1] <- "Effect (1st)"
      
    }
    
    if(weights == 2){
      
      row.names(combined.dat)[1] <- "Effect (2nd)"
      
    }
    
    if(weights == 3){
      
      row.names(combined.dat)[1] <- "Effect (Mod.2nd)"
      
    }
    
    
    print(combined.dat)
    
    cat("\n")
    
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
    cat("Number of iterations =", Bhat1.Iterative$count)
    cat("\n")
    
  }
  
  out_data<-data.frame(r_input[,1],r_input[,6],r_input[,7],r_input[,8])
  out_data$Wj<-Wj
  out_data$BetaWj<-BetaWj
  out_data<-out_data[c(1,5,6,2,3,4)]
  names(out_data)<-c("SNP","Wj","BetaWj","Qj","Qj_Chi","Outliers")
  
  multi_return <- function() {
    Out_list <- list("coef" = EstimatesIVW$coef,"qstatistic"= Total_Q, "df" = length(r_input[,2])-1, "outliers" = outtab, "data" = out_data, "confint" = confint(IVW.Model),
                     "it.coef"=combined.dat[2,],"fe.coef"=combined.dat[3,],"re.coef" = combined.dat[4,1], "it.confint"= Bhat1.Iterative$It.CI,"fe.confint" = FCI$CI, "re.confint" = RCI, "meanF"= mf)
    class(Out_list)<-"IVW"
    
    return(Out_list) 
  }
  OUT<-multi_return()
  
  
}

