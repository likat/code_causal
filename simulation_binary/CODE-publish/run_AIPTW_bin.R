### function to calculate the doubly robust coefficients from the estimating equations 
###dataInput-the dataset
###propenVarlist-the variables in the propensity model
###treat.varname-names of the treatment 
###outcome.varname-name of the outcome variable
##outcomeVarList0-variables included in the outcome model for Y0
##outcomeVarList1-varialbes included in the outcome model for Y1
##original-1 for calculating treatment effect estimate on the original dataset, instead of bootstrap samples (for calculating standard error)
###numCarlo-number of monte carlo runs in g computation and estimating the augmentation term in AIPTW
###baseList-variables at baseline, draw from baseline empirical distributions (in the simulations there are 3 baseline covariates for one time point treatment)

# source("../Simulation/CODE/gcomputeFunc_application.R")
AIPTW = 
  function(
    dataInput, 
    sampwts,
    # baselineVar, 
    # linear, 
    # propenVarList, 
    # outcomeVarList0, 
    # outcomeVarList1, 
    covariateXnames,
    treat.varname=trt.varname, 
    outcome.varname=y_varname, 
    propen.model2a = ptrtmod_aiptw,
    modely1 = outcome_mod_aiptw,
    modely0 = outcome_mod_aiptw,
    # data = dat,
    numCarlo=2000) { 
  
  # tryCatch ( 
  #   {
      
      data=NULL
      # if (original==1){
        # 
        data=dataInput
        
      # } else {
        
        ###stratified bootstrap based on treatment groups
        treatID=which(dataInput[, treat.varname]==1)
        controlID=which(dataInput[, treat.varname]==0)
        # data = dataInput[c( sample(treatID,replace=T), sample(controlID,replace=T) ),]  ##stratified bootstraps
        
      # }
      
      treatInd = data[, treat.varname]  ###indicator of treated subjects
      Yobs = data[, outcome.varname]
      
      
      ######## MSM weight########################################################### 
      ###propensity model on numerator
      # (1a) Fit null propensity model on the sample data?
      propen.model1a=formulaF(varList=c(1), y.name=treat.varname)
      model1a=glm(propen.model1a, data=data, family="binomial") 
      
      # phat(A|X) = A(Phat(A=1)) + (1-A)(1-Phat(A=1))
      temp1a=treatInd * model1a$fitted.values + (1-treatInd) * (1-model1a$fitted.values)
      
      ###propensity model on denominator
      # (1b) Then fit full propensity model on the sample data
      # propen.model2a=formulaF(varList=propenVarList, y.name=treat.varname)
      model2a=glm(propen.model2a, data=data, family="binomial") 
      temp2a=treatInd * model2a$fitted.values + (1-treatInd) * (1-model2a$fitted.values)
      
      weight=temp1a/temp2a * sampwts
      
      # summary(weight)
      
      #### use g-computation to compute beta(Qn) 
      ### (2) Prediction model for Y1 and Y0
      # modely0 = formulaF(varList=outcomeVarList0, y.name=outcome.varname)
      # modely1 = formulaF(varList=outcomeVarList1, y.name=outcome.varname)
      
      mod3_1=glm(modely1, data=data[treatInd==1,] %>% as.data.frame(), family = "binomial") 
      mod3_0=glm(modely0, data=data[treatInd==0,]%>% as.data.frame(), family = "binomial") 
      mod3_des = glm(modely1, data=data%>% as.data.frame(), family = "binomial")
      ## numCarlo number of runs in G-computation 
      ############################################################### 
      ###simulate potential outcomes under each treatment A0=1
      ## What is baselineVar? Baseline covariates.. baseline X?
      ## Basically the empirical distribution of X in this case is drawn from a N(mean(X), sd(X)),
      ## each of the baseline covariates are drawn independently... not sure if that's reasonable?
      baselineVar <- data %>% select(paste0(covariateXnames))
      Gcomp1=gcomputeFunc(final.mod=mod3_1, data=data, des = mod3_des,
                          baselineVar=baselineVar, numCarlo=numCarlo, 
                          treat.varname=treat.varname, outcome.varname=outcome.varname, 
                          treat=1)
      # Notes: 
        # This function outputs all the MC baseline covariates, treatment group, and MC Y draws 
      
      ############################################################### 
      ###simulate potential outcomes under each treatment A0=0
      Gcomp0=gcomputeFunc(final.mod=mod3_0, data=data,  des = mod3_des,
                          baselineVar=baselineVar, numCarlo=numCarlo, 
                          treat.varname=treat.varname, outcome.varname=outcome.varname, 
                          treat=0)
      
      ## Bind it all together!
      Gcomp=rbind(Gcomp1, Gcomp0) 
      # dim(Gcomp) 
      
      ##get the treatment coefficients in the MSM model #### 
      designMat=cbind(Gcomp[, which(names(Gcomp)==treat.varname)])
      designY=Gcomp[, which(names(Gcomp)==outcome.varname)]
      
      # (3) Fit MSM on MC Y and the treatment group using the merged imputed items
      msmModel=glm(designY ~ designMat, data=Gcomp,family="binomial")   
      # summary(msmModel)
      
      
      d=matrix(0, nrow=dim(dataInput)[1], ncol=2) 
    
  #### THERE ARE 2 LOOPS:
      # Loop 1 estimates E(D_IPTW|A, X)
      # Loop 2 estimates E(D_IPTW|X)
      # After the loops, you use NR algorithm to find the treatment effect
  for(u in 1:dim(dataInput)[1]) { 
        ## for each row in the simulated data
        
        currentData=NULL ##store baseline coariate
        ### generate L0  from empirical distribution of L0 obtained from the data
        for(ind in 1:ncol(baselineVar)){
          ## for each variable in baseline variables
          draw_g=NULL
          ## define draw_g for observation u as matrix of baseline variable "ind" repeated numCarlo (2000) times
          draw_g=cbind(draw_g, rep(data[u, covariateXnames[ind] ], numCarlo))
          # colnames(draw_g)=baselineVar[ind]
          colnames(draw_g) = covariateXnames[ind]
          currentData=cbind(currentData, draw_g)
          ## Result: if baseline vars are {X1, X2, X3}, then currentData is numCarlo x numBaselineVar 
          ## matrix with [X1, X2, X3] for u
        }
        
        
        ####squared terms ###############
        
        # currentData[ , "age35-65:wiceligYes"] = currentData[,"35-65"]*currentData[,"Yes"]
        # currentData[ , "age65+:wiceligYes"] = currentData[,"65+"]*currentData[,"Yes"]
        
        ## append the treatment group assigned to u 
        currentData=data.frame(currentData, A0_g=data[u, treat.varname])
        names(currentData)[which(names(currentData) == "A0_g")]=treat.varname
        # currentData$hhsize_scaled <- as.numeric(currentData$hhsize_scaled)
        
        # currentData[ , "L1_sq"] = currentData[, "L1"]^2
        # currentData[ , "L2_sq"] = currentData[, "L2"]^2
        # currentData[ , "L3_sq"] = currentData[, "L3"]^2
        # currentData[, "L1L2"] = currentData[, "L1"] * currentData[, "L2"]
        
        ####potential outcome under control
        # draw y^i_mc
        y_mean0 = as.vector(predict(mod3_0, newdata=currentData))
        # y_draw0=rnorm(numCarlo, mean=y_mean0, sd=summary(mod3_0)$sigma)
        y_draw0=rbinom(n=numCarlo, size=1, prob=expit(y_mean0))
        
        ###potential outcome under treatment
        y_mean1 = as.vector(predict(mod3_1, newdata=currentData))
        # y_draw1=rnorm(numCarlo, mean=y_mean1, sd=summary(mod3_1)$sigma)
        y_draw1=rbinom(n=numCarlo, size=1, prob=expit(y_mean1))
        
        ####
        firstTreat = currentData[, treat.varname]
        y_draw=y_draw0 * (1 - firstTreat) + y_draw1 * firstTreat
        
        ## Append MC y^i_mc to current dataset
        currentData[, outcome.varname]=y_draw
        
        ## obtain 
        temp1a_mc_=predict(model1a, currentData, type="response") 
        temp1a_mc=firstTreat * temp1a_mc_ + (1 - firstTreat) * (1-temp1a_mc_)
        
        temp2a_mc_=predict(model2a, currentData, type="response") 
        temp2a_mc=firstTreat * temp2a_mc_ + (1 - firstTreat)* (1-temp2a_mc_)
        
        weight_mc=temp1a_mc/temp2a_mc* rep(dataInput$wts[u], numCarlo)
        
        temp=weight_mc * ( y_draw - expit(cbind(rep(1, numCarlo), firstTreat) %*% msmModel$coef) ) 
        d[u,]=c( mean(temp), mean(firstTreat * temp) ) 
        
        # print(u)
      } 
      
      ########################################################################################################### 
      ### first expectation given history up to time poin 0 (baseline) 
      d_=matrix(0, nrow=dim(dataInput)[1], ncol=2) 
      
      for (u in 1:dim(dataInput)[1]) { 
        
        currentData=NULL ##store baseline covariate
        ### generate L0 from empirical distribution of L0 obtained from the data
        for(ind in 1:length(baselineVar)){
          draw_g=NULL
          # draw_g=cbind(draw_g, rep(as.numeric(data[u, baselineVar[ind] ]), numCarlo))
          draw_g=cbind(draw_g, rep(data[u, covariateXnames[ind] ], numCarlo))
          colnames(draw_g) = covariateXnames[ind]
          currentData=cbind(currentData, draw_g)
          # colnames(draw_g)=baselineVar[ind]
          # currentData=cbind(currentData, draw_g)
        }
        
        currentData=data.frame(currentData, A0_g=data[u, treat.varname])
        names(currentData)[which(names(currentData) == "A0_g")]=treat.varname
        # currentData$hhsize_scaled <- as.numeric(currentData$hhsize_scaled)
        
        ####squared terms ###############
        # currentData[ , "L1_sq"] = currentData[, "L1"]^2
        # currentData[ , "L2_sq"] = currentData[, "L2"]^2
        # currentData[ , "L3_sq"] = currentData[, "L3"]^2
        # currentData[, "L1L2"] = currentData[, "L1"] * currentData[, "L2"]
        
        prob = predict(model2a, currentData, type = "response")  ###probability of treatment at first time point
        
        ## draw a^mc_i
        A0_draw=rbinom(numCarlo, 1, prob ) ### draw treatment at time point 2 
        currentData[, treat.varname]=A0_draw
        
        ####potential outcome under control
        # y_mean0 = as.vector(predict(mod3_0, newdata=currentData))
        # y_draw0=rnorm(numCarlo, mean=y_mean0, sd=summary(mod3_0)$sigma)
        y_mean0 = as.vector(predict(mod3_0, newdata=currentData))
        y_draw0=rbinom(n=numCarlo, size=1, prob=expit(y_mean0))
        
        ###potential outcome under treatment
        # y_mean1 = as.vector(predict(mod3_1, newdata=currentData))
        # y_draw1=rnorm(numCarlo, mean=y_mean1, sd=summary(mod3_1)$sigma)
        y_mean1 = as.vector(predict(mod3_1, newdata=currentData))
        y_draw1=rbinom(n=numCarlo, size=1, prob=expit(y_mean1))
        
        ####
        firstTreat = currentData[, treat.varname]
        y_draw=y_draw0 * (1 - firstTreat) + y_draw1 * firstTreat
        
        currentData[, outcome.varname]=y_draw
        
        # Obtain potential outcomes for the MC values
        temp1a_mc_=predict(model1a, currentData, type="response") 
        temp1a_mc=firstTreat * temp1a_mc_ + (1 - firstTreat) * (1-temp1a_mc_)
        
        temp2a_mc_=predict(model2a, currentData, type="response") 
        temp2a_mc=firstTreat * temp2a_mc_ + (1 - firstTreat)* (1-temp2a_mc_)
        
        ## This is the h(A_i)/P(A_i|X_i) = P(A_i)/P(A_i|X_i)
        weight_mc=temp1a_mc/temp2a_mc * rep(dataInput$wts[u], numCarlo)
        
        ## Calculate D^i_mc
        temp=weight_mc * ( y_draw - expit(cbind(rep(1, numCarlo), firstTreat) %*% msmModel$coef) ) 
        
        ## Record point estimate for D^i_mc overall and for treated only
        d_[u,]=c( mean(temp), mean(firstTreat * temp) ) 
        
        # print(u)
      } 
      
      ######### Newton-Raphson to solve the estimating equations 
      # install.packages("rootSolve")
      # library("rootSolve")
      
      firstTreat = NULL
      ## Treatment group membership
      firstTreat = data[, treat.varname]
      ## m = [1, A]
      sampleSize = nrow(data)
      m=cbind(rep(1,sampleSize), firstTreat)
      ## m2 = pihat = E(D_IPTW|A, X) - E(D_IPTW|X)
      ## *****not sure why we need this extra matrix operation here
      m2= t(d-d_)%*% rep(1,sampleSize)
      
      ## Function for NR algorithm to solve
      D<-function(beta){
        ## m = [1, A] %*% [b0, b1]
        m<-expit(cbind( rep(1,sampleSize), firstTreat ) %*% beta)
        return(
          ## [1,A]^T %*% (P(A)/P(A|X) * (Y-b0-b1A)) - pihat
          ## Basically computing D_i^mc?
          ## Why is cbind(1, A) multiplying though? seems like a duplicate
          t(cbind(1, firstTreat)) %*% (weight*(data[, outcome.varname] - m))-m2
          # (weight*(data[, outcome.varname] - m))-m2
          
        )
      }
      
      oldmybeta<-c(0,1)
      
      repeat{
        mybeta <- oldmybeta - solve( gradient( D, oldmybeta) ) %*% D(oldmybeta)
        if (max(abs(oldmybeta-mybeta))<1e-8) break
        oldmybeta<-mybeta
      }
      
      ## measure of PATE
      result=expit(c(1,1)%*%mybeta)-expit(c(1,0)%*%mybeta)
      
      return(result) 
      
    # }, error=function(e) return(NA) )
  
} 

