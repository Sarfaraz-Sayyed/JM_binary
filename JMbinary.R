
EM_MI_Anal = function(n1 = ,
                      n2 = ,
                      MAR = 10, 
                      MNAR = 0,
                      CMP = 15,
                      rho = 0.8,
                      rhorint = 0.4,
                      seed_addon = 1,
                      comptyp = "C",
                      imputtyp = 'EM',
                      resp_type= 'C',
                      prob1=0.8,
                      prob2=0.5,
                      threshold = 0.7)

{
  #Test drug
  test1 <- Sim_Mis_Compl(n = n1, MARper = MAR, MNARper = MNAR,Compper = CMP,
                         Median_compl = 0.9, perc_10_compl = 0.3, comp_type= comptyp,
                         pos=prob1,
                         corr=rho,
                         seed = seed_addon,
                         cseed = seed_addon)
  
  #Comparator
  test2 <- Sim_Mis_Compl(n = n2, MARper = MAR, MNARper = MNAR,Compper = CMP,
                         Median_compl = 0.8, perc_10_compl = 0.3,comp_type= comptyp,
                         pos=prob2,
                         corr=rho,
                         seed = seed_addon,
                         cseed =seed_addon)
  
  test1cmp <- subset(test1[[2]], select = c(C1, C2, C3, C4, C5, C6) ) 
  test1cmpB <- subset(test1[[2]], select = c(C11, C12, C13, C14, C15, C16) )
  test1rsp <- subset(test1[[2]], select = c(X1, X2, X3, X4, X5, X6) ) 
  imp.cmpdraw1 <- EM_out(dat_in = test1cmp,
                         seedin = (1e5*seed_addon)+6)
  imp.cmpbdraw1 <- EM_out(dat_in = test1cmpB,
                          seedin = (1e5*seed_addon)+6)
  imp.resdraw1 <- EM_out(dat_in = test1rsp,
                         seedin = (1e5*seed_addon)+6)
  imp.draw1 <- cbind(imp.resdraw1,imp.cmpdraw1,imp.cmpbdraw1) 
  
  test2cmp <- subset(test2[[2]], select = c(C1, C2, C3, C4, C5, C6) ) 
  test2cmpB <- subset(test2[[2]], select = c(C11, C12, C13, C14, C15, C16) )
  test2rsp <- subset(test2[[2]], select = c(X1, X2, X3, X4, X5, X6) ) 
  imp.cmpdraw2 <- EM_out(dat_in = test2cmp,
                         seedin = (1e5*seed_addon)+7)
  imp.cmpbdraw2 <- EM_out(dat_in = test2cmpB,
                          seedin = (1e5*seed_addon)+6)
  imp.resdraw2 <- EM_out(dat_in = test2rsp,
                         seedin = (1e5*seed_addon)+7)
  imp.draw2 <- cbind(imp.resdraw2,imp.cmpdraw2,imp.cmpbdraw2) 
  
  imp.draw <- rbind(cbind(data.frame(imp.draw1),"Drug" = rep(1, n1)), cbind(data.frame(imp.draw2),"Drug" = rep(0, n2)))
  imp.draw$subject <- seq.int(sum(n1, n2))
  
  imp.draw_1 <- rbind(cbind(data.frame(test1[[1]]),"Drug" = rep(1, n1)), cbind(data.frame(test2[[1]]),"Drug" = rep(0, n2)))
  imp.draw_1$subject <- seq.int(sum(n1, n2))
  
  imp.draw_2 <- rbind(cbind(data.frame(test1[[2]]),"Drug" = rep(1, n1)), cbind(data.frame(test2[[2]]),"Drug" = rep(0, n2)))
  imp.draw_2$subject <- seq.int(sum(n1, n2))
  
  tresp <- preanal_prep(datain = imp.draw,
                        idvars = c("subject", "Drug"),
                        measurevar = c("X1", "X2", "X3", "X4","X5", "X6"),
                        compliancevars = c("C1", "C2", "C3", "C4", "C5", "C6"),
                        compliancebvars= c("C11", "C12", "C13", "C14", "C15", "C16"),
                        timepoints =c(1, 2, 3, 4, 5, 6))
  tresp1 <- preanal_prep(datain = imp.draw_1,
                         idvars = c("subject", "Drug"),
                         measurevar = c("X1", "X2", "X3", "X4","X5", "X6"),
                         compliancevars = c("C1", "C2", "C3", "C4", "C5", "C6"),
                         compliancebvars= c("C11", "C12", "C13", "C14", "C15", "C16"),
                         timepoints =c(1, 2, 3, 4, 5, 6))
  tresp2 <- preanal_prep(datain = imp.draw_2,
                         idvars = c("subject", "Drug"),
                         measurevar = c("X1", "X2", "X3", "X4","X5", "X6"),
                         compliancevars = c("C1", "C2", "C3", "C4", "C5", "C6"),
                         compliancebvars= c("C11", "C12", "C13", "C14", "C15", "C16"),
                         timepoints =c(1, 2, 3, 4, 5, 6))
  
  max_iter <- 100
  tolerance <- 1e-6
  mu <- c(0, 0)
  Sigmaint <- matrix(c(1, rhorint, rhorint, 1), ncol = 2, byrow = TRUE)
  z1z2 <- MASS::mvrnorm(n = 1200, mu = rep(0, length(mu)), Sigma = Sigmaint)
  tresp2$z1 <- z1z2[,1]
  tresp2$z2 <- z1z2[,2]
  utau <- 0.3
  
  
  # Initialize parameters
  beta <- rep(0.001, 5)  # Intercept and time coefficient for binary response
  alpha <- rep(0.001, 4)
  gamma <- 0  # Coefficient for continuous covariate
  theta <- rep(0.001, n1)  # Random effects for binary response
  eta <- rep(0.001, n2)  # Random effects for continuous covariate
  
  # tresp2 <- within(tresp2, Drug <- relevel(as.factor(Drug), ref = 0))
  
  for (iter in 1:max_iter) {
    
    if (comptyp == "C")
    {
      # E-step: Impute missing continuous covariate using current parameter estimates
      tresp2$expected_covariate <- alpha[1]+alpha[2]*as.numeric(tresp2$Drug)+alpha[3]*tresp2$Time + eta + alpha[4]*tresp2$z2
      tresp2$imputed_covariate <- ifelse(is.na(tresp2$compliance), tresp2$expected_covariate, tresp2$compliance)
      
      # M-step: Update parameter estimates for continuous covariate using imputed data
      old_alpha <- alpha
      old_eta <- eta
      
      # Fit a linear regression for normal using imputed data
      # model_normal <- lm(tresp2$imputed_covariate ~ tresp2$Drug + tresp2$Time + tresp2$z2)
      # alpha <- coef(model_normal)
      model_normal <- lme4::lmer(
        imputed_covariate ~ Drug + Time + z2 + (1|subject),
        data = tresp2
      )
      alpha <- summary(model_normal)$coefficients[,1]
      
      # Update random effects for complaince using imputed data
      eta <- mean(tresp2$imputed_covariate - alpha[1] - alpha[2] * as.numeric(tresp2$Drug) - alpha[3] * tresp2$Time - alpha[4]*tresp2$z2)
      tresp2$resid <- tresp2$imputed_covariate - alpha[1] - alpha[2] * as.numeric(tresp2$Drug) - alpha[3] * tresp2$Time  - alpha[4]*tresp2$z2                 
    } else
    { 
      # E-step: Impute missing binary covariate using current parameter estimates
      tresp2$expected_covariate <- exp(alpha[1]+alpha[2]*as.numeric(tresp2$Drug)+alpha[3]*tresp2$Time + eta + alpha[4]*tresp2$z2)/
        (1+exp(alpha[1]+alpha[2]*as.numeric(tresp2$Drug)+alpha[3]*tresp2$Time + eta + alpha[4]*tresp2$z2))
      tresp2$imputed_covariate <- ifelse(is.na(tresp2$compliance), tresp2$expected_covariate, tresp2$compliance)
      
      # M-step: Update parameter estimates for binary covariate using imputed data
      old_alpha <- alpha
      old_eta <- eta
      
      # Fit a logistic regression for beta using imputed data
      # model_beta <- glm(tresp2$imputed_response ~ tresp2$Drug + tresp2$Time + tresp2$resid + tresp2$U, family = "binomial")
      modelc_beta <- lme4::glmer(
        imputed_covariate ~ Drug + Time + z2 + (1|subject),
        data = tresp2,
        family = binomial(link = "logit")
      )
      # beta <- coef(model_beta)
      alpha <- summary(modelc_beta)$coefficients[,1]
      # Update random effects for binary response using imputed data
      tresp2$pred <- exp(alpha[1] + alpha[2] * as.numeric(tresp2$Drug) + alpha[3] * tresp2$Time + alpha[4]*tresp2$z2)
      eta <- mean(tresp2$imputed_covariate - tresp2$pred/(1+tresp2$pred))
      tresp2$resid <- tresp2$imputed_covariate - tresp2$pred/(1+tresp2$pred) 
    }
    
    phi <- 1.0/sqrt(1+utau*utau*(3/pi*pi))
    tresp2$U <- (1/phi)*log (sin (pi*pnorm(tresp2$z1)*phi)/sin (phi*pi*(1-pnorm(tresp2$z1))))
    
    # E-step: Impute missing binary response using current parameter estimates
    tresp2$expx <- exp(beta[1] + beta[2]*as.numeric(tresp2$Drug) + beta[3]*tresp2$Time  + beta[4]*tresp2$resid + beta[5]*tresp2$U + theta)
    tresp2$expected_response <- tresp2$expx/(1+tresp2$expx)
    tresp2$imputed_response <- ifelse(is.na(tresp2$respval), tresp2$expected_response, tresp2$respval)
    tresp2$imputed_response <- ifelse(tresp2$imputed_response < 0 , 0, tresp2$imputed_response)
    tresp2$imputed_response <- ifelse(tresp2$imputed_response > 1 , 1, tresp2$imputed_response)
    
    # M-step: Update parameter estimates for binary response using imputed data
    old_beta <- beta
    old_theta <- theta
    
    # Fit a logistic regression for beta using imputed data
    # model_beta <- glm(tresp2$imputed_response ~ tresp2$Drug + tresp2$Time + tresp2$resid + tresp2$U, family = "binomial")
    model_beta <- lme4::glmer(
      imputed_response ~ Drug + Time + resid + U + (1|subject),
      data = tresp2,
      family = binomial(link = "logit")
    )
    # beta <- coef(model_beta)
    beta <- summary(model_beta)$coefficients[,1]
    # Update random effects for binary response using imputed data
    tresp2$expx1 <- exp(beta[1] + beta[2] * as.numeric(tresp2$Drug) + beta[3] * tresp2$Time + beta[4]*tresp2$resid + beta[5]*tresp2$U)
    theta <- mean(tresp2$imputed_response - tresp2$expx1/(1+tresp2$expx1))
    
    # Check for convergence
    if (max(abs(beta - old_beta)) < tolerance && max(abs(theta - old_theta)) < tolerance && 
        max(abs(alpha - old_alpha)) < tolerance && max(abs(eta - old_eta)) < tolerance ) {
      break
    }
    # return(list(beta, alpha,theta, eta))
  }
  # glm(tresp1$respval ~ tresp2$Drug + tresp2$Time , family = "binomial")
  # glm(tresp2$imputed_response ~ tresp2$Drug + tresp2$Time , family = "binomial")
  # glm(ifelse(is.na(tresp2$respval),0,tresp2$respval) ~ tresp2$Drug + tresp2$Time , family = "binomial")
  
  compl <- lme4::glmer(
    respval  ~ Drug + Time + (1|subject),
    data = tresp1,
    family = binomial(link = "logit")
  )
  tresp1$respval_pred = ifelse (exp(predict(compl))/(1+exp(predict(compl))) > threshold, 1, 0)
  tresp1$classif <- dplyr::case_when(tresp1$respval == 1 & tresp1$respval_pred == 1 ~ 'TP',
                                     tresp1$respval == 0 & tresp1$respval_pred == 1 ~ 'FP',
                                     tresp1$respval == 0 & tresp1$respval_pred == 0 ~ 'TN',
                                     tresp1$respval == 1 & tresp1$respval_pred == 0 ~ 'FN'
  )
  comp_confusion.matrix <- plyr::count(tresp1, 'classif')
  compl.precision <- comp_confusion.matrix[4,2]/(comp_confusion.matrix[4,2]+comp_confusion.matrix[2,2])
  compl.recall <- comp_confusion.matrix[4,2]/(comp_confusion.matrix[4,2]+comp_confusion.matrix[1,2])
  compl.f1score <- (2*compl.precision*compl.recall)/(compl.precision + compl.recall)
  compl.TPR <- comp_confusion.matrix[4,2]/(comp_confusion.matrix[4,2]+comp_confusion.matrix[1,2])
  compl.FPR <- comp_confusion.matrix[2,2]/(comp_confusion.matrix[2,2]+comp_confusion.matrix[3,2])
  compl.ones <- tresp1[tresp1$respval==1, ] # Subset ones
  compl.zeros <- tresp1[tresp1$respval==0, ] # Subsetzeros
  compl.totalPairs <- nrow (compl.ones) * nrow (compl.zeros) # calculate total number of pairs to check
  compl.conc <- sum(c(vapply(compl.ones$respval_pred, function(x) {((x > compl.zeros$respval_pred))}, FUN.VALUE=logical(nrow(compl.zeros)))), na.rm=T)
  compl.disc <- sum(c(vapply(compl.ones$respval_pred, function(x) {((x < compl.zeros$respval_pred))}, FUN.VALUE = logical(nrow(compl.zeros)))), na.rm = T)
  compl.concordance <- compl.conc/compl.totalPairs
  compl.discordance <- compl.disc/compl.totalPairs
  compl.tiesPercent <- (1-compl.concordance-compl.discordance)
  compl.AUC = compl.concordance + 0.5*compl.tiesPercent
  # pROC::auc(tresp1$respval, tresp1$respval_pred)
  compl.AUC_lower <- pROC::ci.auc(tresp1$respval, tresp1$respval_pred, cinf.level = 0.95, method=c("delong"))[1]
  compl.AUC_median <- pROC::ci.auc(tresp1$respval, tresp1$respval_pred, cinf.level = 0.95, method=c("delong"))[2]
  compl.AUC_upper <- pROC::ci.auc(tresp1$respval, tresp1$respval_pred, cinf.level = 0.95, method=c("delong"))[3]
  # pROC::ci.auc(tresp1$respval, tresp1$respval_pred, cinf.level = 0.95, method=c("bootstrap"), 
  #              boot.n = 2000)[1]
  # pROC::ci.auc(tresp1$respval, tresp1$respval_pred, cinf.level = 0.95, method=c("bootstrap"), 
  #              boot.n = 2000)[3]
  
  
  imput <- lme4::glmer(
    imputed_response  ~ Drug + Time + (1|subject),
    data = tresp2,
    family = binomial(link = "logit")
  )
  tresp2$respval_pred = ifelse (exp(predict(imput))/(1+exp(predict(imput))) > threshold, 1, 0)
  tresp2$respval_o <- tresp1$respval
  tresp2$classif <- dplyr::case_when(tresp2$respval_o == 1 & tresp2$respval_pred == 1 ~ 'TP',
                                     tresp2$respval_o == 0 & tresp2$respval_pred == 1 ~ 'FP',
                                     tresp2$respval_o == 0 & tresp2$respval_pred == 0 ~ 'TN',
                                     tresp2$respval_o == 1 & tresp2$respval_pred == 0 ~ 'FN'
  )
  imput_confusion.matrix <- plyr::count(tresp2, 'classif')
  imput.precision <- imput_confusion.matrix[4,2]/(imput_confusion.matrix[4,2]+imput_confusion.matrix[2,2])
  imput.recall <- imput_confusion.matrix[4,2]/(imput_confusion.matrix[4,2]+imput_confusion.matrix[1,2])
  imput.f1score <- (2*imput.precision*imput.recall)/(imput.precision + imput.recall)
  imput.TPR <- imput_confusion.matrix[4,2]/(imput_confusion.matrix[4,2]+imput_confusion.matrix[1,2])
  imput.FPR <- imput_confusion.matrix[2,2]/(imput_confusion.matrix[2,2]+imput_confusion.matrix[3,2])
  imput.ones <- tresp2[tresp2$respval_o==1, ]  # Subset ones
  imput.zeros <- tresp2[tresp2$respval_o==0, ]  # Subsetzeros
  imput.totalPairs <- nrow (imput.ones) * nrow (imput.zeros) # calculate total number of pairs to check
  imput.conc <- sum(c(vapply(imput.ones$respval_pred, function(x) {((x > imput.zeros$respval_pred))}, FUN.VALUE=logical(nrow(imput.zeros)))), na.rm=T)
  imput.disc <- sum(c(vapply(imput.ones$respval_pred, function(x) {((x < imput.zeros$respval_pred))}, FUN.VALUE = logical(nrow(imput.zeros)))), na.rm = T)
  imput.concordance <- imput.conc/imput.totalPairs
  imput.discordance <- imput.disc/imput.totalPairs
  imput.tiesPercent <- (1-imput.concordance-imput.discordance)
  imput.AUC = imput.concordance + 0.5*imput.tiesPercent
  imput.AUC_lower <- pROC::ci.auc(tresp2$respval_o, tresp1$respval_pred, cinf.level = 0.95, method=c("delong"))[1]
  imput.AUC_median <- pROC::ci.auc(tresp2$respval_o, tresp1$respval_pred, cinf.level = 0.95, method=c("delong"))[2]
  imput.AUC_upper <- pROC::ci.auc(tresp2$respval_o, tresp1$respval_pred, cinf.level = 0.95, method=c("delong"))[3]
  
  missn <- lme4::glmer(
    ifelse(is.na(tresp2$respval),0,tresp2$respval)  ~ Drug + Time + (1|subject),
    data = tresp2,
    family = binomial(link = "logit")
  )
  tresp2$respval_pred1 = ifelse (exp(predict(missn))/(1+exp(predict(missn))) > threshold, 1, 0)
  tresp2$classif1 <- dplyr::case_when(tresp2$respval_o == 1 & tresp2$respval_pred1 == 1 ~ 'TP',
                                      tresp2$respval_o == 0 & tresp2$respval_pred1 == 1 ~ 'FP',
                                      tresp2$respval_o == 0 & tresp2$respval_pred1 == 0 ~ 'TN',
                                      tresp2$respval_o == 1 & tresp2$respval_pred1 == 0 ~ 'FN'
  )
  missn_confusion.matrix <- plyr::count(tresp2, 'classif1')
  missn.precision <- missn_confusion.matrix[4,2]/(missn_confusion.matrix[4,2]+missn_confusion.matrix[2,2])
  missn.recall <- missn_confusion.matrix[4,2]/(missn_confusion.matrix[4,2]+missn_confusion.matrix[1,2])
  missn.f1score <- (2*missn.precision*missn.recall)/(missn.precision + missn.recall)
  missn.TPR <- missn_confusion.matrix[4,2]/(missn_confusion.matrix[4,2]+missn_confusion.matrix[1,2])
  missn.FPR <- missn_confusion.matrix[2,2]/(missn_confusion.matrix[2,2]+missn_confusion.matrix[3,2])
  missn.ones <- tresp2[tresp2$respval_o==1, ]  # Subset ones
  missn.zeros <- tresp2[tresp2$respval_o==0, ]  # Subsetzeros
  missn.totalPairs <- nrow (missn.ones) * nrow (missn.zeros) # calculate total number of pairs to check
  missn.conc <- sum(c(vapply(missn.ones$respval_pred1, function(x) {((x > missn.zeros$respval_pred1))}, FUN.VALUE=logical(nrow(missn.zeros)))), na.rm=T)
  missn.disc <- sum(c(vapply(missn.ones$respval_pred1, function(x) {((x < missn.zeros$respval_pred1))}, FUN.VALUE = logical(nrow(missn.zeros)))), na.rm = T)
  missn.concordance <- missn.conc/missn.totalPairs
  missn.discordance <- missn.disc/missn.totalPairs
  missn.tiesPercent <- (1-missn.concordance-missn.discordance)
  missn.AUC = missn.concordance + 0.5*missn.tiesPercent
  missn.AUC_lower <- pROC::ci.auc(tresp2$respval_o, tresp2$respval_pred1, cinf.level = 0.95, method=c("delong"))[1]
  missn.AUC_median <- pROC::ci.auc(tresp2$respval_o, tresp2$respval_pred1, cinf.level = 0.95, method=c("delong"))[2]
  missn.AUC_upper <- pROC::ci.auc(tresp2$respval_o, tresp2$respval_pred1, cinf.level = 0.95, method=c("delong"))[3]
  
  parms_mis <- summary(missn)$coefficients[,]
  deviance_mis <- summary(missn)$AIC
  
  parms_comp <- summary(compl)$coefficients[,]
  deviance_comp <- summary(compl)$AIC
  
  parms_imp <- summary(imput)$coefficients[,]
  deviance_imp <- summary(imput)$AIC
  
  parms <- data.frame(rbind(cbind(parms_comp,categ='com',f1score = compl.f1score,
                                  TP = comp_confusion.matrix[4,2], FP = comp_confusion.matrix[2,2],
                                  TN = comp_confusion.matrix[3,2], FN = comp_confusion.matrix[1,2],
                                  TPR = compl.TPR, FPR = compl.FPR, AUC = compl.AUC, AUC_median = compl.AUC_median,
                                  AUC_lower = compl.AUC_lower,AUC_upper = compl.AUC_upper), 
                            cbind(parms_mis,categ='mis',f1score=missn.f1score,
                                  TP = missn_confusion.matrix[4,2], FP = missn_confusion.matrix[2,2],
                                  TN = missn_confusion.matrix[3,2], FN = missn_confusion.matrix[1,2],
                                  TPR = missn.TPR,FPR = missn.FPR, AUC = missn.AUC, AUC_median = missn.AUC_median,
                                  AUC_lower = missn.AUC_lower,AUC_upper = missn.AUC_upper), 
                            cbind(parms_imp,categ='imp',f1score=imput.f1score,
                                  TP = imput_confusion.matrix[4,2], FP = imput_confusion.matrix[2,2],
                                  TN = imput_confusion.matrix[3,2], FN = imput_confusion.matrix[1,2],
                                  TPR = imput.TPR,FPR = imput.FPR,AUC = imput.AUC, AUC_median = imput.AUC_median,
                                  AUC_lower = imput.AUC_lower,AUC_upper = imput.AUC_upper)))
  deviance <- data.frame(rbind(cbind(deviance_comp,categ='com'), cbind(deviance_mis,categ='mis'), cbind(deviance_imp,categ='imp')))
  
  return(list(parms, deviance))
} 

