
Sim_Mis_Compl = function(
  MARper = 10,
  MNARper = 5,
  Compper = 15,
  comp_type= 'C',
  Median_compl = 0.9,
  perc_10_compl = 0.3,
  comp_cov = ,
  comp_mean = ,
  mu = c(0, 0, 0, 0, 0, 0),
  sigma = cov,
  pos = pos,
  ncol = length(mu),
  byrow = TRUE,
  corr = 0.8,
  n = 100,
  seed = 1e5,
  cseed = 1e6
)

{
  
  #Compliance 
  sigmacomp <- comp_cov
    
  bg <-
    (qnorm(perc_10_compl) - qnorm(Median_compl)) / qnorm(perc_10_compl)
  ag <- qnorm(Median_compl)
  set.seed(cseed)
  Cijt_ <-
    ag + bg*(MASS::mvrnorm(
      n,
      comp_mean,
      sigmacomp,
      tol = 1e-6,
      empirical = FALSE,
      EISPACK = FALSE
    ))
  Cijt <- data.frame(round(pnorm(Cijt_), 3))
  colnames(Cijt) <- paste0("C",seq.int(1,length(mu),1)) 
  
  set.seed(seed)
  x0 <-
    MASS::mvrnorm(
      n,
      mu,
      cov,
      tol = 1e-6,
      empirical = FALSE,
      EISPACK = FALSE
    )
  
  for (i in 2:length(mu))
  {
    x0[, i] <- x0[, 1] + Cijt[,i] * (x0[, i] - x0[, 1])
  }
  
  if (MNARper > 0)
  {
    rslt_MNAR <-
      mice::ampute(
        x0,
        prop = MNARper / 100,
        mech = "MNAR",
        freq = my_freq,
        patterns = my_pat,
        bycases = TRUE
      )
    rsltMNR <- rslt_MNAR$amp
  } else
  {
    rsltMNR <- data.frame(x0)
  }
  rslt_MAR <-
    mice::ampute(
      rsltMNR[stats::complete.cases(rsltMNR[, paste0("X",seq.int(1,length(mu),1))]), ],
      prop = MARper, 
      mech = "MAR",
      freq = my_freq,
      patterns = my_pat,
      bycases = TRUE
    )
  
  rslt <- rbind(rslt_MAR$amp,
                rsltMNR[!stats::complete.cases(rsltMNR), ])
  rslt_missid <- which(is.na(rslt))
  brslt<- x0[NA,]
  
  for (i in 1:length(mu))
  {
    brslt[,i] <- ifelse(x0[,i] >= quantile(x0[,i],1-pos),1,0)
  }
  
  brslt_miss <- brslt
  brslt_miss[rslt_missid] <- NA
  brslt_miss <- data.frame(brslt_miss)
  
  rslt_CMP <-
    mice::ampute(
      Cijt,
      prop = Compper / 100,
      mech = "MAR",
      freq = my_freq,
      patterns = my_pat,
      bycases = TRUE
    )
  
  rsltCMP0 <- rslt_CMP$amp
  rsltCMP1 <- data.frame(matrix(as.numeric(rsltCMP0>0.7),ncol=length(mu)))
  rsltCMP2 <- data.frame(matrix(as.numeric(Cijt>0.7),ncol=length(mu)))
  
  colnames(rsltCMP0) <- paste0("C",seq.int(1,length(mu),1)) 
  colnames(rsltCMP1) <- paste0("C1",seq.int(1,length(mu),1)) 
  colnames(rsltCMP2) <- paste0("C1",seq.int(1,length(mu),1)) 
  
  dat_with_miss <-  cbind(brslt_miss, rsltCMP0,rsltCMP1)
  colnames(brslt) <- paste0("X",seq.int(1,length(mu),1)) 
  Original_data <- cbind(brslt,Cijt,rsltCMP2)
  
  return(list(Original_data,dat_with_miss))
}
