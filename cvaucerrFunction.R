cv.auc_err <- function(jointFit = jointFit1, n = 5, iter = 500L, startTime, PredTime ){
  data <- jointFit$Data[[2]]    ## extract the survival data from jointfit object
  dd <- data[sample(nrow(data), nrow(data), replace = F),]    ## randomly sample the survival data
  dat.splited <- split(dd, rep(1:n, each=nrow(dd)/(n)))   ## assign data to n lists
  auc.vec_jointFit1 <- NULL
  pred.vec_jointFit1 <-NULL
  
  for (i in 1:n) {
    ## creat training and validation dataset for survival process and longitidinal process
    trainSurv <- dat.splited[-i]
    trainSurv <- do.call(rbind, trainSurv)
    trainSurv <- trainSurv[order(trainSurv$id),]
    validSurv <- dat.splited[[i]]
    AoValv <- jointFit$Data[[1]]
    trainLME <- AoValv[AoValv$id %in% trainSurv$id,]
    validLME <- AoValv[AoValv$id %in% validSurv$id,]
    ## fit lme and survival with form extracted from the joint fit model
    lmeFit <- lme(jointFit1$Forms$formYx, data = trainLME, random = list(id = pdDiag(form = jointFit1$Forms$formYz)))
    survFit <- coxph(jointFit1$Forms$formT, data = trainSurv, x = TRUE)
    ## fit the jointfit model
    jointFit1 <- jointModelBayes(lmeFit, survFit, timeVar = "time", n.iter = iter)
    ## calculate auc and prediction error
    auc.vec_jointFit1[i] <- aucJM(jointFit1, newdata = validLME,  Tstart = startTime, Thoriz = PredTime)$auc
    pred.vec_jointFit1[i] <- prederrJM(jointFit1, newdata = validLME,  Tstart = startTime, Thoriz = PredTime)$prederr
    print(i)
  }
  tbl <- data.frame(auc.vec_jointFit1, 
                    pred.vec_jointFit1)
  print(tbl)
  tbl
} 
