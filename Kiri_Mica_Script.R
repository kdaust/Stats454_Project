##Midterm Project R script
##Mica Grant-Hagen, Kiri Daust
library(data.table)
library(tidyverse)
library(caret)
library(naivebayes)
library(GA)
library(ranger)
library(HiClimR)
library(MTPS)
library(kernlab)

###Data processing and cleaning###
##load data
dat <- readRDS("project556.rds")
preds <- t(dat$count)
##remove zero variance columns
temp <- apply(preds,2,FUN = function(x) min(x) == max(x))
preds <- preds[,!temp]

##correlation matrix. This takes a WHILE
corMat <- fastCor(preds,nSplit = 2, upperTri = T, optBLAS = T, verbose = T)
temp <- which(corMat > 0.8, arr.ind = T) ##remove corr > 0.8
preds <- preds[,-temp[,1]]
##add condition
cond <- dat$condt
cond <- fifelse(cond == "Control",0,1)
preds <- cbind(cond, preds)

##usually just start here
outcome <- as.factor(dat$celltype)
#load("CleanedFeatures.RData")
##use random forest to cut down variables
set.seed(0)
rf1 <- ranger(x = preds, y = outcome,num.trees = 501,importance = "impurity")
varImp <- importance(rf1)
plot(varImp,pch = 16,col = "#FF0000AA",ylab = "Variable Importance",xlab = "Variable")
abline(h = 1, col = "black")

##histograms for top 4 variables
var2 <- varImp[varImp > 35]
var2 <- names(var2)
par(mfrow = c(1,4))
for(var in var2){
  hist(preds[,var],xlab = "Value", main = var)
}

toUse <- names(varImp[varImp > 1]) ##importance cutoff - kind of arbitrary
preds <- preds[,toUse]

#############################################
### Naieve Bayes
############################################
##Genetic Algorithm
##fitness function
calc_fitness <- function(vars, data_x, data_y){
  cvFold <- createFolds(y = data_y, k = 3, list = F)
  varNames <- colnames(data_x)
  names_2 <- varNames[vars==1]
  # get the columns of the current solution
  data_sol=data_x[, names_2]
  
  acc <- numeric(length = 3L)
  for(i in 1:3){
    YTrain <- data_y[cvFold != i]
    YTest <- data_y[cvFold == i]
    XTrain <- data_sol[cvFold != i,]
    XTest <- data_sol[cvFold == i,]
    nbfit <- naive_bayes(x = XTrain, y = YTrain,usepoisson = T)
    pred <- predict(nbfit, XTest)
    diff <- pred == YTest
    acc[i] <- length(diff[diff])/length(diff)
  }
  # time for your magic
  #fitness_value=roc_value/q_vars
  
  return(mean(acc))
}

data_x <- preds
data_y <- outcome
# GA parameters
param_nBits=ncol(data_x)
col_names=colnames(data_x)

# initial_names <- names(varImp[varImp > 5])
# initial <- fifelse(col_names %in% initial_names, 1,0)
# ititial <- matrix(data = initial, nrow = 1)

ga_1 <- ga(fitness = function(vars) calc_fitness(vars = vars, 
                                                 data_x =  data_x, 
                                                 data_y = data_y), # custom fitness function
           type = "binary", # optimization data type
           crossover=gabin_uCrossover,  # cross-over method
           elitism = 3, # number of best ind. to pass to next iteration
           pmutation = 0.08, # mutation rate prob
           popSize = 50, # the number of indivduals/solutions
           nBits = param_nBits, # total number of variables
           names=col_names, # variable name
           run=20, # max iter without improvement (stopping criteria)
           maxiter = 50, # total runs or generations
           monitor = plot, # plot the result at each iteration
           keepBest = TRUE, # keep the best solution at the end
           parallel = F # allow parallel procesing
)

vars <- ga_1@solution
vars <- vars[1,]
varNames <- names(vars[vars == 1])

###run cross validation to get accuracy and metrics
pred_bayes <- array(NA, c(length(outcome),10))
allBayes <- data.table()
for(i in 1:10){
  cat(".")
  cvFolds <- createFolds(outcome,k = 5, list = F)
  for(j in 1:5){
    YTrain <- outcome[cvFolds != j]
    YTest <- outcome[cvFolds == j]
    XTrain <- data_x[cvFolds != j,]
    XTest <- data_x[cvFolds == j,]
    nbFit <- naive_bayes(x = XTrain, y = YTrain,usepoisson = F)
    nbPred <- predict(nbFit,XTest)
    pred_bayes[cvFolds == j,i] <- nbPred
    diff <- nbPred == YTest
    temp <- data.table(Rep = i, Fold = j, Acc = length(diff[diff == T])/length(diff))
    allBayes <- rbind(allBayes,temp, fill = T)
  }
}

tp_lda <- array(NA, c(10,5))
fn_lda <- array(NA, c(10,5))
fp_lda <- array(NA, c(10,5))
tn_lda <- array(NA, c(10,5))

for (i in 1:10){
  temp_table <- table(pred_bayes[,i],outcome)
  for (j in 1:5){
    tp_lda[i,j] <- temp_table[j,j]
    fn_lda[i,j] <- sum(temp_table[,j])-tp_lda[i,j]
    fp_lda[i,j] <- sum(temp_table[j,])-tp_lda[i,j]
    tn_lda[i,j] <- sum(temp_table)-tp_lda[i,j]-fn_lda[i,j]-fp_lda[i,j]
  }
}

#Calculating the criteria
miscls <- (fn_lda+fp_lda)/(tp_lda+fp_lda+tn_lda+fn_lda)
precision <- tp_lda/(tp_lda+fp_lda)
recall <- tp_lda/(tp_lda+fn_lda)
f1 <- 2*(precision*recall)/(precision+recall)

celltypenames <- array(c("AT1", "B_Cells", "Basal", "Fibroblasts", "Mast"))

boxplot(f1, main = "F1",  names =celltypenames)
par(mfrow = c(1,2))
boxplot(precision, main = "Precision", names = celltypenames, las = 1)
boxplot(recall, main = "Recall", names = celltypenames, las = 1)

####################################################################
### Elastic Net Multinomial Regression
#####################################################################
###use cv to determine best alpha
cvFolds <- createFolds(outcome,k = 5, list = F)
eNetcv <- data.table()
##this part takes a couple of hours
alphaVals <- seq(0,1,by = 0.1)
for(j in 1:5){
  YTrain <- outcome[cvFolds != j]
  YTest <- outcome[cvFolds == j]
  XTrain <- preds[cvFolds != j,]
  XTest <- preds[cvFolds == j,]
  for(a in 1:11){
    cat(".")
    enet.fit <- cv.glmnet(x = XTrain,y = YTrain,alpha = alphaVals[a],
                          family = "multinomial", type.multinomial = "ungrouped")
    
    enet.pred <- predict(enet.fit, XTest,s = "lambda.1se",type = "class")[,1]
    diff <- enet.pred == YTest
    temp <- data.table(Fold = j, Alpha = alphaVals[a], 
                       Acc = length(diff[diff == T])/length(diff))
    eNetcv <- rbind(eNetcv,temp,fill = T)
  }
}

eNet_sum <- eNetcv[,.(Acc = mean(Acc)), by = .(Alpha)]
eNetcv[,Alpha := as.factor(Alpha)]
ggplot(eNet_sum,aes(x = Alpha, y = Acc))+
  geom_line() +
  labs(y = "Mean Accuracy")+
  theme_classic()

###cross validiation for stats
set.seed(0)
toUse <- names(varImp[varImp > 2]) ##importance cutoff - kind of arbitrary
data_x <- preds[,toUse]
pred_enet <- array(NA, c(length(outcome),10))
allEnet <- data.table()
for(i in 1:10){
  cvFolds <- createFolds(outcome,k = 5, list = F)
  for(j in 1:5){
    cat(".")
    YTrain <- outcome[cvFolds != j]
    YTest <- outcome[cvFolds == j]
    XTrain <- data_x[cvFolds != j,]
    XTest <- data_x[cvFolds == j,]
    nbFit <- cv.glmnet(x = XTrain,y = YTrain,alpha = 0.1,
                       family = "multinomial", type.multinomial = "ungrouped")
    nbPred <- predict(nbFit, XTest,s = "lambda.1se",type = "class")[,1]
    pred_enet[cvFolds == j,i] <- nbPred
    diff <- nbPred == YTest
    temp <- data.table(Rep = i, Fold = j, Acc = length(diff[diff == T])/length(diff))
    
    allEnet <- rbind(allEnet,temp, fill = T)
    
  }
}

tp_lda <- array(NA, c(10,5))
fn_lda <- array(NA, c(10,5))
fp_lda <- array(NA, c(10,5))
tn_lda <- array(NA, c(10,5))

for (i in 1:10){
  temp_table <- table(pred_enet[,i],outcome)
  for (j in 1:5){
    tp_lda[i,j] <- temp_table[j,j]
    fn_lda[i,j] <- sum(temp_table[,j])-tp_lda[i,j]
    fp_lda[i,j] <- sum(temp_table[j,])-tp_lda[i,j]
    tn_lda[i,j] <- sum(temp_table)-tp_lda[i,j]-fn_lda[i,j]-fp_lda[i,j]
  }
}

#Calculating the criteria
miscls <- (fn_lda+fp_lda)/(tp_lda+fp_lda+tn_lda+fn_lda)
precision <- tp_lda/(tp_lda+fp_lda)
recall <- tp_lda/(tp_lda+fn_lda)
f1 <- 2*(precision*recall)/(precision+recall)

celltypenames <- array(c("AT1", "B_Cells", "Basal", "Fibroblasts", "Mast"))

boxplot(f1, main = "F1",  names =celltypenames)
par(mfrow = c(1,2))
boxplot(precision, main = "Precision", names = celltypenames, las = 1)
boxplot(recall, main = "Recall", names = celltypenames, las = 1)

#########################################################3
### LDA
##########################################################
##Genetic Algorithm for lda
calc_fitness <- function(vars, data_x, data_y){
  cvFold <- createFolds(y = data_y, k = 3, list = F)
  varNames <- colnames(data_x)
  names_2 <- varNames[vars==1]
  # get the columns of the current solution
  data_sol=data_x[, names_2]
  
  acc <- numeric(length = 3L)
  for(i in 1:3){
    YTrain <- data_y[cvFold != i]
    YTest <- data_y[cvFold == i]
    XTrain <- data_sol[cvFold != i,]
    XTest <- data_sol[cvFold == i,]
    temp_dat <- data.frame(yy = YTrain, XTrain)
    ldafit <- lda(YTrain~XTrain, type="response")
    pred <- predict(ldafit, newdata =  as.data.frame(XTest))
    diff <- pred == YTest
    acc[i] <- length(diff[diff])/length(diff)
  }
  return(mean(acc))
}

data_x <- preds
data_y <- outcome
# GA parameters
param_nBits=ncol(data_x)
col_names=colnames(data_x)

ga_1 <- ga(fitness = function(vars) calc_fitness(vars = vars, 
                                                 data_x =  data_x, 
                                                 data_y = data_y), # custom fitness function
           type = "binary", # optimization data type
           crossover=gabin_uCrossover,  # cross-over method
           elitism = 3, # number of best ind. to pass to next iteration
           pmutation = 0.08, # mutation rate prob
           popSize = 50, # the number of indivduals/solutions
           nBits = param_nBits, # total number of variables
           names=col_names, # variable name
           run=20, # max iter without improvement (stopping criteria)
           maxiter = 50, # total runs or generations
           monitor = plot, # plot the result at each iteration
           keepBest = TRUE, # keep the best solution at the end
           parallel = F # allow parallel procesing
)

vars <- ga_1@solution
vars <- vars[1,]
varNames <- names(vars[vars == 1])
data_x <- preds[,varNames]

y.pred_LDA <- array(NA, c(length(outcome),10))
allLDA <- data.table()
for(i in 1:10){
  cat(".")
  cvFolds <- createFolds(outcome,k = 5, list = F)
  for(j in 1:5){
    YTrain <- outcome[cvFolds != j]
    YTest <- outcome[cvFolds == j]
    XTrain <- data_x[cvFolds != j,]
    XTest <- data_x[cvFolds == j,]
    temp_dat <- data.frame(yy = YTrain, XTrain)
    lda.fit <- lda(x = XTrain,grouping = YTrain)
    pred <- predict(lda.fit, XTest)$class
    y.pred_LDA[cvFolds == j,i] <- pred
    diff <- pred == YTest
    temp <- data.table(Rep = i, Fold = j, Acc = length(diff[diff == T])/length(diff))
    
    allLDA <- rbind(allLDA,temp, fill = T)
  }
}

tp_lda <- array(NA, c(10,5))
fn_lda <- array(NA, c(10,5))
fp_lda <- array(NA, c(10,5))
tn_lda <- array(NA, c(10,5))

#table(y.pred_LDA[,1],outcome)

for (i in 1:10){
  temp_table <- table(y.pred_LDA[,i],outcome)
  for (j in 1:5){
    
    tp_lda[i,j] <- temp_table[j,j]
    fn_lda[i,j] <- sum(temp_table[,j])-tp_lda[i,j]
    fp_lda[i,j] <- sum(temp_table[j,])-tp_lda[i,j]
    tn_lda[i,j] <- sum(temp_table)-tp_lda[i,j]-fn_lda[i,j]-fp_lda[i,j]
    
  }
}

#Calculating the criteria
miscls <- (fn_lda+fp_lda)/(tp_lda+fp_lda+tn_lda+fn_lda)
precision <- tp_lda/(tp_lda+fp_lda)
recall <- tp_lda/(tp_lda+fn_lda)
f1 <- 2*(precision*recall)/(precision+recall)

celltypenames <- array(c("AT1", "B_Cells", "Basal", "Fibroblasts", "Mast_Cells"))

boxplot(f1, main = "Fi values for classes for LDA",  names =celltypenames)
boxplot(miscls, main = "Misclassification values for classes for LDA", names  = celltypenames)
boxplot(precision, main = "Precision values for classes for LDA", names = celltypenames)
boxplot(recall, main = "Recall values for classes for LDA", names = celltypenames)

##########################################################
### KNN
##########################################################
              
##Genetic Algorithm
##fitness function
calc_fitness <- function(vars, data_x, data_y){
  cvFold <- createFolds(y = data_y, k = 3, list = F)
  varNames <- colnames(data_x)
  names_2 <- varNames[vars==1]
  # get the columns of the current solution
  data_sol=data_x[, names_2]
  
  acc <- numeric(length = 3L)
  for(i in 1:3){
    cat(".")
    YTrain <- data_y[cvFold != i]
    YTest <- data_y[cvFold == i]
    XTrain <- data_sol[cvFold != i,]
    XTest <- data_sol[cvFold == i,]
    temp_dat <- data.frame(yy = YTrain, XTrain)
    knn.fit1 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=3)
    knn.fit2 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=5)
    knn.fit3 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=7)
    
    diff1 <- knn.fit1 == YTest
    diff2 <- knn.fit2 == YTest
    diff3 <- knn.fit3 == YTest
    diff <- max(c(length(diff1), length(diff2), length(diff3)))
    
    acc[i] <- length(diff[diff])/length(diff)
  }

  
  return(mean(acc))
}

data_x <- preds
data_y <- outcome
# GA parameters
param_nBits=ncol(data_x)
col_names=colnames(data_x)


ga_1_knn <- ga(fitness = function(vars) calc_fitness(vars = vars, 
                                                     data_x =  data_x, 
                                                     data_y = data_y), # custom fitness function
               type = "binary", # optimization data type
               crossover=gabin_uCrossover,  # cross-over method
               elitism = 3, # number of best ind. to pass to next iteration
               pmutation = 0.08, # mutation rate prob
               popSize = 50, # the number of indivduals/solutions
               nBits = param_nBits, # total number of variables
               names=col_names, # variable name
               run=20, # max iter without improvement (stopping criteria)
               maxiter = 50, # total runs or generations
               monitor = plot, # plot the result at each iteration
               keepBest = TRUE, # keep the best solution at the end
               parallel = F # allow parallel procesing
)

vars_knn <- ga_1_knn@solution
vars_knn <- vars_knn[1,]
varNames_knn <- names(vars_knn[vars_knn == 1])

data_x <- preds[,varNames_knn]


##cross validation of accuracy

y.pred_Knn <- array(NA, c(length(outcome),10,3))
allKNN <- data.table()
for(i in 1:10){
  cat(".")
  cvFolds <- createFolds(outcome,k = 5, list = F)
  for(j in 1:5){
    YTrain <- outcome[cvFolds != j]
    YTest <- outcome[cvFolds == j]
    XTrain <- data_x[cvFolds != j,]
    XTest <- data_x[cvFolds == j,]
    knn.fit1 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=3)
    knn.fit2 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=5)
    knn.fit3 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=7)
    

    y.pred_Knn[cvFolds == j,i,1] <- knn.fit1
    y.pred_Knn[cvFolds == j,i,2] <- knn.fit2
    y.pred_Knn[cvFolds == j,i,3] <- knn.fit3
    #diff <- pred == YTest
    diff1 <- knn.fit1 == YTest
    diff2 <- knn.fit2 == YTest
    diff3 <- knn.fit3 == YTest
    temp <- data.table(Rep = i, Fold = j, Acc_3 = length(diff1[diff1 == T])/length(diff1), Acc_5 = length(diff2[diff2 == T])/length(diff2), Acc_7 = length(diff3[diff3 == T])/length(diff3))
    
    allKNN <- rbind(allKNN,temp, fill = T)
  }
}
knn_mean_acc <- array(NA, c(10,3))
knn_mean_acc[,1] <- allKNN[,.(MSE = mean(Acc_3)), by = .(Rep)]$MSE
knn_mean_acc[,2] <- allKNN[,.(MSE = mean(Acc_5)), by = .(Rep)]$MSE
knn_mean_acc[,3] <- allKNN[,.(MSE = mean(Acc_7)), by = .(Rep)]$MSE
best_k <- sort(colMeans(knn_mean_acc),index.return = TRUE, decreasing = TRUE)$ix[1]
tp_knn <- array(NA, c(10,5))
fn_knn <- array(NA, c(10,5))
fp_knn <- array(NA, c(10,5))
tn_knn <- array(NA, c(10,5))


for (i in 1:10){
  print(i)
  temp_table <- table(y.pred_Knn[,i,best_k],outcome)
  for (j in 1:5){
    
    tp_knn[i,j] <- temp_table[j,j]
    fn_knn[i,j] <- sum(temp_table[,j])-tp_knn[i,j]
    fp_knn[i,j] <- sum(temp_table[j,])-tp_knn[i,j]
    tn_knn[i,j] <- sum(temp_table)-tp_knn[i,j]-fn_knn[i,j]-fp_knn[i,j]
    
  }
}


#Calculating the criteria
miscls <- (fn_knn+fp_knn)/(tp_knn+fp_knn+tn_knn+fn_knn)
precision <- tp_knn/(tp_knn+fp_knn)
recall <- tp_knn/(tp_knn+fn_knn)
f1 <- 2*(precision*recall)/(precision+recall)
celltypenames <- array(c("AT1", "B_Cells", "Basal", "Fibroblasts", "Mast_Cells"))

boxplot(f1, main = "Fi values for classes for KNN",  names =celltypenames)
boxplot(miscls, main = "Misclassification values for classes for KNN", names  = celltypenames)
boxplot(precision, main = "Precision values for classes for KNN", names = celltypenames)
boxplot(recall, main = "Recall values for classes for KNN", names = celltypenames)
