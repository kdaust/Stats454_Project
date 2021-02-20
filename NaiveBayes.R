library(data.table)
library(tidyverse)
library(caret)
library(e1071)
library(naivebayes)
library(GA)
library(Rfast)
library(ranger)
library(rpart)
library(HiClimR)


dat <- readRDS("project556.rds")
preds <- t(dat$count)
temp <- apply(preds,2,FUN = function(x) min(x) == max(x))
preds <- preds[,!temp]

corMat <- fastCor(preds,nSplit = 2, upperTri = T, optBLAS = T, verbose = T)
temp <- which(corMat > 0.8, arr.ind = T)
preds <- preds[,-temp[,1]]
cond <- dat$condt
cond <- fifelse(cond == "Control",0,1)
preds <- cbind(cond, preds)
outcome <- as.factor(dat$celltype)

fit <- rpart(Class ~ .,data = data_x, method = "class")
rf1 <- ranger(x = preds, y = outcome,num.trees = 201,importance = "impurity")

nb <- naiveBayes(x = preds, y = outcome)

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

initial_names <- names(varImp[varImp > 5])
initial <- fifelse(col_names %in% initial_names, 1,0)
ititial <- matrix(data = initial, nrow = 1)

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
             run=5, # max iter without improvement (stopping criteria)
             maxiter = 50, # total runs or generations
             monitor = plot, # plot the result at each iteration
             keepBest = TRUE, # keep the best solution at the end
             parallel = F # allow parallel procesing
)
