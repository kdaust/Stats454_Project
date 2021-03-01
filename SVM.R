library(LiblineaR)

##usually just start here
dat <- readRDS("project556.rds")
outcome <- as.factor(dat$celltype)
load("CleanedFeatures.RData")
##use random forest to cut down variables
set.seed(0)
rf1 <- ranger(x = preds, y = outcome,num.trees = 501,importance = "impurity")
varImp <- importance(rf1)
plot(varImp,pch = 16,col = "#FF0000AA",ylab = "Variable Importance",xlab = "Variable")
abline(h = 1, col = "black")

var2 <- varImp[varImp > 35]
var2 <- names(var2)
par(mfrow = c(1,4))
for(var in var2){
  hist(preds[,var],xlab = "Value", main = var)
}

toUse <- names(varImp[varImp > 1]) ##importance cutoff - kind of arbitrary
preds <- preds[,toUse]

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
    fit.svm <- LiblineaR(data = XTrain, target = YTrain,type = 2, cost = 1)
    pred <- predict(fit.svm,XTest)$predictions
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
