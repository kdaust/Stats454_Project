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
library(MTPS)

##load data
dat <- readRDS("project556.rds")
preds <- t(dat$count)
##remove zero variance columns
temp <- apply(preds,2,FUN = function(x) min(x) == max(x))
preds <- preds[,!temp]

##correlation matrix. This takes a WHILE
corMat <- fastCor(preds,nSplit = 2, upperTri = T, optBLAS = T, verbose = T)
temp <- which(corMat > 0.7, arr.ind = T) ##remove corr > 0.7
preds <- preds[,-temp[,1]]
##add condition
cond <- dat$condt
cond <- fifelse(cond == "Control",0,1)
preds <- cbind(cond, preds)

##usually just start here
outcome <- as.factor(dat$celltype)
load("CleanedFeatures.RData")
##use random forest to cut down variables
rf1 <- ranger(x = preds, y = outcome,num.trees = 501,importance = "impurity")
varImp <- importance(rf1)
plot(varImp)
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

##varnames from genetic alg
varNames <- c("cond", "PLA2G2A", "RPL11", "SFN", "LAPTM5", "RPS8", "TSPAN1", 
              "CYP4B1", "LEPR", "RPL5", "CD53", "S100A6", "RPS27", "LY9", "RCSD1", 
              "PRRX1", "GLUL", "RGS1", "ELF3", "FCMR", "LAMB3", "C1orf198", 
              "EPCAM", "SPTBN1", "CNRIP1", "CAPG", "IL1RL1", "FHL2", "CXCR4", 
              "ZEB2", "IGFBP2", "SP140", "COL6A3", "NFKBIZ", "HCLS1", "CPA3", 
              "ARHGEF26", "RPL35A", "SLC34A2", "LIMCH1", "HOPX", "ANXA3", "ARHGAP24", 
              "SPARCL1", "HPGDS", "ADH1B", "RPS3A", "SFRP2", "HPGD", "WWC2", 
              "ITGA2", "MEF2C", "HLA-A", "HLA-E", "DDR1", "HLA-B", "LTB", "AGER", 
              "HLA-DMA", "RPS18", "DLK2", "VEGFA", "LAMA4", "AKAP12", "TSPAN13", 
              "SFRP4", "STK17A", "CLDN4", "UPK3B", "CAV2", "MET", "TMSB4X", 
              "TIMP1", "CHRDL1", "FHL1", "ARHGEF6", "TUSC3", "ASAH1", "CLU", 
              "SCARA5", "RAB11FIP1", "COL14A1", "ALDH1A1", "ANXA1", "GSN", 
              "RPL12", "PTGDS", "CD151", "CTSD", "RPS13", "CD44", "MDK", "SERPING1", 
              "MS4A1", "MS4A15", "FTH1", "FOLR1", "UCP2", "MMP7", "POU2AF1", 
              "CADM1", "MPZL2", "RPS25", "VIM", "ALOX5", "RTKN2", "SRGN", "SPOCK2", 
              "RPS24", "PAPSS2", "C1R", "SLC2A3", "ARHGDIB", "MGST1", "KRT5", 
              "IGFBP6", "CD63", "CPM", "LUM", "TSC22D1", "ITM2B", "KLF5", "LMO7", 
              "SCEL", "COL4A2", "IFI27", "FBN1", "HDC", "LINC00926", "ANXA2", 
              "RPLP1", "RPS2", "CORO1A", "PRSS8", "CRNDE", "MMP2", "LPCAT2", 
              "MT1X", "PLLP", "COTL1", "SERPINF1", "PMP22", "TNFRSF13B", "RPL19", 
              "IGFBP4", "KRT19", "KRT17", "SCPEP1", "CD79B", "ITGB4", "CHST9", 
              "RAB27B", "DSTN", "ID1", "MYL9", "EYA2", "PLPP2", "FSTL3", "RPS15", 
              "ICAM1", "KLF2", "FXYD3", "FXYD5", "HSPB6", "GMFG", "RPS16", 
              "LTBP4", "CD79A", "POU2F2", "CXCL17", "CD37", "RPS9", "RPL28", 
              "VPREB3", "TIMP3", "RAC2", "RPL3", "TNFRSF13C", "SAMSN1", "BACE2", 
              "COL6A1", "MT-CO1")

data_x <- preds[,varNames]

##cross validation of accuracy
allBayes <- data.table()
for(i in 1:10){
  cat(".")
  cvFolds <- createFolds(outcome,k = 5, list = F)
  for(j in 1:5){
    YTrain <- outcome[cvFolds != j]
    YTest <- outcome[cvFolds == j]
    XTrain <- data_x[cvFolds != j,]
    XTest <- data_x[cvFolds == j,]
    nbFit <- naive_bayes(x = XTrain, y = YTrain,usepoisson = T)
    nbPred <- predict(nbFit,XTest)
    diff <- nbPred == YTest
    temp <- data.table(Rep = i, Fold = j, Acc = length(diff[diff == T])/length(diff))
    allBayes <- rbind(allBayes,temp, fill = T)
  }
}

bayesAvg <- allBayes[,.(MSE = mean(Acc)), by = .(Rep)]

lda.fit <- lda(x = XTrain,grouping = YTrain)
pred <- predict(lda.fit, XTest)$class

cvFolds <- createFolds(outcome,k = 5, list = F)
j = 1
YTrain <- outcome[cvFolds != j]
YTest <- outcome[cvFolds == j]
XTrain <- preds[cvFolds != j,]
XTest <- preds[cvFolds == j,]

tic()
enet.fit <- cv.glmnet(x = XTrain,y = YTrain,alpha = 0.8,
                   family = "multinomial", type.multinomial = "ungrouped")
toc()
enet.pred <- predict(enet.fit, XTest,type = "class")[,1]
diff <- enet.pred == YTest
table(diff)

library(LiblineaR)
fit.svm <- LiblineaR(data = XTrain, target = YTrain,type = 2, cost = 1)
pred <- predict(fit.svm,XTest)$predictions
diff <- pred == YTest
table(diff)
