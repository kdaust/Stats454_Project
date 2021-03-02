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
set.seed(0)
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
print(varNames_knn)

varNames_LDA <- c(
  "cond",       "PRKCZ",      "PLA2G2A",    "SH3BGRL3",   "CD52",       "SMAP2" ,     "PLPP3",      "TACSTD2",   
  "CYR61",      "OLFML3",     "CTSK",       "S100A10",    "S100A4",     "S100A2",     "S100A16",    "S100A13"   ,
  "EFNA1" ,     "RAB25",      "NTRK1",      "FCER1A",     "ITLN2",      "FCER1G",     "FCGR2B",     "FCRLA"     ,
  "ATP1B1" ,    "SELL",       "PRRX1",      "RGS1",       "ELF3",       "PRELP",      "FCMR",       "CD34"      ,
  "LAMB3"   ,   "GCSAML",     "HPCAL1",     "SPTBN1",     "TMSB10",     "CAPG",       "VAMP8",      "IL1RL1"    ,
  "EPB41L5",    "GYPC",       "NCKAP5",     "ARHGAP15",   "NR4A2",      "CYTIP",      "TFPI",       "COL3A1"    ,
  "STK17B",     "CFLAR",      "EEF1B2",     "FN1",        "IGFBP2",     "TNS1",       "COL4A4",     "SP140"     ,
  "PTMA",       "COL6A3",     "BHLHE40",    "FBLN2",      "RPL15",      "RPSA" ,      "SEMA3B",     "CACNA2D2"  ,
  "ALCAM",      "CCDC80",     "GATA2",      "CLDN18",     "CPA3",       "GPR87" ,     "VEPH1" ,     "RARRES1"   ,
  "TNFSF10",    "CLDN1" ,     "ATP13A4",    "APOD",       "SPON2",      "WFS1",       "LIMCH1",     "IGFBP7"    ,
  "AREG",       "ANXA3",      "SPARCL1",    "ADH7",       "NPNT",       "SFRP2",      "MARCH1",     "CPE"       ,
  "GPM6A",      "BTF3",       "MEF2C",      "TNFAIP8",    "CYSTM1",     "PDGFRB",     "CD74",       "SLIT3"     ,
  "RNF130",     "CD83",       "GMPR",       "HLA-C",      "GPSM3",      "MRPL14",     "CLIC5",      "ADGRF5"    ,
  "DST",        "COL12A1",    "HMGN3",      "CD24",       "LAMA4",      "LAMA2",      "TCF21",      "TNFAIP3"   ,
  "PERP",       "AKAP12",     "SYTL3",      "EZR",        "TAGAP",      "AHR" ,       "GPNMB",      "SFRP4"     ,
  "STK17A",     "EGFR",       "ELN",        "HSPB1",      "GNG11",      "BRI3" ,      "CAV2",       "MET"       ,
  "ANOS1",      "SRPX",       "MAOB" ,      "TIMP1",      "ITM2A" ,     "TSC22D3",    "ACSL4",      "CHRDL1"    ,
  "PLS3",       "FHL1",       "ARHGEF6",    "HMGB3",      "DLC1",       "PEBP4" ,     "ADAM28",     "DPYSL2"    ,
  "CLU" ,       "SARAF",      "SDCBP",      "PABPC1",     "COL14A1" ,   "DENND3",     "SLC1A1",     "IL33"      ,
  "PLIN2",      "AQP3",       "CCDC107",    "ALDH1A1",    "ANXA1" ,     "IFITM3",     "CD151",      "TSPAN4"    ,
  "CTSD"  ,     "SYT8",       "OLFML1",     "CAT",        "EHF" ,       "CD44",       "MDK",        "SPI1"      ,
  "MYRF",       "FTH1" ,      "EEF1G",      "LRRN4CL",    "PRDX5",      "NEAT1",      "MALAT1",     "CST6"      ,
  "CATSPER1",   "GSTP1" ,     "ALDH3B1",    "CCND1",      "POU2AF1",    "NNMT",       "CADM1",      "TMPRSS4"   ,
  "CXCR5",      "THY1",       "SFTA1P",     "VIM",        "APBB1IP",    "ALOX5" ,     "ADIRF",      "PAPSS2"    ,
  "PDLIM1",     "CD9" ,       "CD27",       "A2M",        "CD69" ,      "MGP",        "MGST1",      "LDHB"      ,
  "KRT7",       "KRT8" ,      "IGFBP6",     "MYL6",       "NACA",       "BTG1",       "RPLP0",      "RPL21"     ,
  "MEDAG",      "RGCC",       "TSC22D1",    "TPT1",       "LCP1",       "ITM2B",      "SPRYD7",     "KLF5"      ,
  "SCEL",       "ABCC4",      "COL4A2",     "NKX2-1",     "NPC2" ,      "FBLN5",      "B2M",        "FGF7"      ,
  "HDC" ,       "LINC00926",  "ANXA2",      "DAPK2",      "BCL2A1",     "MFGE8",      "MSLN",       "TPSB2"     ,
  "NSMCE1",     "CD19",       "CORO1A",     "PRSS8",      "PLLP",       "CDH3",       "CRISPLD2",   "TNFRSF13B" ,
  "ALDH3A1",    "IGFBP4",     "TNS4",       "KRT15",      "ABCA6",      "ST6GALNAC1", "LAMA3",      "CHST9"     ,
  "DSC3",       "TCF4",       "SERPINB5",   "FLRT3",      "CST3",       "SLPI",       "PTGIS",      "FSTL3"     ,
  "OAZ1",       "ICAM1",      "KLF2" ,      "ISYNA1",     "FXYD3",      "FXYD5",      "DMKN",       "HSPB6"     ,
  "SPINT2",     "LGALS7",     "CEACAM6",    "BCAM",       "FTL",        "CD37",       "RPL13A",     "NAPSA"     ,
  "CLEC11A",    "SIGLEC6",    "LGALS1",     "GRAP2",      "SAMSN1",     "CXADR",      "APP",        "ETS2"      ,
  "BACE2")
##varnames from genetic alg
varNames_NB <- c("cond", "PLA2G2A", "RPL11", "SFN", "LAPTM5", "RPS8", "TSPAN1",
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

varNames_knn <- c("cond",       "MMP23B",     "PLA2G2A",    "LAPTM5",     "FHL3",       "TMEM125",    "TSPAN1",     "PODN",      
                  "PLPP3",      "ST6GALNAC5", "ATP1A1",     "CTSK",       "S100A10",    "S100A4",     "S100A2",     "S100A16",   
                  "EFNA1"  ,    "RAB25",      "FCER1A",     "CD48" ,      "FCER1G",     "FCGR2B",     "SELL",       "GLUL",      
                  "PTGS2"   ,   "RGS2",       "PTPRC",      "ELF3" ,      "PRELP",      "CD55",       "CD34",       "EPHX1",     
                  "C1orf198" ,  "GCSAML",     "LAPTM4A",    "EPCAM" ,     "SPTBN1",     "CNRIP1",     "CAPG" ,      "SFTPB" ,    
                  "AFF3",       "FHL2",       "NCKAP5",     "ARHGAP15",   "CYTIP" ,     "TFPI",       "STK17B",     "ANKRD44",   
                  "CFLAR",      "EEF1B2",     "COL4A4",     "SP140",      "SP100",      "PTMA",       "COL6A3",     "BHLHE40" ,  
                  "FBLN2" ,     "RPL15",      "RHOA" ,      "COL8A1",     "NFKBIZ",     "IFT57",      "CCDC80",     "PARP15",    
                  "MGLL"  ,     "GATA2",      "CPA3",       "VEPH1",      "RARRES1",    "TNIK",       "CLDN1",      "ATP13A4",   
                  "APOD",       "HTRA3",      "SOD3",       "SLC34A2",    "LIMCH1",     "KIT",        "ANXA3",      "ARHGAP24",  
                  "SPARCL1",    "ADH7",       "BANK1",      "NPNT",       "CPE" ,       "PDLIM3",     "ITGA2",      "PLPP1",     
                  "TNFAIP8" ,   "HINT1",      "CD74",       "ANXA6",      "SLIT3",      "LCP2",       "NPM1" ,      "LY86",      
                  "DSP",        "GMPR",       "DDR1" ,      "LTB",        "GPSM3",      "HLA-DMB",    "VEGFA",      "MRPL14",    
                  "CLIC5",      "ADGRF5",     "COL12A1",    "HMGN3",      "CD24" ,      "LAMA4",      "LAMA2",      "ARHGAP18",  
                  "TCF21",      "TNFAIP3",    "PERP",       "STX11",      "TAGAP",      "AHR",        "GPNMB",      "ELN",       
                  "CAV2",       "ANOS1",      "MAOB",       "IL2RG",      "ITM2A",      "ACSL4",      "CHRDL1",     "PLS3",      
                  "FHL1" ,      "CLIC2",      "BLK",        "DLC1",       "ASAH1",      "ADAM28" ,    "CLU" ,       "RAB11FIP1", 
                  "IL33",       "PLIN2",      "FAM189A2",   "ALDH1A1",    "ANXA1",      "PTGDS",      "NPDC1",      "TSPAN4" ,   
                  "CTSD",       "EHF",        "CD44",       "SPI1",       "FTH1" ,      "AHNAK",      "EEF1G" ,     "LRRN4CL",   
                  "NEAT1",      "MALAT1",     "ALDH3B1",    "CCND1",      "FOLR1",      "POU2AF1",   "NNMT",       "VSIG2",     
                  "NTM" ,       "KLF6",       "VIM" ,       "NEBL",       "ALOX5",      "ARID5B",     "RTKN2",      "DDIT4",     
                  "PLAC9",      "ADIRF",      "PAPSS2",     "CD9" ,       "A2M" ,       "CD69" ,      "LDHB" ,      "KRT7" ,     
                  "KRT5",       "CD63",       "NACA",       "LRP1",       "CPM" ,       "LUM",        "BTG1",       "IGF1",      
                  "RPL21",      "MEDAG",      "TSC22D1",    "SCEL" ,      "NDFIP2",     "GAS6",       "NKX2-1",     "CRIP2",     
                  "IGHM",       "HDC",        "LINC00926",  "ANXA2",      "ISLR" ,      "BCL2A1",     "IL16",       "AKAP13",    
                  "ISG20",      "MFGE8",      "TPSB2",      "TPSD1",      "SOCS1",      "NUPR1" ,     "CORO1A",     "PRSS8",     
                  "LPCAT2",     "CDH3",       "CRISPLD2",   "IRF8",       "CYBA" ,      "PMP22",      "TNFRSF13B",  "IGFBP4",    
                  "CD79B" ,     "ABCA6",      "ITGB4",      "UNC13D",     "ST6GALNAC1", "MYL12A",     "MYL12B",     "DSC3",      
                  "RAB27B",     "NEDD4L",     "SERPINB5",   "CYB5A",      "ID1" ,       "STK4",       "SLPI" ,      "WFDC2",     
                  "PLTP",       "ADGRE2",     "ISYNA1",     "FXYD1",      "COX7A1",     "SPINT2",     "RPS19",      "BCAM",      
                  "EMP3" ,      "FTL",        "NAPSA",      "CLEC11A",    "SIGLEC6" ,   "GGT1",       "TIMP3",      "LGALS1",    
                  "GRAP2" ,     "TNFRSF13C",  "ETS2" )

data_x <- preds[,varNames_knn]
#data_x <- data_x[,varNames_knn]


##cross validation of accuracy
# allBayes <- data.table()
# for(i in 1:10){
#   cat(".")
#   cvFolds <- createFolds(outcome,k = 5, list = F)
#   for(j in 1:5){
#     YTrain <- outcome[cvFolds != j]
#     YTest <- outcome[cvFolds == j]
#     XTrain <- data_x[cvFolds != j,]
#     XTest <- data_x[cvFolds == j,]
#     nbFit <- naive_bayes(x = XTrain, y = YTrain,usepoisson = T)
#     nbPred <- predict(nbFit,XTest)
#     diff <- nbPred == YTest
#     temp <- data.table(Rep = i, Fold = j, Acc = length(diff[diff == T])/length(diff))
#     allBayes <- rbind(allBayes,temp, fill = T)
#   }
# }

#bayesAvg <- allBayes[,.(MSE = mean(Acc)), by = .(Rep)]

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
    #temp_dat <- data.frame(yy = YTrain, XTrain)
    knn.fit1 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=3)
    knn.fit2 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=5)
    knn.fit3 <- knn(as.data.frame(XTrain), as.data.frame(XTest), YTrain, k=7)
    
    # diff1 <- knn.fit1 == YTest
    # diff2 <- knn.fit2 == YTest
    # diff3 <- knn.fit3 == YTest
    # diff <- max(c(length(diff1), length(diff2), length(diff3)))
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

#table(y.pred_LDA[,1],outcome)

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

