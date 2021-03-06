##Stats 454 Final project
##Mica Grant-Hagen and Kiri Daust
##Note that the overall structure of this script is 
##based on Dr. Zhang's elastic net script

library(MTPS)
require("xgboost")
require(superml)
require(caret)
require(data.table)
library(ggplot2)
#read data
data.list=readRDS("final556.rds")

##just some functions for doing tests
recall <- function(true, predicted){
  predicted <- as.numeric(as.character(predicted))
  tp <- predicted[true == 1]
  tp <- sum(tp)
  condPos <- sum(true)
  return(tp/condPos)
}

precision <- function(true, predicted){
  predicted <- as.numeric(as.character(predicted))
  tp <- predicted[true == 1]
  tp <- sum(tp)
  predPos <- sum(predicted)
  return(tp/predPos)
}

F1 <- function(true,predicted){
  predicted <- as.numeric(as.character(predicted))
  RC <- recall(true,predicted)
  PR <- precision(true,predicted)
  res <- (RC*PR)/(RC+PR)
  return(res*2)
}


##########set up the dataset##########
#change fabriclung to binary
FabricLung=data.frame(FabricLung=ifelse(data.list$condt=="ILD",1,0))
#change cell type to binary 
#binomial
cell_type=data.frame(B_Cells=ifelse(data.list$celltype=="B_Cells",1,0),
                     Mesothelial_Cells=ifelse(data.list$celltype=="Mesothelial_Cells",1,0),
                     Myofibroblasts=ifelse(data.list$celltype=="Myofibroblasts",1,0),
                     pDCs=ifelse(data.list$celltype=="pDCs",1,0),
                     Smooth_Muscle_Cells=ifelse(data.list$celltype=="Smooth_Muscle_Cells",1,0))
#cell type label 
cell_factor=data.frame(celltype=factor(data.list$celltype))

#form matrix for each cell type
dat.B_Cells=data.frame(yy=cell_type[,1],FabricLung,t(data.list$count))
dat.Mesothelial_Cells=data.frame(yy=cell_type[,2],FabricLung,t(data.list$count))
dat.Myofibroblasts=data.frame(yy=cell_type[,3],FabricLung,t(data.list$count))
dat.pDCs=data.frame(yy=cell_type[,4],FabricLung,t(data.list$count))
dat.Smooth_Muscle_Cells=data.frame(yy=cell_type[,5],FabricLung,t(data.list$count))

##########generate 10 different version of 5-fold splits##########
set.seed(0)
nrep=10 # replication 10 times
nfolds=5 # 5-fold Cross valiation 

###"idx_cv.csv" can be read by read.csv()
idx.cv=read.csv(file="idx_cv.csv",row.names = 1)

###hyperparameter tuning: we only present one celltype example of parameter tuning here
##all others were conducted in a similar way

##tune max number of iterations
response <- dat.Smooth_Muscle_Cells$yy##specify celltype
xmat <- dat.Smooth_Muscle_Cells[,-1]
xmat <- as.matrix(xmat)
datCurr <- dat.Smooth_Muscle_Cells
init.parms <- list(learning_rate =0.05,
                   n_estimators=1000,
                   max_depth=8,
                   min_child_weight=1,
                   gamma=0.1,
                   subsample=0.8,
                   colsample_bytree=0.8,
                   objective= 'binary:logistic',
                   scale_pos_weight=1)
##run cv function to find optimal number  of iterations
cvres <- xgb.cv(init.parms,data = xmat, label = response,
                nrounds = 500, nfold = 5, metrics = "logloss",
                early_stopping_rounds = 50, print_every_n = 10)

##now setup to do grid search over max_depth and min child weight
plist <- list(max_depth = seq(2,25,by = 5),
              min_child_weight = seq(1,6,2))

##for simplicity, we used the superml package structure here
##which has nice functions for cross-validation grid search
xg <- XGBTrainer$new(learning_rate =0.05,
                     n_estimators=328, ##this is based on the results from above
                     gamma=0.1,
                     subsample=0.8,
                     colsample_bytree=0.8,
                     objective= 'binary:logistic')
gst <- GridSearchCV$new(trainer = xg,
                        parameters = plist,
                        n_folds = 3,
                        scoring = "f1")
gst$fit(datCurr,"yy")##run grid search
gst$best_iteration()##print results

##now try max depth with a finer range of values around the previous best
plist2 <- list(max_depth = c(3,4,5))
gst2 <- GridSearchCV$new(trainer = xg,
                         parameters = plist2,
                         n_folds = 3,
                         scoring = "f1")
gst2$fit(datCurr,"yy")
gst2$best_iteration()

##set parameters based on grid search results
xg$max_depth = 3
xg$min_child_weight = 1

##finally, tune gamma
gst3 <- GridSearchCV$new(trainer = xg,
                                 parameters = list(
                                   gamma = c(0,0.1,0.2,0.3,0.4,0.5)),
                                 n_folds = 3,
                                 scoring = "f1")
gst3$fit(datCurr,"yy")
gst3$best_iteration()


##########create matrix that will save the predict result##########
#predict probability for each cell type
y.pred=matrix(NA,nrow=nrow(idx.cv),ncol=ncol(cell_type))
colnames(y.pred)=c("B_Cells",
                   "Mesothelial_Cells","Myofibroblasts",
                   "pDCs","Smooth_Muscle_Cells")

#predict label based on the predict probability 
cell.pred=matrix(NA,nrow = nrow(idx.cv))

#F1 score 
F1=matrix(NA,nrow = nrep,ncol =ncol(cell_type))
colnames(F1)=c("B_Cells",
               "Mesothelial_Cells","Myofibroblasts",
               "pDCs","Smooth_Muscle_Cells")

##########fit models with tune hyperparameters##########
for(ii in 1:nrep){
  for(jj in 1:nfolds){
    x.train=as.matrix(dat.B_Cells[idx.cv[,ii]!=jj,-1])
    x.test=as.matrix(dat.B_Cells[idx.cv[,ii]==jj,-1])
    #B_Cells
    y.train.B_Cells=as.matrix(as.matrix(dat.B_Cells[idx.cv[,ii]!=jj,1]))
    #y.test.B_Cells=as.matrix(dat.B_Cells[idx.cv[,ii]==jj,1])
    model.ela.B_Cells=xgboost(data = x.train, label = y.train.B_Cells, 
                              max.depth = 4, gamma = 0, min_child_weight = 1, nrounds = 67,
                              objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,1]=predict(model.ela.B_Cells,x.test,type = "response")
    #Mesothelial_Cells
    y.train.Mesothelial_Cells=as.matrix(as.matrix(dat.Mesothelial_Cells[idx.cv[,ii]!=jj,1]))
    #y.test.Mesothelial_Cells=as.matrix(dat.Mesothelial_Cells[idx.cv[,ii]==jj,1])
    model.ela.Mesothelial_Cells=xgboost(data = x.train, label = y.train.Mesothelial_Cells, 
                                        max.depth = 4, gamma = 0, min_child_weight = 1, nrounds = 75,
                                        objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,2]=predict(model.ela.Mesothelial_Cells,x.test,type = "response")
    #Myofibroblasts
    y.train.Myofibroblasts=as.matrix(as.matrix(dat.Myofibroblasts[idx.cv[,ii]!=jj,1]))
    #y.test.Myofibroblasts=as.matrix(dat.Myofibroblasts[idx.cv[,ii]==jj,1])
    model.ela.Myofibroblasts=xgboost(data = x.train, label = y.train.Myofibroblasts, 
                                     max.depth = 4, gamma = 0, min_child_weight = 1, nrounds = 136,
                                     objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,3]=predict(model.ela.Myofibroblasts,x.test,type = "response")
    #pDCs
    y.train.pDCs=as.matrix(as.matrix(dat.pDCs[idx.cv[,ii]!=jj,1]))
    #y.test.pDCs=as.matrix(dat.pDCs[idx.cv[,ii]==jj,1])
    model.ela.pDCs=xgboost(data = x.train, label = y.train.pDCs, 
                           max.depth = 4, gamma = 0, min_child_weight = 1, nrounds = 27,
                           objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,4]=predict(model.ela.pDCs,x.test,type = "response")
    #Smooth_Muscle_Cells
    y.train.Smooth_Muscle_Cells=as.matrix(as.matrix(dat.Smooth_Muscle_Cells[idx.cv[,ii]!=jj,1]))
    #y.test.Smooth_Muscle_Cells=as.matrix(dat.Smooth_Muscle_Cells[idx.cv[,ii]==jj,1])
    model.ela.Smooth_Muscle_Cells=xgboost(data = x.train, label = y.train.Smooth_Muscle_Cells, learning_rate = 0.08,
                                          max.depth = 3, gamma = 0, min_child_weight = 2, nrounds = 250,
                                          objective = "binary:logistic",subsample=0.9,
                                          colsample_bytree=0.9, eta = 0.8,booster = "gbtree")
    y.pred[idx.cv[,ii]==jj,5]=predict(model.ela.Smooth_Muscle_Cells,x.test,type = "response")
  }
  
  #the predict label is based on the highest probability of each cell type
  cell.pred[apply(y.pred,1,which.max)==1]="B_Cells"
  cell.pred[apply(y.pred,1,which.max)==2]="Mesothelial_Cells"
  cell.pred[apply(y.pred,1,which.max)==3]="Myofibroblasts"
  cell.pred[apply(y.pred,1,which.max)==4]="pDCs"
  cell.pred[apply(y.pred,1,which.max)==5]="Smooth_Muscle_Cells"
  
  #criteria for each 
  #for B_Cells
  tp=length(which(cell_factor=="B_Cells"&cell.pred=="B_Cells"))#true positive
  fp=length(which(cell_factor!="B_Cells"&cell.pred=="B_Cells"))#false positive
  fn=length(which(cell_factor=="B_Cells"&cell.pred!="B_Cells"))#false negative
  precision.binary=tp/(tp+fp)#precision.binary
  recall.binary=tp/(tp+fn)#recall.binary
  F1[ii,1]=2*(precision.binary*recall.binary)/(precision.binary+recall.binary)#F1.B_Cells score
  #for Mesothelial_Cells
  tp=length(which(cell_factor=="Mesothelial_Cells"&cell.pred=="Mesothelial_Cells"))#true positive
  fp=length(which(cell_factor!="Mesothelial_Cells"&cell.pred=="Mesothelial_Cells"))#false positive
  fn=length(which(cell_factor=="Mesothelial_Cells"&cell.pred!="Mesothelial_Cells"))#false negative
  precision.binary=tp/(tp+fp)#precision.binary
  recall.binary=tp/(tp+fn)#recall.binary
  F1[ii,2]=2*(precision.binary*recall.binary)/(precision.binary+recall.binary)#F1.Mesothelial_Cells score
  #for Myofibroblasts
  tp=length(which(cell_factor=="Myofibroblasts"&cell.pred=="Myofibroblasts"))#true positive
  fp=length(which(cell_factor!="Myofibroblasts"&cell.pred=="Myofibroblasts"))#false positive
  fn=length(which(cell_factor=="Myofibroblasts"&cell.pred!="Myofibroblasts"))#false negative
  precision.binary=tp/(tp+fp)#precision.binary
  recall.binary=tp/(tp+fn)#recall.binary
  F1[ii,3]=2*(precision.binary*recall.binary)/(precision.binary+recall.binary)#F1.Myofibroblasts score
  #for pDCs
  tp=length(which(cell_factor=="pDCs"&cell.pred=="pDCs"))#true positive
  fp=length(which(cell_factor!="pDCs"&cell.pred=="pDCs"))#false positive
  fn=length(which(cell_factor=="pDCs"&cell.pred!="pDCs"))#false negative
  precision.binary=tp/(tp+fp)#precision.binary
  recall.binary=tp/(tp+fn)#recall.binary
  F1[ii,4]=2*(precision.binary*recall.binary)/(precision.binary+recall.binary)#F1.pDCs score
  #for Smooth_Muscle_Cells
  tp=length(which(cell_factor=="Smooth_Muscle_Cells"&cell.pred=="Smooth_Muscle_Cells"))#true positive
  fp=length(which(cell_factor!="Smooth_Muscle_Cells"&cell.pred=="Smooth_Muscle_Cells"))#false positive
  fn=length(which(cell_factor=="Smooth_Muscle_Cells"&cell.pred!="Smooth_Muscle_Cells"))#false negative
  precision.binary=tp/(tp+fp)#precision.binary
  recall.binary=tp/(tp+fn)#recall.binary
  F1[ii,5]=2*(precision.binary*recall.binary)/(precision.binary+recall.binary)#F1.Smooth_Muscle_Cells score
}

#save the F1 result in a csv file 
fwrite(F1,file="F1Xgboost_tuned.csv")

##read in the elastic net F1 scores
compF1 <- fread("F1.binary.csv")
compF1[,V1 := NULL]

##make difference plot
diff <- F1 - compF1
diff2 <- melt(diff)
setnames(diff2, c("Celltype","F1_Difference"))
ggplot(diff2, aes(x = Celltype, y = F1_Difference, fill = Celltype))+
  geom_boxplot() +
  geom_abline(slope = 0, intercept = 0, lty = 2)+
  theme_light()+
  theme(legend.position = "n")+
  labs(y = "F1 Difference") +
  ggtitle("Fig 1: XGBoost F1 - Elastic Net F1")

##make XGBoost F1 plot
F1_2<-melt(F1)
compF1_2<-melt(compF1)
setnames(F1_2, c("Celltype","F1"))
ggplot(F1_2, aes(x = Celltype, y = F1, fill = Celltype))+
  geom_boxplot() +
  geom_abline(slope = 0, intercept = 0, lty = 2)+
  theme_light()+
  theme(legend.position = "n")+
  ggtitle("Fig 2: F1 Scores for XGBoost")

##make elastic nett F1 plot
setnames(compF1_2, c("Celltype","F1"))
ggplot(compF1_2, aes(x = Celltype, y = F1, fill = Celltype))+
  geom_boxplot() +
  geom_abline(slope = 0, intercept = 0, lty = 2)+
  theme_light()+
  theme(legend.position = "n")+
  ggtitle("Fig 3: F1 Scores for Elastic Net")

##make overall boxplot
enet <- melt(compF1)
xgb <- melt(F1)
compdat <- data.table(Enet = enet$value, XGboost = xgb$value)
#colMeans(compdat)
compdat <- melt(compdat)
setnames(compdat,c("Model","F1"))
ggplot(compdat, aes(x = Model, y = F1, fill = Model))+
  geom_boxplot() +
  theme_light()+
  theme(legend.position = "n")+
  ggtitle("Fig 4: All F1 values")

##make table and run wilcox test for each celltype
tab1 <- data.table(Celltype = enet$variable,
                   el_net = enet$value,
                   xgboost = xgb$value)
tab1 <- melt(tab1, id.vars = "Celltype")
tab2 <- dcast(tab1,Celltype ~ variable, fun.aggregate = mean)
tab2[,Wilcox_pval := NA_real_]

for(cell in tab2$Celltype){
  wil.test <- wilcox.test(compF1[[cell]],F1[[cell]],paired = T)
  tab2[Celltype == cell, Wilcox_pval := wil.test$p.value]
}

tab2
##run overall wilcox test
compdat <- data.table(Enet = enet$value, XGboost = xgb$value)
wilcox.test(compdat$Enet,compdat$XGboost,paired = T)
###end of script

##note that we also tried feature selection, but abandoned it because
##it didn't improve performance

##feature selection
# response <- dat.Smooth_Muscle_Cells$yy
# xmat <- dat.Smooth_Muscle_Cells[,-1]
# xmat <- as.matrix(xmat)
# 
# smallMod <- xgboost(data = xmat, label = response, learning_rate = 0.05,
#                     max.depth = 3, gamma = 0, min_child_weight = 1, nrounds = 328,
#                     objective = "binary:logistic",subsample=0.5,
#                     colsample_bytree=0.5)
# imp <- xgb.importance(model = smallMod)
# smallPred <- dat.Smooth_Muscle_Cells[,c("yy",imp$Feature)]
# smallPred2 <- smallPred[,-1]
# smallpred2 <- model.matrix(~(TPM2+CALD1+HES4+VCAN+S100A10)^2,smallPred2)
# smallpred2 <- smallpred2[,-c(1:6)]
# smallPred <- cbind(smallPred,smallpred2)
# 
# ii = 1
# res <- numeric(5)
# for(jj in 1:5){
#   x.train=as.matrix(smallPred[idx.cv[,ii]!=jj,-1])
#   x.test=as.matrix(smallPred[idx.cv[,ii]==jj,-1])
#   y.train.Smooth_Muscle_Cells=as.matrix(as.matrix(smallPred[idx.cv[,ii]!=jj,1]))
#   ytest=smallPred$yy[idx.cv[,ii]==jj]
#   #y.test.Smooth_Muscle_Cells=as.matrix(dat.Smooth_Muscle_Cells[idx.cv[,ii]==jj,1])
#   model.ela.Smooth_Muscle_Cells=xgboost(data = x.train, label = y.train.Smooth_Muscle_Cells, learning_rate = 0.08,
#                                         max.depth = 3, gamma = 0, min_child_weight = 2, nrounds = 250,
#                                         objective = "binary:logistic",subsample=0.9,
#                                         colsample_bytree=0.9, eta = 0.8,booster = "gbtree")
#   pred=predict(model.ela.Smooth_Muscle_Cells,x.test,type = "response")
#   pred <- fifelse(pred > 0.5,1,0)
#   res[jj] <- F1(ytest,predicted = pred)
# }
# res
# mean(res)
# var.imp <- xgb.importance(model = model.ela.Smooth_Muscle_Cells)