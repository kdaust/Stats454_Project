#clean
rm(list = ls())

start_time=Sys.time()
#library package
library(MTPS)
require("xgboost")
#read data
data.list=readRDS("final556.rds")

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

##########find the fix lambda and alpha for binary elastic net for each cell type ##########
##tune max number of iterations
response <- dat.Smooth_Muscle_Cells$yy
xmat <- dat.Smooth_Muscle_Cells[,-1]
xmat <- as.matrix(xmat)
init.parms <- list(learning_rate =0.1,
                   n_estimators=1000,
                   max_depth=8,
                   min_child_weight=1,
                   gamma=0.1,
                   subsample=0.8,
                   colsample_bytree=0.8,
                   objective= 'binary:logistic',
                   scale_pos_weight=1)
cvres <- xgb.cv(init.parms,data = xmat, label = response,
                nrounds = 500, nfold = 5, metrics = "error",
                early_stopping_rounds = 50, print_every_n = 10)
##104 iterations
plist <- list(max_depth = seq(5,25,by = 5),
              min_child_weight = seq(1,6,2))

library(superml)
library(data.table)
xg <- XGBTrainer$new(learning_rate =0.1,
                     n_estimators=104,
                     gamma=0.1,
                     subsample=0.8,
                     colsample_bytree=0.8,
                     objective= 'binary:logistic')
gst <- GridSearchCV$new(trainer = xg,
                        parameters = plist,
                        n_folds = 3,
                        scoring = "f1")
gst$fit(dat.Smooth_Muscle_Cells,"yy")

plist2 <- list(max_depth = c(4,5,6))
gst2 <- GridSearchCV$new(trainer = xg,
                         parameters = plist2,
                         n_folds = 3,
                         scoring = "f1")
gst2$fit(dat.Smooth_Muscle_Cells,"yy")
gst2$best_iteration()

xg$max_depth = 4
xg$min_child_weight = 1

gst3 <- GridSearchCV$new(trainer = xg,
                                 parameters = list(
                                   gamma = c(0,0.1,0.2,0.3,0.4,0.5)),
                                 n_folds = 3,
                                 scoring = "f1")
gst3$fit(dat.Smooth_Muscle_Cells,"yy")
gst3$best_iteration()

##max_depth = 4, min_child_weight = 1, gamma = 0
xgFinal <- XGBTrainer$new(learning_rate = 0.1,
                          n_estimators = 104,
                          gamma=0,
                          max_depth = 4,
                          min_child_weight = 1,
                          subsample=0.8,
                          colsample_bytree=0.8,
                          objective= 'binary:logistic')

library(caret)
split <- createDataPartition(y = response, p = 0.8)
xTrain <- dat.Smooth_Muscle_Cells[split$Resample1,]
xTest <- dat.Smooth_Muscle_Cells[-split$Resample1,]

xgFinal$fit(X = xTrain, y = "yy")
test1 <- xgFinal$predict(df = xTest[,-1])
test2 <- fifelse(test1 < 0.5,0,1)
F1(xTest$yy,test2)
##hyperparameter tuning for xgboost

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

##########fit model##########
for(ii in 1:nrep){
  for(jj in 1:nfolds){
    x.train=as.matrix(dat.B_Cells[idx.cv[,ii]!=jj,-1])
    x.test=as.matrix(dat.B_Cells[idx.cv[,ii]==jj,-1])
    #B_Cells
    y.train.B_Cells=as.matrix(as.matrix(dat.B_Cells[idx.cv[,ii]!=jj,1]))
    #y.test.B_Cells=as.matrix(dat.B_Cells[idx.cv[,ii]==jj,1])
    model.ela.B_Cells=xgboost(data = x.train, label = y.train.B_Cells, max.depth = 20, eta = 1, nrounds = 20, nthread = 2, objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,1]=predict(model.ela.B_Cells,x.test,type = "response")
    #Mesothelial_Cells
    y.train.Mesothelial_Cells=as.matrix(as.matrix(dat.Mesothelial_Cells[idx.cv[,ii]!=jj,1]))
    #y.test.Mesothelial_Cells=as.matrix(dat.Mesothelial_Cells[idx.cv[,ii]==jj,1])
    model.ela.Mesothelial_Cells=xgboost(data = x.train, label = y.train.Mesothelial_Cells, max.depth = 20, eta = 1, nrounds = 20, nthread = 2, objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,2]=predict(model.ela.Mesothelial_Cells,x.test,type = "response")
    #Myofibroblasts
    y.train.Myofibroblasts=as.matrix(as.matrix(dat.Myofibroblasts[idx.cv[,ii]!=jj,1]))
    #y.test.Myofibroblasts=as.matrix(dat.Myofibroblasts[idx.cv[,ii]==jj,1])
    model.ela.Myofibroblasts=xgboost(data = x.train, label = y.train.Myofibroblasts, max.depth = 20, eta = 1, nrounds = 20, nthread = 2, objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,3]=predict(model.ela.Myofibroblasts,x.test,type = "response")
    #pDCs
    y.train.pDCs=as.matrix(as.matrix(dat.pDCs[idx.cv[,ii]!=jj,1]))
    #y.test.pDCs=as.matrix(dat.pDCs[idx.cv[,ii]==jj,1])
    model.ela.pDCs=xgboost(data = x.train, label = y.train.pDCs, max.depth = 20, eta = 1, nrounds = 20, nthread = 2, objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,4]=predict(model.ela.pDCs,x.test,type = "response")
    #Smooth_Muscle_Cells
    y.train.Smooth_Muscle_Cells=as.matrix(as.matrix(dat.Smooth_Muscle_Cells[idx.cv[,ii]!=jj,1]))
    #y.test.Smooth_Muscle_Cells=as.matrix(dat.Smooth_Muscle_Cells[idx.cv[,ii]==jj,1])
    model.ela.Smooth_Muscle_Cells=xgboost(data = x.train, label = y.train.Smooth_Muscle_Cells, max.depth = 20, eta = 1, nrounds = 20, nthread = 2, objective = "binary:logistic")
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
write.csv(F1,file="F1.xgboost_3.csv")

end_time=Sys.time()

#length to run all of code
start_time-end_time