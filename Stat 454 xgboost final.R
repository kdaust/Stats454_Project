#clean
rm(list = ls())

start_time=Sys.time()
#library package
library(MTPS)
require("xgboost")
#read data
data.list=readRDS("final556.rds")


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

# #B_Cells
# xx=as.matrix(dat.B_Cells[,-1])
# yy=as.matrix(dat.B_Cells[,1])
# cv.B_Cells=cv.glmnet2(xx,yy,foldid=NULL,alpha=seq(1,9,by=1)/10,family="binomial",type.measure="class")
# lamdba.fix.B_Cells=cv.B_Cells$lambda.1se
# alpha.fix.B_Cells=cv.B_Cells$alpha
# 
# #Mesothelial_Cells
# xx=as.matrix(dat.Mesothelial_Cells[,-1])
# yy=as.matrix(dat.Mesothelial_Cells[,1])
# #cv.Mesothelial_Cells=cv.glmnet2(xx,yy,foldid=NULL,alpha=seq(1,9,by=1)/10,family="binomial",type.measure="class")
# #lamdba.fix.Mesothelial_Cells=cv.Mesothelial_Cells$lambda.1se
# #alpha.fix.Mesothelial_Cells=cv.Mesothelial_Cells$alpha
# 
# #Myofibroblasts
# xx=as.matrix(dat.Myofibroblasts[,-1])
# yy=as.matrix(dat.Myofibroblasts[,1])
# # cv.Myofibroblasts=cv.glmnet2(xx,yy,foldid=NULL,alpha=seq(1,9,by=1)/10,family="binomial",type.measure="class")
# # lamdba.fix.Myofibroblasts=cv.Myofibroblasts$lambda.1se
# # alpha.fix.Myofibroblasts=cv.Myofibroblasts$alpha
# 
# #pDCs
# xx=as.matrix(dat.pDCs[,-1])
# yy=as.matrix(dat.pDCs[,1])
# cv.pDCs=cv.glmnet2(xx,yy,foldid=NULL,alpha=seq(1,9,by=1)/10,family="binomial",type.measure="class")
# lamdba.fix.pDCs=cv.pDCs$lambda.1se
# alpha.fix.pDCs=cv.pDCs$alpha
# 
# #Smooth_Muscle_Cells
# xx=as.matrix(dat.Smooth_Muscle_Cells[,-1])
# yy=as.matrix(dat.Smooth_Muscle_Cells[,1])
# cv.Smooth_Muscle_Cells=cv.glmnet2(xx,yy,foldid=NULL,alpha=seq(1,9,by=1)/10,family="binomial",type.measure="class")
# lamdba.fix.Smooth_Muscle_Cells=cv.Smooth_Muscle_Cells$lambda.1se
# alpha.fix.Smooth_Muscle_Cells=cv.Smooth_Muscle_Cells$alpha

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
    model.ela.B_Cells=xgboost(data = x.train, label = y.train.B_Cells, max.depth = 20, eta = 1, nrounds = 2, nthread = 2, objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,1]=predict(model.ela.B_Cells,x.test,type = "response")
    #Mesothelial_Cells
    y.train.Mesothelial_Cells=as.matrix(as.matrix(dat.Mesothelial_Cells[idx.cv[,ii]!=jj,1]))
    #y.test.Mesothelial_Cells=as.matrix(dat.Mesothelial_Cells[idx.cv[,ii]==jj,1])
    model.ela.Mesothelial_Cells=xgboost(data = x.train, label = y.train.Mesothelial_Cells, max.depth = 20, eta = 1, nrounds = 2, nthread = 2, objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,2]=predict(model.ela.Mesothelial_Cells,x.test,type = "response")
    #Myofibroblasts
    y.train.Myofibroblasts=as.matrix(as.matrix(dat.Myofibroblasts[idx.cv[,ii]!=jj,1]))
    #y.test.Myofibroblasts=as.matrix(dat.Myofibroblasts[idx.cv[,ii]==jj,1])
    model.ela.Myofibroblasts=xgboost(data = x.train, label = y.train.Myofibroblasts, max.depth = 20, eta = 1, nrounds = 2, nthread = 2, objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,3]=predict(model.ela.Myofibroblasts,x.test,type = "response")
    #pDCs
    y.train.pDCs=as.matrix(as.matrix(dat.pDCs[idx.cv[,ii]!=jj,1]))
    #y.test.pDCs=as.matrix(dat.pDCs[idx.cv[,ii]==jj,1])
    model.ela.pDCs=xgboost(data = x.train, label = y.train.pDCs, max.depth = 20, eta = 1, nrounds = 2, nthread = 2, objective = "binary:logistic")
    y.pred[idx.cv[,ii]==jj,4]=predict(model.ela.pDCs,x.test,type = "response")
    #Smooth_Muscle_Cells
    y.train.Smooth_Muscle_Cells=as.matrix(as.matrix(dat.Smooth_Muscle_Cells[idx.cv[,ii]!=jj,1]))
    #y.test.Smooth_Muscle_Cells=as.matrix(dat.Smooth_Muscle_Cells[idx.cv[,ii]==jj,1])
    model.ela.Smooth_Muscle_Cells=xgboost(data = x.train, label = y.train.Smooth_Muscle_Cells, max.depth = 20, eta = 1, nrounds = 2, nthread = 2, objective = "binary:logistic")
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
write.csv(F1,file="F1.xgboost_2.csv")

end_time=Sys.time()

#length to run all of code
start_time-end_time