library(e1071)
library(limma)
args <- commandArgs(trailingOnly = TRUE)
print(args)

RMA=read.delim("training.txt", header = TRUE, row.names = 1)
negative=-1.6
myrows=c(6:85)
totalclasses=array(1,dim(RMA)[1])
totalclasses[which(RMA$l2fc<negative)]=-1

predictor<-svm( as.matrix(RMA[,myrows]), as.factor(totalclasses),probability=T)
write.svm(predictor, svm.file="wang.model")
set1=read.delim(args[1], header = TRUE, row.names = 1)
set1.pred <- slot(predict(predictor, as.matrix(set1),probability=T),"probabilities")[,2]
write.table(cbind(rownames(set1),set1.pred),args[2],sep="\t",row.names=F,col.names=F)
