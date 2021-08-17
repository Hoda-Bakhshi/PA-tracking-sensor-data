
load(file="Final_features.Rda")


# 
# Classification model using accelerometer only
# for(i in 1:Nuser){
#   
#   for(j in seq_along(devices))
#   {
#     # remove GPS features
#     Final_features[[i]][[j]]<- Final_features[[i]][[j]][,-c(86,87)]
#     
#   }}

##############
#Total data

cbindlist <- function(list) {
  n <- length(list)
  res <- c(1:nrow(list[[1]]))
  for (i in seq(n)) res <- cbind(res, list[[i]])
  return(res)
}

library(data.table)
trainx=c()
for(i in 1:Nuser){
  trainx[[i]]=cbindlist(Final_features[[i]])
}

for(i in 1:Nuser){
  
  trainx[[i]][["res"]]=NULL
}

for(i in 1:Nuser){
  trainx[[i]]$activity2=NULL
  trainx[[i]]$activity3=NULL
  trainx[[i]]$activity4=NULL
  trainx[[i]]$activity5=NULL
}

#Create training and test datasets for Leave one subject out cross validation
library(cvTools)
folds<-cvFolds(length(trainx),K=Nuser,type = c( "consecutive"))

train_dataset <- as.list(c(1:Nuser))
test_dataset<- as.list(c(1:Nuser))

for(i in 1:Nuser){
  for(j in 1:length(folds$subsets[folds$which != i]))
  {train_index<-folds$subsets[folds$which != i]
  train_dataset[[i]]  <- rbind(train_dataset[[i]],(trainx[[train_index[j]]]))#Set the training set
  }
  test_dataset[[i]]<-trainx[[folds$subsets[folds$which == i]]]
}

for(i in 1:Nuser){
  # e=which(is.na(train_dataset[[i]][["activity"]]))
  # if(length(e)>0)
  {train_dataset[[i]]<-train_dataset[[i]][-1,]}
}


for(i in 1:Nuser){
  train_dataset[[i]]=na.omit(train_dataset[[i]])
  test_dataset[[i]]=na.omit(test_dataset[[i]])
}

for(i in 1:Nuser){
  train_dataset[[i]]$activity1<-as.factor(train_dataset[[i]]$activity1)
  test_dataset[[i]]$activity1<-as.factor(test_dataset[[i]]$activity1)}



#Random forest classification

m<-as.list(c(1:Nuser))

#install.packages("tictoc")

library(tictoc)
tic()
system.time( 
  for(i in 1:Nuser){
    
    library(ranger)
    
    train_dataset[[i]] <- train_dataset[[i]][!is.na(train_dataset[[i]]$activity1), ]
    m[[i]] <- ranger(activity1~.,data=as.data.frame(train_dataset[[i]]) ,num.trees = 200,num.threads=24,mtry=20, importance = "impurity")
    message("This is step:", i)
    
  })
toc()


model_cm=c()
model_acc=c()
Prediction_Result=c()
Prediction_cm=c()
Prediction_acc<-c()
Prediction_Result2=c()
obs=c()
Prediction_acc2=c()
for(i in 1:Nuser)
{library(caret)
  model_cm[[i]] =confusionMatrix(train_dataset[[i]]$activity,m[[i]][["predictions"]])
  model_acc[[i]]=model_cm[[i]][["overall"]][["Accuracy"]]
  
  Prediction_Result[[i]] <- predict(m[[i]],as.data.frame(test_dataset[[i]]),type='response')
  Prediction_cm[[i]] <- confusionMatrix(test_dataset[[i]]$activity,Prediction_Result[[i]][["predictions"]], mode = "everything")
  Prediction_acc[[i]] <- Prediction_cm[[i]][["overall"]][["Accuracy"]]
  
  Prediction_Result2[[i]] <- predict(m[[i]],as.data.frame(test_dataset[[i]]),type='response')
  obs[[i]] <- test_dataset[[i]]$activity
  Prediction_acc2[[i]] <- length(which(obs[[i]] == Prediction_Result2[[i]][["predictions"]])) / length(obs[[i]])
}
kappa=c()
for(i in 1:Nuser)
  kappa=rbind(kappa,Prediction_cm[[i]][["overall"]])



options(digits=2)
mean(model_acc)
mean(Prediction_acc)
mean(Prediction_acc2)

colMeans(kappa)
kappa
#####################
# Gradient boosting  classification

require(xgboost)
#load libraries
require(data.table)

train_dataset2=train_dataset
test_dataset2= test_dataset



labels=c()
ts_label=c()
new_tr=c()
new_ts=c()
dtrain=as.list(c(1:Nuser))
dtest=as.list(c(1:Nuser))
for(i in 1:Nuser){
  train_dataset[[i]]<-data.table(train_dataset[[i]])
  test_dataset[[i]]<-data.table(test_dataset[[i]])
  #using one hot encoding 
  labels[[i]] <- train_dataset[[i]]$activity1
  ts_label[[i]] <- test_dataset[[i]]$activity1
  
  new_tr[[i]] <- model.matrix(~.+0,data = train_dataset[[i]][,-c("activity1"),with=F]) 
  new_ts[[i]] <- model.matrix(~.+0,data = test_dataset[[i]][,-c("activity1"),with=F])
  
  #convert factor to numeric 
  labels [[i]]<- as.numeric(labels[[i]])-1
  ts_label[[i]] <- as.numeric(ts_label[[i]])-1
  
  #preparing matrix 
  dtrain[[i]] <- xgb.DMatrix(data = new_tr[[i]],label = labels[[i]]) 
  dtest[[i]] <- xgb.DMatrix(data = new_ts[[i]],label=ts_label[[i]])
  
}

#default parameters
params <- list(booster = "gbtree", objective = "multi:softmax", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1,num_class =6)


gc()
tic()
m<-as.list(c(1:Nuser))


library(tictoc)
tic()
system.time( 
  for(i in 1:Nuser){
    m[[i]] <-  xgb.train (params = params, data = dtrain[[i]], nrounds = 200)
    message("This is step:", i)
  })
toc()



model_cm=c()
model_acc=c()
Prediction_Result=c()
Prediction_cm=c()
Prediction_acc<-c()
Prediction_Result2=c()
obs=c()
Prediction_acc2=c()
for(i in 1:Nuser)
{library(caret)
  model_cm[[i]] =confusionMatrix(train_dataset[[i]]$activity,m[[i]][["predictions"]])
  model_acc[[i]]=model_cm[[i]][["overall"]][["Accuracy"]]
  
  Prediction_Result[[i]] <- predict(m[[i]],as.data.frame(test_dataset[[i]]),type='response')
  Prediction_cm[[i]] <- confusionMatrix(test_dataset[[i]]$activity,Prediction_Result[[i]][["predictions"]], mode = "everything")
  Prediction_acc[[i]] <- Prediction_cm[[i]][["overall"]][["Accuracy"]]
  
  Prediction_Result2[[i]] <- predict(m[[i]],as.data.frame(test_dataset[[i]]),type='response')
  obs[[i]] <- test_dataset[[i]]$activity
  Prediction_acc2[[i]] <- length(which(obs[[i]] == Prediction_Result2[[i]][["predictions"]])) / length(obs[[i]])
}
kappa=c()
for(i in 1:Nuser)
  kappa=rbind(kappa,Prediction_cm[[i]][["overall"]])



options(digits=2)
mean(model_acc)
mean(Prediction_acc)
mean(Prediction_acc2)

colMeans(kappa)
kappa

