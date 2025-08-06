

cond_data = data.frame()
#cond_2 = data.frame()

cond_raw_data = read.csv("/Users/zhangyiliang/Desktop/superconduct/train.csv",stringsAsFactors = FALSE)
#cond_data_2 = unlist(read.table("/Users/zhangyiliang/Desktop/Geographical Original of Music/default_plus_chromatic_features_1059_tracks.txt",stringsAsFactors = FALSE))
#strsplit(geo_data[1], ",")

cond_data = cond_raw_data

cond_data[,1:81] <- scale(cond_data[,1:81])
#rand2 = sample(1:10000,100)


error_xgb = c()
error_rf = c()
error_rf_half = c()
error_rf_recover = c()
error_rf_mrf = c()
error_ls = c()
error_lasso = c()
error_ridge = c()
error_rf_mice = c()
error_rf_missForest = c()
error_rf_Hmisc = c()
error_rf_mi = c()
variance = c()


set.seed(888)
rand2 = sample(1:10000,100)

for(i in 1:10){
set.seed(rand2[i])
# select 759 data as training sample, rest to be test sample
test_index = sample(1:21263, 3263)

cond_data_test = cond_data[test_index,]
cond_data_training = cond_data[-test_index,]
cond_data_training = cond_data_training[sample(nrow(cond_data_training)),]


x_train = as.matrix(cond_data_training[,1:81])
x_test = as.matrix(cond_data_test[,1:81])


#pca_train <- prcomp(x_train_p, center = TRUE,scale. = TRUE, rank = 50)

#x_train = x_train_p %*% pca_train$rotation
#x_test = x_test_p %*% pca_train$rotation



cor_XY <- abs(cor(cond_data_training[,1:81], cond_data_training[,82]))
for(j in 1:81){
  
  x_train[,rank(-cor_XY)[j]] = cond_data_training[,j]
  x_test[,rank(-cor_XY)[j]] = cond_data_test[,j]
  
}

y_train = cond_data_training[,82];
y_test = cond_data_test[,82];


#ls to drop out outliers

#fit <- lsfit(x = x_train, y = y_train)
#box = boxplot(fit$residuals)

#x_train = x_train[-which(fit$residuals %in% box$out),]
#y_train = y_train[-which(fit$residuals %in% box$out)]





x_train_half = x_train[1:5000,]
y_train_half = y_train[1:5000]

x_train_recover = x_train
y_train_recover = y_train

n_lost = 50;

for(j in 1:n_lost){
  aa = mean(x_train_recover[,j])
  for(k in 1:(nrow(x_train)-5000)){
    x_train_recover[k+5000,j] = aa#runif(1, min(x_train[,j]), max(x_train[,j]))
  }
}

#x_train_recover = x_train_p %*% pca_train$rotation






# regression for two coordinates respectively


# xgboost
library(xgboost)
# para tuned
xgb <- xgboost(data = x_train, label = y_train, max.depth = 2, eta = 1, nthread = 2, nrounds = 2)
#or
#xgb_cv <- xgb.cv(data = x_train, label = y_train, max.depth = 2, eta = 1, nthread = 2, nrounds = 2, nfold = 10)

pred <- predict(xgb, x_test)
#pred <- predict(xgb_cv, x_test)

error_xgb[i] = mean((pred - y_test)^2)


#random forest
library(randomForest)

rf <- randomForest(x = x_train, y = y_train, importance = T)
pred <- predict(rf, x_test)
error_rf[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train_half, y = y_train_half, importance = T)
pred <- predict(rf, x_test)
error_rf_half[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train_recover, y = y_train_recover, importance = T)
pred <- predict(rf, x_test)
error_rf_recover[i] = mean((pred - y_test)^2)


#logistic

model <- nnet::multinom(V55 ~., data = cond_data_training)
pred <- predict(model,newdata = cond_data_test)
error_ls[i] = mean((pred - y_test)^2)

#lasso, Ridge and E-Net

lambda.max <- cv.glmnet(x_train, y_train, family = "gaussian", type.measure="mse", alpha = 1)$lambda[1]
lambda.min <- lambda.max*0.1
lambda <- exp(seq(log(lambda.max), log(lambda.min), length = 1000))


library(glmnet)
fit_lasso <- cv.glmnet(x_train, y_train, alpha = 1, family = "gaussian", lambda = lambda/20)

pred = predict(fit_lasso, newx = x_test,type = "lambda.1se")        
error_lasso[i] = mean((pred - y_test)^2)



fit_ridge <- cv.glmnet(x_train, y_train, alpha = 0, family = "gaussian", lambda = lambda/1)

pred = predict(fit_ridge, newx = x_test, s = "lambda.1se")        
error_ridge[i] = mean((pred - y_test)^2)

variance[i] = var(y_test)


}

# calculate the mse
mean(mse_xgb)
mean(mse_rf)
mean(mse_rf_half)
mean(mse_rf_recover)
mean(mse_ls)
mean(mse_lasso)
mean(mse_ridge)
mean(variance)






























# transport the data outside and use MRF

for(i in 1:20){
  set.seed(rand2[i])

  test_index = sample(1:21263, 3263)
  
  cond_data_test = cond_data[test_index,]
  cond_data_training = cond_data[-test_index,]
  cond_data_training = cond_data_training[sample(nrow(cond_data_training)),]
  
  
  x_train = as.matrix(cond_data_training[,1:81])
  x_test = as.matrix(cond_data_test[,1:81])
  
  
  #aggregation
  cor_XY <- abs(cor(cond_data_training[,1:81], cond_data_training[,82]))
  for(j in 1:81){
    
    x_train[,rank(-cor_XY)[j]] = cond_data_training[,j]
    x_test[,rank(-cor_XY)[j]] = cond_data_test[,j]
    
  }
  
  y_train = cond_data_training[,82];
  y_test = cond_data_test[,82];

  # use logistic to produce missing value (MAR)
  missing = c();
  alpha = c(0,0.1,-0.1,0.1,-0.1,0.2,-0.2,0.2,-0.2,0.3,-0.3)#c(0,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2)
  
  for(j in 1:nrow(x_train)){
    missing[j] = rbern(1,exp(alpha%*%c(1,x_train[j,11:20]))/(1+exp(alpha%*%c(1,x_train[j,11:20]))))
  }

  
  
  x_train_half = x_train[which(missing == 1),]
  y_train_half = y_train[which(missing == 1)]
  
  #generate NA in missing data
  n_lost = 10;
  
  x_train_missing = x_train
  x_train_missing[which(missing == 0),1:n_lost] = NA
  
  x_train_recover = x_train_missing
  y_train_recover = y_train  

  #x_train_missing = x_train_missing[c(which(missing == 0),which(missing == 1)),]

  
  for(j in 1:n_lost){
    aa = mean(x_train[which(missing == 1),j])
    x_train_recover[which(missing == 0),j] = aa #runif(1, min(x_train[,j]), max(x_train[,j]))
  }

  
    

z_train = data.frame(x_train[which(missing == 1),11:21],y_train[which(missing == 1)],x_train[which(missing == 1),1:10]) #,y_train_half,y_train[20001:40000]
z_test = data.frame(x_train[which(missing == 0),11:21],y_train[which(missing == 0)])
  
write.table(z_train,file = "/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/training_data.txt",row.names = F,col.names = F)
write.table(z_test,file = "/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/test_data.txt",row.names = F,col.names = F)




cond_mrf_data = read.table("/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/recovered_data.txt",stringsAsFactors = FALSE)

x_train_mrf = x_train
y_train_mrf = y_train

for(k in 1:length(which(missing == 0))){
  for(j in 1:n_lost){
    x_train_mrf[which(missing == 0)[k],j] = cond_mrf_data[k,12+j]
  }
}


#rf <- randomForest(x = x_train_mrf[1:17999,], y = y_train_mrf[1:17999], importance = T)
#pred <- predict(rf, x_test)
#error_rf_mrf[i] = mean((pred - y_test)^2)
apply(x_train[,1:10],2,mean)
apply(x_train_half[,1:10],2,mean)
apply(x_train_recover[,1:10],2,mean)
apply(x_train_mrf[,1:10],2,mean)


rf <- randomForest(x = x_train, y = y_train, importance = T)
pred <- predict(rf, x_test)
error_rf[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train_half, y = y_train_half, importance = T)
pred <- predict(rf, x_test)
error_rf_half[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train_recover, y = y_train_recover, importance = T)
pred <- predict(rf, x_test)
error_rf_recover[i] = mean((pred - y_test)^2)

#mice
library(mice)
data_train_missing = cbind(x_train_missing,y_train)



#for(j in 1:n_lost){
#  for(k in 1:(nrow(x_train)-20000)){
#    data_train_missing[k+20000,j] = NA#runif(1, min(x_train[,j]), max(x_train[,j]))
#  }
#}
#colnames(data_train_missing)[55:61] = c("V55","V56","V57","V58","V59","V60","V61")

data_train_mice <- mice(data_train_missing, m=5, maxit = 5, method = 'pmm', seed = 500)
complete = complete(data_train_mice)
apply(complete[,1:10],2,mean)


x_train_mice = data_train_mice[,1:ncol(x_train)]
y_train_mice = y_train

rf <- randomForest(x = x_train_mice, y = y_train_mice, importance = T)
pred <- predict(rf, x_test)
error_rf_mice[i] = mean((pred - y_test)^2)

  
#missForest
library(missForest)
data_train_missForest = missForest(data_train_missing,maxiter = 5,ntree = 5)
complete_forest = data_train_missForest$ximp
apply(complete_forest[,1:10],2,mean)

x_train_missForest = data_train_missForest[,1:ncol(x_train)]
y_train_missForest = y_train

rf <- randomForest(x = x_train_missForest, y = y_train_missForest, importance = T)
pred <- predict(rf, x_test)
error_rf_missForest[i] = mean((pred - y_test)^2)

#Hmisc
library(Hmisc)



#mi
library(mi)
}






###### Diagnostic ######

# correlation heat map
head(x_train)
cormat <- round(cor(x_train),2)
melted_cormat <- melt(cormat)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# diagnostic of the model

# take log, does not work out
y_train = log(geo_data_training[,69]);
y_test = log(geo_data_test[,69]);

# merge 4 into 1, worse than using 68 features
x_train = geo_data_training[,1:17]
x_test = geo_data_test[,1:17]

for(i in 1:17){
  x_train[,i] = geo_data_training[,i] + geo_data_training[,i+17] + geo_data_training[,i+2*17] + geo_data_training[,i+3*17];
  x_test[,i] = geo_data_test[,i] + geo_data_test[,i+17] + geo_data_test[,i+2*17] + geo_data_test[,i+3*17];
}

# PCA does not work out, drop out outlier does not help either
pca_train <- prcomp(x_train, center = TRUE,scale. = TRUE)

x_train = x_train %*% pca_train$rotation
x_test = x_test %*% pca_train$rotation


# change the column
for(j in 1:17){
  for(k in 1:4){
    x_train[4*(j-1) + k] = geo_data_training[17*(k-1) + j];
    x_test[4*(j-1) + k] = geo_data_test[17*(k-1) + j];
  }
}
#x_train = geo_data_training[,1:68];
y_train = geo_data_training[,69];

#x_test = geo_data_test[,1:68];
y_test = geo_data_test[,69];



# check the collinearity by calculating VIF







# check the correlation, turns out to be stable
cor_XY <- abs(cor(x_train, y_train))
sort(cor_XY, decreasing = T)[1:20]

cor(x_train, y_train)[which(abs(cor(x_train, y_train)) >= sort(cor_XY, decreasing = T)[20])]
cor(x_test, y_test)[which(abs(cor(x_train, y_train)) >= sort(cor_XY, decreasing = T)[20])]
sort(abs(cor(x_test, y_test)), decreasing = T)[1:20]






