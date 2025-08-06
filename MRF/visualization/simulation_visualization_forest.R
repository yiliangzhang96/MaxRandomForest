

forest_data = data.frame()
#forest_2 = data.frame()

forest_raw_data = unlist(read.table("/Users/zhangyiliang/Desktop/covtype.txt",stringsAsFactors = FALSE))
#forest_data_2 = unlist(read.table("/Users/zhangyiliang/Desktop/Geographical Original of Music/default_plus_chromatic_features_1059_tracks.txt",stringsAsFactors = FALSE))
#strsplit(geo_data[1], ",")

for(i in 1:50000){
  for(j in 1:55){
    forest_data[i,j] = as.numeric(unlist(strsplit(forest_raw_data[i], ",")))[j]
  }
}

forest_data[,55] = as.factor(forest_data[,55])
forest_data = forest_data[1:50000,]

forest_data[,1:54] <- scale(forest_data[,1:54])
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

for(i in 1:10){
set.seed(rand2[i])
# select 759 data as training sample, rest to be test sample
test_index = sample(1:50000, 10000)

forest_data_test = forest_data[test_index,]
forest_data_training = forest_data[-test_index,]


x_train = as.matrix(forest_data_training[,1:54])
x_test = as.matrix(forest_data_test[,1:54])


#pca_train <- prcomp(x_train_p, center = TRUE,scale. = TRUE, rank = 50)

#x_train = x_train_p %*% pca_train$rotation
#x_test = x_test_p %*% pca_train$rotation


y_train = forest_data_training[,55];
y_test = forest_data_test[,55];


#ls to drop out outliers

#fit <- lsfit(x = x_train, y = y_train)
#box = boxplot(fit$residuals)

#x_train = x_train[-which(fit$residuals %in% box$out),]
#y_train = y_train[-which(fit$residuals %in% box$out)]





x_train_half = x_train[1:20000,]
y_train_half = y_train[1:20000]

x_train_recover = x_train
y_train_recover = y_train

n_lost = 20;

for(j in 1:n_lost){
  aa = mean(x_train_recover[,j])
  for(k in 1:(nrow(x_train)-20000)){
    x_train_recover[k+20000,j] = aa#runif(1, min(x_train[,j]), max(x_train[,j]))
  }
}

#x_train_recover = x_train_p %*% pca_train$rotation






# regression for two coordinates respectively


# xgboost
library(xgboost)
# para tuned
xgb <- xgboost(data = x_train, label = as.character(y_train), max.depth = 2, eta = 1, nthread = 2, nrounds = 2, objective = "multi:softmax",num_class = 8)
#or
#xgb_cv <- xgb.cv(data = x_train, label = y_train, max.depth = 2, eta = 1, nthread = 2, nrounds = 2, nfold = 10)

pred <- predict(xgb, x_test)
#pred <- predict(xgb_cv, x_test)

error_xgb[i] = length(which(pred != y_test))


#random forest
library(randomForest)

rf <- randomForest(x = x_train, y = y_train, importance = T)
pred <- predict(rf, x_test)
error_rf[i] = length(which(pred != y_test))

rf <- randomForest(x = x_train_half, y = y_train_half, importance = T)
pred <- predict(rf, x_test)
error_rf_half[i] = length(which(pred != y_test))

rf <- randomForest(x = x_train_recover, y = y_train_recover, importance = T)
pred <- predict(rf, x_test)
error_rf_recover[i] = length(which(pred != y_test))


#logistic

model <- nnet::multinom(V55 ~., data = forest_data_training)
pred <- predict(model,newdata = forest_data_test)
error_ls[i] = length(which(pred != y_test))

#lasso, Ridge and E-Net

lambda.max <- 1#cv.glmnet(x_train, y_train, family = "multinomial", type.measure="mse", alpha = 1)$lambda[1]
lambda.min <- lambda.max*0.1
lambda <- exp(seq(log(lambda.max), log(lambda.min), length = 1000))


library(glmnet)
fit_lasso <- cv.glmnet(x_train, y_train, alpha = 1, family = "multinomial", lambda = lambda/20)

pred = predict(fit_lasso, newx = x_test,type = "class")        
error_lasso[i] = length(which(pred != y_test))



fit_ridge <- cv.glmnet(x_train, y_train, alpha = 0, family = "multinomial", lambda = lambda/1)

pred = predict(fit_ridge, newx = x_test, s = "lambda.1se")        
error_ridge[i] = length(which(pred != y_test))

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
  # select 759 data as training sample, rest to be test sample
  test_index = sample(1:50000, 10000)
  
  forest_data_test = forest_data[test_index,]
  forest_data_training = forest_data[-test_index,]
  forest_data_training = forest_data_training[sample(nrow(forest_data_training)),]
  
  x_train = as.matrix(forest_data_training[,1:54])
  x_test = as.matrix(forest_data_test[,1:54])
  
  
  #pca_train <- prcomp(x_train_p, center = TRUE,scale. = TRUE, rank = 50)
  
  #x_train = x_train_p %*% pca_train$rotation
  #x_test = x_test_p %*% pca_train$rotation
  
  
  y_train = forest_data_training[,55];
  y_test = forest_data_test[,55];
  
  # change the column (group the correlated features)
  #for(j in 1:17){
  #  for(k in 1:4){
  #    x_train[,4*(j-1) + k] = geo_data_training[,17*(k-1) + j];
  #    x_test[,4*(j-1) + k] = geo_data_test[,17*(k-1) + j];
  #  }
  #}

  
  x_train_half = x_train[1:20000,]
  y_train_half = y_train[1:20000]
  
  x_train_recover = x_train
  y_train_recover = y_train
  
  n_lost = 10;
  
  for(j in 1:n_lost){
    aa = mean(x_train_recover[,j])
    for(k in 1:(nrow(x_train)-20000)){
      x_train_recover[k+20000,j] = aa#runif(1, min(x_train[,j]), max(x_train[,j]))
    }
  }
  
  # apply one-hot excoding scheme
  y_train_hot = as.data.frame(matrix(0,ncol = 7,nrow = 20000))
  for(j in 1:20000){
    y_train_hot[j,as.numeric(y_train[j])] = 1
  }
  
  y_test_hot = as.data.frame(matrix(0,ncol = 7,nrow = 20000))
  for(j in 1:20000){
    y_test_hot[j,as.numeric(y_train[j+20000])] = 1
  }

z_train = data.frame(y_train_hot,x_train_half[,1:10]) #,y_train_half,y_train[20001:40000]
z_test = data.frame(y_test_hot)
  
write.table(z_train,file = "/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/training_data.txt",row.names = F,col.names = F)
write.table(z_test,file = "/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/test_data.txt",row.names = F,col.names = F)




forest_mrf_data = read.table("/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/recovered_data.txt",stringsAsFactors = FALSE)

x_train_mrf = x_train
y_train_mrf = y_train

for(k in 1:20000){
  for(j in 1:n_lost){
    x_train_mrf[k+20000,j] = forest_mrf_data[k,7+j]
  }
}


rf <- randomForest(x = x_train_mrf[1:39999,], y = y_train_mrf[1:39999], importance = T)
pred <- predict(rf, x_test)
error_rf_mrf[i] = length(which(pred != y_test))

rf <- randomForest(x = x_train, y = y_train, importance = T)
pred <- predict(rf, x_test)
error_rf[i] = length(which(pred != y_test))

rf <- randomForest(x = x_train_half, y = y_train_half, importance = T)
pred <- predict(rf, x_test)
error_rf_half[i] = length(which(pred != y_test))

rf <- randomForest(x = x_train_recover, y = y_train_recover, importance = T)
pred <- predict(rf, x_test)
error_rf_recover[i] = length(which(pred != y_test))

#mice
library(mice)
data_train_missing = cbind(x_train,y_train_hot)



for(j in 1:n_lost){
  for(k in 1:(nrow(x_train)-20000)){
    data_train_missing[k+20000,j] = NA#runif(1, min(x_train[,j]), max(x_train[,j]))
  }
}
colnames(data_train_missing)[55:61] = c("V55","V56","V57","V58","V59","V60","V61")

data_train_mice <- mice(data_train_missing, m=5, maxit = 50, method = 'pmm', seed = 500)

x_train_mice = data_train_mice[,1:ncol(x_train)]
y_train_mice = y_train

rf <- randomForest(x = x_train_mice, y = y_train_mice, importance = T)
pred <- predict(rf, x_test)
error_rf_mice[i] = length(which(pred != y_test))

  
#missForest
library(missForest)
data_train_missForest = missForest(data_train_missing)

x_train_missForest = data_train_missForest[,1:ncol(x_train)]
y_train_missForest = y_train

rf <- randomForest(x = x_train_missForest, y = y_train_missForest, importance = T)
pred <- predict(rf, x_test)
error_rf_missForest[i] = length(which(pred != y_test))

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






