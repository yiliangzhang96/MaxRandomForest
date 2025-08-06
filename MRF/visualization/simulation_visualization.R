

# R script for simulations of Max-Ranodm Forest

# marginal density


setwd("/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/visualization")

data = unlist(read.table("data.txt",stringsAsFactors = FALSE))
data_true = unlist(read.table("pdf_true.txt",stringsAsFactors = FALSE))
data_estimate = unlist(read.table("pdf_estimate.txt",stringsAsFactors = FALSE))

plot(data_true~data,type = "l")
lines(data_estimate~data,type = "l")








# time for other CDE methods







# KDE

bw <- npcdensbw(formula = data~data)






data("Italy")
attach(Italy)
# First, compute the bandwidths... note that this may take a minute or
# two depending on the speed of your computer.
bw <- npcdensbw(formula=gdp~ordered(year))
# Next, compute the condensity object...
fhat <- npcdens(bws=bw)
# The object fhat now contains results such as the estimated conditional
# density function (fhat$condens) and so on...
summary(fhat)
# Call the plot() function to visualize the results (<ctrl>-C will
# interrupt on *NIX systems, <esc> will interrupt on MS Windows
# systems).
plot(bw)
detach(Italy)
# EXAMPLE 1 (INTERFACE=DATA FRAME): For this example, we load Giovanni
# Baiocchi's Italian GDP panel (see Italy for details), and compute the
# likelihood cross-validated bandwidths (default) using a second-order
# Gaussian kernel (default). Note - this may take a minute or two
# depending on the speed of your computer.
data("Italy")
attach(Italy)
# First, compute the bandwidths... note that this may take a minute or
# two depending on the speed of your computer.
# Note - we cast 􏰀X' and 􏰀y' as data frames so that plot() can
# automatically grab names (this looks like overkill, but in
# multivariate settings you would do this anyway, so may as well get in
# the habit).
X <- data.frame(year=ordered(year))
y <- data.frame(gdp)
bw <- npcdensbw(xdat=X, ydat=y)
# Next, compute the condensity object...
fhat <- npcdens(bws=bw)
# The object fhat now contains results such as the estimated conditional
# density function (fhat$condens) and so on...
summary(fhat)
# Call the plot() function to visualize the results (<ctrl>-C will
# interrupt on *NIX systems, <esc> will interrupt on MS Windows systems).
plot(bw)
detach(Italy)









# see real data
geo_data = data.frame()
geo_2 = data.frame()

geo_raw_data = unlist(read.table("/Users/zhangyiliang/Desktop/Geographical Original of Music/default_features_1059_tracks.txt",stringsAsFactors = FALSE))
geo_data_2 = unlist(read.table("/Users/zhangyiliang/Desktop/Geographical Original of Music/default_plus_chromatic_features_1059_tracks.txt",stringsAsFactors = FALSE))
#strsplit(geo_data[1], ",")

for(i in 1:1059){
  for(j in 1:70){
    geo_data[i,j] = as.numeric(unlist(strsplit(geo_raw_data[i], ",")))[j]
  }
}

rand2 = sample(1:10000,100)


mse_xgb = c()
mse_rf = c()
mse_rf_half = c()
mse_rf_recover = c()
mse_rf_mrf = c()
mse_ls = c()
mse_lasso = c()
mse_ridge = c()
variance = c()

for(i in 1:10){
set.seed(rand2[i])
# select 759 data as training sample, rest to be test sample
test_index = sample(1:1059, 300)

geo_data_test = as.matrix(geo_data[test_index,])
geo_data_training = as.matrix(geo_data[-test_index,])

x_train = geo_data_training[,1:68]
x_test = geo_data_test[,1:68]


#pca_train <- prcomp(x_train_p, center = TRUE,scale. = TRUE, rank = 50)

#x_train = x_train_p %*% pca_train$rotation
#x_test = x_test_p %*% pca_train$rotation


y_train = geo_data_training[,70];
y_test = geo_data_test[,70];


#ls to drop out outliers

#fit <- lsfit(x = x_train, y = y_train)
#box = boxplot(fit$residuals)

#x_train = x_train[-which(fit$residuals %in% box$out),]
#y_train = y_train[-which(fit$residuals %in% box$out)]





x_train_half = x_train[1:380,]
y_train_half = y_train[1:380]

x_train_recover = x_train
y_train_recover = y_train

n_lost = 20;

for(k in 1:(nrow(x_train)-380)){
  for(j in 1:n_lost){
    x_train_recover[k+380,j] = mean(x_train_recover[,j])#runif(1, min(x_train[,j]), max(x_train[,j]))
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

mse_xgb[i] = mean((pred - y_test)^2)


#random forest
library(randomForest)

rf <- randomForest(x = x_train, y = y_train, importance = T)
pred <- predict(rf, x_test)
mse_rf[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train_half, y = y_train_half, importance = T)
pred <- predict(rf, x_test)
mse_rf_half[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train_recover, y = y_train_recover, importance = T)
pred <- predict(rf, x_test)
mse_rf_recover[i] = mean((pred - y_test)^2)


#ls

fit <- lsfit(x = x_train, y = y_train)
pred <- x_test %*% fit$coefficients[-1] + fit$coefficients[1]
mse_ls[i] = mean((pred - y_test)^2)

#lasso, Ridge and E-Net

lambda.max <- cv.glmnet(x_train, y_train, family = "gaussian", type.measure="mse", alpha = 1)$lambda[1]
lambda.min <- lambda.max*0.1
lambda <- exp(seq(log(lambda.max), log(lambda.min), length = 1000))


library(glmnet)
fit_lasso <- cv.glmnet(x_train, y_train, alpha = 1, family = "gaussian", lambda = lambda/20)

pred = predict(fit_lasso, newx = x_test, s = "lambda.1se")        
mse_lasso[i] = mean((pred - y_test)^2)  



fit_ridge <- cv.glmnet(x_train, y_train, alpha = 0, family = "gaussian", lambda = lambda/1)

pred = predict(fit_ridge, newx = x_test, s = "lambda.1se")        
mse_ridge[i] = mean((pred - y_test)^2)  

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
  test_index = sample(1:1059, 300)
  
  geo_data_test = as.matrix(geo_data[test_index,])
  geo_data_training = as.matrix(geo_data[-test_index,])
  
  x_train = geo_data_training[,1:68];
  y_train = geo_data_training[,70];
  
  x_test = geo_data_test[,1:68];
  y_test = geo_data_test[,70];
  
  # change the column (group the correlated features)
  #for(j in 1:17){
  #  for(k in 1:4){
  #    x_train[,4*(j-1) + k] = geo_data_training[,17*(k-1) + j];
  #    x_test[,4*(j-1) + k] = geo_data_test[,17*(k-1) + j];
  #  }
  #}
  #x_train = geo_data_training[,1:68];
  #y_train = geo_data_training[,69];
  
  #x_test = geo_data_test[,1:68];
  #y_test = geo_data_test[,69];
  
  x_train_half = x_train[1:380,]
  y_train_half = y_train[1:380]
  
  x_train_recover = x_train
  y_train_recover = y_train
  
  n_lost = 10;
  
  for(k in 1:379){
    for(j in 1:n_lost){
      x_train_recover[k+380,j] = mean(x_train_recover[,j])#runif(1, min(x_train[,j]), max(x_train[,j]))
    }
  }

z_train = data.frame(x_train_half[,11:68],y_train_half,x_train_half[,1:10])
z_test = data.frame(x_train[381:759,11:68],y_train[381:759])
  
write.table(z_train,file = "/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/training_data.txt",row.names = F,col.names = F)
write.table(z_test,file = "/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/test_data.txt",row.names = F,col.names = F)




geo_mrf_data = read.table("/Users/zhangyiliang/Desktop/Max random Forest/Final/Data/conditional_density/regression/recovered_data.txt",stringsAsFactors = FALSE)

x_train_mrf = x_train
y_train_mrf = y_train

for(k in 1:378){
  for(j in 1:n_lost){
    x_train_mrf[k+380,j] = geo_mrf_data[k,59+j]
  }
}

rf <- randomForest(x = x_train_mrf, y = y_train_mrf, importance = T)
pred <- predict(rf, x_test)
mse_rf_mrf[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train, y = y_train, importance = T)
pred <- predict(rf, x_test)
mse_rf[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train_half, y = y_train_half, importance = T)
pred <- predict(rf, x_test)
mse_rf_half[i] = mean((pred - y_test)^2)

rf <- randomForest(x = x_train_recover, y = y_train_recover, importance = T)
pred <- predict(rf, x_test)
mse_rf_recover[i] = mean((pred - y_test)^2)


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






