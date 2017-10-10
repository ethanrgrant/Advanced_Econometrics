solve_matrix = matrix(c(1, 2, 3, 4,2,1,2,3,3,2,1,2,4,3,2,1),nrow=4,ncol=4)
solve(solve_matrix, matrix(c(5,6,7,8),nrow=4,ncol=1))

ran_sample = rnorm(20, mean=1, sd=2)
mean(ran_sample)
var(ran_sample)
hist(ran_sample)

my_mean <- function(y){
  L = length(y)
  m=0
  for(i in 1:L){
    m = m+y[i]/L
  }
  return(m)
}

my_mean(ran_sample)

my_mean_vector <- function(y){
  L = length(y)
  y = matrix(y,1,L)
  v = matrix(1,L,1)
  m = y %*% v/L
  return(m)
}

my_mean_vector(ran_sample)
install.packages('foreign')
library(foreign)
setwd('C:/Users/Ethan/Desktop/Columbia/Senior_Fall/Advanced_Econometrics/')

data=read.dta('TeachingRatings.dta')
head(data)
Y = data$course_eval
X= cbind(1,data$beauty,data$age)
colnames(X)=c('intercept','beauty','age')
plot(Y, X[,'beauty'])
beta_hat <- ginv(t(X)%*%X)%*%t(X)%*%Y
residual <- Y - X%*%beta_hat
mean(residual)
var(residual)
#test <- lm(Y ~ X[,'beauty']+X[,'age'])
M <- diag(length(Y))-X%*%ginv(t(X)%*%X)%*%t(X)

#iii
var(Y)
SS_res <- sum(residual^2)
SS_total <- sum((y-mean(y)^2))
r_squared <- 1-SS_res/SS_total


#iv
X_age <- cbind(1, data$age)
beta_hat_age <- ginv(t(X_age)%*%X_age)%*%t(X_age)%*%Y
residuals_age <- Y-X_age%*%beta_hat_age
beta_hat_beauty <- ginv(t(X_age)%*%X_age)%*%t(X_age)%*%X[,'beauty']
residsuals_beauty <- X[,'beauty']-X_age%*%beta_hat_beauty
beta_2 <- ginv(t(residsuals_beauty)%*%residsuals_beauty)%*%t(residsuals_beauty)%*%residuals_age
#sign indicates that greater age increases rating by a smalla mount for each extra year of age

#v
Y = data$course_eval
X1 = cbind(1, data$age)
X0 = data$beauty
I = diag(length(y))
M1 = I - X1%*%ginv(t(X1)%*%X1)%*%t(X1)
beta = ginv(t(M1%*%X0)%*%M1%*%X0)%*%t(M1%*%X0)%*%M1%*%Y
