n = 1000
B = 1000
VarianceRatio <- function(rho){
  e = rep(NA, n)
  b1 = rep(NA,B) #stores OLS
  b2 = rep(NA, B) #stores F-GLS
  for (b in 1:B){
    e[1] = rnorm(1)
    for(t in 2:n){
      e[t] = rho*e[t-1]+rnorm(1)
    }
    x = rnorm(n) #simulate x
    y = x+e #simulate y
    b1[b] = solve(t(x)%*%x)%*%t(x)%*%y
    e_hat = y - b1[b]%*%x
    rho_hat = solve(t(e_hat[1:n-1])%*%e_hat[1:n-1])%*%t(e_hat[1:n-1])%*%e_hat[2:n]
    y_ = y[2:n] - rho_hat*y[1:(n-1)] #computes quasi difference
    x_ = x[2:n] - rho_hat*x[1:(n-1)] 
    b2[b] = solve(t(x_)%*%x_)%*%t(x_)%*%y_
  }
  return(var(b1)/var(b2))
}
sapply(c(0,0.1,0.5,0.8,0.9), VarianceRatio)

