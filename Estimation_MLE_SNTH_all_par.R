library("parallel")
library("mvtnorm")
library("LambertW")

N = 200 # sample size
p = 3 # dimension

# True parameters
Psi_bar = matrix(data = c(1,-0.5,0.3,-0.5,1,-0.2,0.3,-0.2,1),nrow = p) 
eta = c(-1.5,2,0.5)
xi = c(0.8,-0.6,1.3)
omega = diag(c(3,5,2))
h = c(0.02,0.08,0.03)

# Generation from SNTH distribution
X = rmvnorm(n = N,mean = rep(0,p),sigma = Psi_bar)
Z = matrix(nrow = N,ncol = p)
for(i in 1:N){
  Z[i,] = eta*abs(rnorm(1)) + X[i,]
}
tau_h = function(x){
  for(i in 1:p){
    x[i] = x[i]*exp(0.5*h[i]*x[i]^2)
  }
  return(x)
}
Y = matrix(nrow = N,ncol = p)
for(i in 1:N){
  Y[i,] = xi + omega%*%tau_h(Z[i,])
}


# marginal MLE estimation 
xi_est = vector(length = p)
omega_est = vector(length = p)
eta_est = vector(length = p)
h_est = vector(length = p)
tau_h_Z_est = matrix(nrow = N,ncol = p)
Z_est = matrix(nrow = N,ncol = p)
for(k in 1:p){
  marginal_llh = function(theta){
    xi = theta[1]
    omega = exp(theta[2])
    eta = theta[3]
    h = exp(theta[4])
    
    g = function(z){
      return(z*exp(-0.5*W(h*z^2)))
    }
    
    llh = 0
    for (i in 1:N) {
      z = (Y[i,k] - xi)/omega
      
      llh = llh + log(1/omega) + 0.5*W(h*z^2) - log(h*z^2 + exp(W(h*z^2))) + log(2) + log(1/sqrt(2*pi)) - 0.5*log(1+eta^2)- 0.5*g(z)^2/(1+eta^2) + log(pnorm(eta*g(z)/sqrt(1+eta^2)))
    }
    return(llh)
  }
  negative_llh = function(theta){
    return(-marginal_llh(theta))
  }
  H_optim = optim(c(mean(Y[,k]),log(sd(Y[,k])),skewness(Y[,k])/sqrt(1+skewness(Y[,k])^2),log(kurtosis(Y[,k])/sqrt(1+kurtosis(Y[,k])^2))),negative_llh)
  xi_est[k] = H_optim$par[1]
  omega_est[k] = exp(H_optim$par[2])
  eta_est[k] = H_optim$par[3]
  h_est[k] = exp(H_optim$par[4])
  
  for(i in 1:N){
    tau_h_Z_est[i,k] = (Y[i,k] - xi_est[k])/omega_est[k]
  }
  
  tukey_h_univariate = function(z,h){
    return(z*exp(0.5*h*z^2))
  }
  
  for(i in 1:N){
    f = function(x){
      return(tau_h_Z_est[i,k] - tukey_h_univariate(x,h_est[k]))
    }
    Z_est[i,k] = uniroot(f,c(-100,100),tol=10^-9)$root
  }
}

## marginal MLE done ##
## estimates are: xi_est, omega_est, eta_est, h_est

## EM for SN ##

llh_SN = function(Psi_0){
  Omega_0 = Psi_0 + eta_est%*%t(eta_est)
  llh = 0
  for(i in 1:N){
    temp_1 = 0.5*as.numeric(t(Z_est[i,])%*%solve(Omega_0)%*%(Z_est[i,]))
    temp_2 = as.numeric(t(eta_est)%*%solve(Psi_0)%*%Z_est[i,]) / sqrt(1 + as.numeric(t(eta_est)%*%solve(Psi_0)%*%eta_est))
    llh = llh + log(2) - p*0.5*log(2*pi) - 0.5*log(det(Omega_0)) - temp_1 + log(pnorm(temp_2))
  }
  return(llh)
}

full_llh = function(xi_0,Psi_bar_0,omega_0,eta_0,h_0){
  llh = 0
  for(i in 1:N){
    y = vector(length = p)
    g = vector(length = p)
    temp_2 = 0
    for(j in 1:p){
      y[j] = (Y[i,j] - xi_0[j])/omega_0[j,j]
      g[j] = y[j]*exp(-0.5*W(h_0[j]*y[j]^2))
      temp_2 = temp_2 + 0.5*W(h_0[j]*y[j]^2) - log(h_0[j]*y[j]^2 + exp(W(h_0[j]*y[j]^2)))
    }
    temp_1 = -sum(log(diag(omega_0)))
    temp_3 = log(2) - 0.5*p*log(2*pi) - 0.5 * log(det(Psi_bar_0+eta_0%*%t(eta_0))) - 0.5 * as.numeric(t(g)%*%solve(Psi_bar_0+eta_0%*%t(eta_0))%*%g)
    con = as.numeric(t(eta_0)%*%solve(Psi_bar_0)%*%g)/sqrt(1+as.numeric(t(eta_0)%*%solve(Psi_bar_0)%*%eta_0))
    temp_4 = log(pnorm(con))
    llh = llh + temp_1+temp_2+temp_3+temp_4
  }
  return(llh)
}

ini_Psi = matrix(nrow = p,ncol = p)
for(i_0 in 1:p){
  ini_Psi[i_0,i_0] = 1
  for (j_0 in (i_0+1):p){
    if(j_0>p){
      break
    }else{
      ini_Psi[j_0,i_0] = cor(Y[,i_0],Y[,j_0])
      ini_Psi[i_0,j_0] = ini_Psi[j_0,i_0]
    }
  }
}

ini_Gamma = solve(ini_Psi)

alpha_sq = as.numeric(t(eta_est)%*%ini_Gamma%*%eta_est)
tau_bar = vector(length = N)
v_1 = vector(length = N)
v_2 = vector(length = N)
for(i in 1:N){
  tau_bar[i] = (1/sqrt(1+alpha_sq))*as.numeric(t(eta_est)%*%ini_Gamma%*%(Z_est[i,]))
  v_1[i] = (1/sqrt(1+alpha_sq))*(tau_bar[i]+(dnorm(tau_bar[i])/pnorm(tau_bar[i])))
  v_2[i] = (1/(1+alpha_sq))*(1 + tau_bar[i]^2 + tau_bar[i]*(dnorm(tau_bar[i])/pnorm(tau_bar[i])))
}

fin_Psi_term = matrix(data = rep(0,p^2),nrow = p)
for(i in 1:N){
  fin_Psi_term = fin_Psi_term + (v_1[i]/N)*(Z_est[i,])%*%t(eta_est) + (v_1[i]/N)*eta_est%*%t(Z_est[i,])
}

fin_Psi = matrix(data = rep(0,p^2),nrow = p)
for(i in 1:N){
  fin_Psi = fin_Psi+ (1/N)*(Z_est[i,])%*%t(Z_est[i,]) + (1/N)*v_2[i]*eta_est%*%t(eta_est)
}
fin_Psi = fin_Psi - fin_Psi_term
h_stop = abs((llh_SN(fin_Psi)/llh_SN(ini_Psi))-1)

while(h_stop>10^-9){
  ini_Psi = fin_Psi
  ini_Gamma = solve(ini_Psi)
  
  alpha_sq = as.numeric(t(eta_est)%*%ini_Gamma%*%eta_est)
  tau_bar = vector(length = N)
  v_1 = vector(length = N)
  v_2 = vector(length = N)
  for(i in 1:N){
    tau_bar[i] = (1/sqrt(1+alpha_sq))*as.numeric(t(eta_est)%*%ini_Gamma%*%(Z_est[i,]))
    v_1[i] = (1/sqrt(1+alpha_sq))*(tau_bar[i]+(dnorm(tau_bar[i])/pnorm(tau_bar[i])))
    v_2[i] = (1/(1+alpha_sq))*(1 + tau_bar[i]^2 + tau_bar[i]*(dnorm(tau_bar[i])/pnorm(tau_bar[i])))
  }
  
  fin_Psi_term = matrix(data = rep(0,p^2),nrow = p)
  for(i in 1:N){
    fin_Psi_term = fin_Psi_term + (v_1[i]/N)*(Z_est[i,])%*%t(eta_est) + (v_1[i]/N)*eta_est%*%t(Z_est[i,])
  }
  
  fin_Psi = matrix(data = rep(0,p^2),nrow = p)
  for(i in 1:N){
    fin_Psi = fin_Psi+ (1/N)*(Z_est[i,])%*%t(Z_est[i,]) + (1/N)*v_2[i]*eta_est%*%t(eta_est)
  }
  fin_Psi = fin_Psi - fin_Psi_term
  
  h_stop = abs((llh_SN(fin_Psi)/llh_SN(ini_Psi))-1)
}
diag_sd_inv = diag((sqrt(diag(fin_Psi)))^(-1))
fin_Psi_bar = diag_sd_inv %*% fin_Psi %*% diag_sd_inv

## EM is done, estimate of Psi_bar is in fin_Psi_bar

## MLE estimation ##

full_llh_2 = function(a){
  xi_0 = a[1:p]
  omega_0 = diag(exp(a[(p+1):(2*p)]))
  eta_0 = a[(2*p+1):(3*p)]
  h_0 = exp(a[(3*p+1):(4*p)])
  Psi_0 = matrix(nrow = p,ncol = p)
  k = 1
  for(i in 1:(p-1)){
    Psi_0[i,i] = 1
    for(j in (i+1):p){
      Psi_0[i,j] = a[4*p+k]/sqrt(1+a[4*p+k]^2)
      Psi_0[j,i] = Psi_0[i,j]
      k = k+1
    }
  }
  Psi_0[p,p] = 1
  
  if(min(eigen(Psi_0)$value)<=0){
    return(10^18)
  }else{
    return(-full_llh(xi_0,Psi_0,omega_0,eta_0,h_0))
  }
}


ini_Psi_mle = vector(length = 0.5*p*(p-1))
k = 1
for(i in 1:(p-1)){
  for(j in (i+1):p){
    ini_Psi_mle[k] = fin_Psi_bar[i,j]/sqrt(1-fin_Psi_bar[i,j]^2)
    k = k+1
  }
}
ini_MLE = c(xi_est,log(omega_est),eta_est,log(h_est),ini_Psi_mle)

H = optim(ini_MLE,full_llh_2,control=list(trace=1),method = "BFGS")

xi_est_mle = H$par[1:p]
omega_est_mle = exp(H$par[(p+1):(2*p)])
eta_est_mle = H$par[(2*p+1):(3*p)]
h_est_mle = exp(H$par[(3*p+1):(4*p)])
Psi_bar_est_mle = matrix(nrow = p,ncol = p)
k = 1
for(i in 1:(p-1)){
  Psi_bar_est_mle[i,i] = 1
  for(j in (i+1):p){
    Psi_bar_est_mle[i,j] = H$par[4*p+k]/sqrt(1+H$par[4*p+k]^2)
    Psi_bar_est_mle[j,i] = Psi_bar_est_mle[i,j]
    k = k+1
  }
}
Psi_bar_est_mle[p,p] = 1

## MLE is done, estimates are in xi_est_mle, omega_est_mle, eta_est_mle, h_est_mle, Psi_bar_est_mle
