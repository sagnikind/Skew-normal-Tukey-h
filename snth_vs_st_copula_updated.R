library("VineCopula")
library("mvtnorm")
library("LambertW")
library("sn")
N = 500
p = 3
# define 3-dimensional R-vine tree structure matrix
Matrix <- c(2, 3, 1,
            0, 3, 1,
            0, 0, 1)
Matrix <- matrix(Matrix, 3, 3)
# define R-vine pair-copula family matrix
family <- c(0, 4, 1,
            0, 0, 3,
            0, 0, 0)
family <- matrix(family, 3, 3)
# define R-vine pair-copula parameter matrix
par <- c(0, 1.9, 0.5,
         0, 0, 4.8,
         0, 0, 0)
par <- matrix(par, 3, 3)
# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 3, 3)
## define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3"))
## simulate a sample of size 500 from the R-vine copula model
set.seed(123)
simdata <- RVineSim(500, RVM)
# ## determine the pair-copula families and parameters
# RVM1 <- RVineCopSelect(simdata, familyset = c(1, 3, 4), Matrix)
# ## see the object's content or a summary
# str(RVM1)
# summary(RVM1)
# ## inspect the fitted model using plots
# ## Not run: plot(RVM1) # tree structure
# contour(RVM1) # contour plots of all pair-copulas


Y = matrix(0,nrow = N,ncol = p)
mu = c(0,0,0)
sig = c(1,1,1)

for(i in 1:N){
  Y[i,1] = mu[1] + sig[1]*qnorm(simdata[i,1])
  Y[i,2] = mu[2] + sig[2]*qt(simdata[i,2],df = 10)
  Y[i,3] = mu[3] + sig[3]*qt(simdata[i,3],df = 1)
}

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
  # print(H_optim$value)
  # print(negative_llh(c(xi[k],log(omega[k,k]),eta[k],log(h[k]))))
  
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
  # print(k)
}

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

full_llh(xi_est,fin_Psi_bar,diag(omega_est),eta_est,h_est)


full_llh_2 = function(a){
  xi_0 = a[1:3]
  omega_0 = diag(exp(a[4:6]))
  eta_0 = a[7:9]
  h_0 = exp(a[10:12])
  Psi_0 = matrix(nrow = p,ncol = p)
  for(i in 1:p){
    Psi_0[i,i] = 1
  }
  Psi_0[1,2] = a[13]/sqrt(1+a[13]^2)
  Psi_0[1,3] = a[14]/sqrt(1+a[14]^2)
  Psi_0[2,3] = a[15]/sqrt(1+a[15]^2)
  Psi_0[2,1] = Psi_0[1,2]
  Psi_0[3,1] = Psi_0[1,3]
  Psi_0[3,2] = Psi_0[2,3]
  
  if(min(eigen(Psi_0)$value)<=0){
    return(10^18)
  }else{
    return(-full_llh(xi_0,Psi_0,omega_0,eta_0,h_0))
  }
}


ini_MLE = c(xi_est,log(omega_est),eta_est,log(h_est),c(fin_Psi_bar[1,2]/sqrt(1-fin_Psi_bar[1,2]^2),fin_Psi_bar[1,3]/sqrt(1-fin_Psi_bar[1,3]^2),fin_Psi_bar[2,3]/sqrt(1-fin_Psi_bar[2,3]^2)))


H = optim(ini_MLE,full_llh_2,control=list(trace=1),method = "BFGS")
x = matrix(rep(1,N),ncol = 1)
skew_t = selm.fit(x,as.matrix(Y),family="ST")
aic_skew_t = 2*13 - 2*skew_t$logL
aic_snth = 2*15 + 2*H$value

xi_est_mle = H$par[1:3]
omega_est_mle = exp(H$par[4:6])
eta_est_mle = H$par[7:9]
h_est_mle = exp(H$par[10:12])
Psi_bar_est_mle = matrix(c(1,H$par[13]/sqrt(1+H$par[13]^2),H$par[14]/sqrt(1+H$par[14]^2),
                           H$par[13]/sqrt(1+H$par[13]^2),1,H$par[15]/sqrt(1+H$par[15]^2),
                           H$par[14]/sqrt(1+H$par[14]^2),H$par[15]/sqrt(1+H$par[15]^2),1),nrow=3)
