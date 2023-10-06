library("mvtnorm")
library("LambertW")
library("sn")

data(wines)
x_wine = which(wines$wine=="Grignolino")

Y = wines[x_wine,c(13,15,24)]

N = nrow(Y)
p = ncol(Y)


xi_est = vector(length = p)
omega_est = vector(length = p)
eta_est = vector(length = p)
h_est = vector(length = p)
tau_h_Z_est = matrix(nrow = N,ncol = p)
Z_est = matrix(nrow = N,ncol = p)

snth_ini_time = Sys.time()
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
    
    z = (Y[,k] - xi)/omega
    llh = -log(omega) + 0.5*W(h*z^2) - log(h*z^2 + exp(W(h*z^2))) + 0.5*log(2) - 0.5*log(pi) - 0.5*log(1+eta^2)- 0.5*g(z)^2/(1+eta^2) + pnorm(eta*g(z)/sqrt(1+eta^2),log.p = TRUE)
    return(sum(llh))
  }
  negative_llh = function(theta){
    return(-marginal_llh(theta))
  }
  H_optim = optim(c(mean(Y[,k]),log(sd(Y[,k])),skewness(Y[,k])/sqrt(1+skewness(Y[,k])^2),log(kurtosis(Y[,k])/sqrt(1+kurtosis(Y[,k])^2))),negative_llh,method = "BFGS")
  xi_est[k] = H_optim$par[1]
  omega_est[k] = exp(H_optim$par[2])
  eta_est[k] = H_optim$par[3]
  h_est[k] = exp(H_optim$par[4])
  
  tau_h_Z_est[,k] = (Y[,k] - xi_est[k])/omega_est[k]
  
  Z_est[,k] = W_delta(tau_h_Z_est[,k], h_est[k])
}

## marginal MLE done ##
## estimates are: xi_est, omega_est, eta_est, h_est

## EM for SN ##
sum_of_square = function(x){
  return(sum(x^2))
}

llh_SN = function(Psi_0){
  Omega_0 = Psi_0 + eta_est%*%t(eta_est)
  Inv_Omega_0 = solve(Omega_0)
  chol_Inv_Omega_0 = t(chol(Inv_Omega_0))
  
  Inv_Psi_0 = solve(Psi_0)
  
  m_1 = Z_est%*%t(chol(Inv_Omega_0))
  temp_1 = 0.5*apply(m_1,1,FUN = sum_of_square)
  temp_2 = (Z_est)%*%Inv_Psi_0%*%eta_est/sqrt(1 + as.numeric(t(eta_est)%*%Inv_Psi_0%*%eta_est))
  
  llh = log(2) - p*0.5*log(2*pi) - 0.5*log(det(Omega_0)) - temp_1 + pnorm(temp_2,log.p = TRUE)
  return(sum(llh))
}

ini_Psi = cor(Y)
ini_Gamma = solve(ini_Psi)

dn_by_pn = function(x){
  if(pnorm(x)==0){
    return(0)
  }else{
    return(dnorm(x)/pnorm(x))
  }
}
dn_by_pn = Vectorize(dn_by_pn)

f_1 = function(A){
  return(A%*%t(eta_est))
}

alpha_sq = as.numeric(t(eta_est)%*%ini_Gamma%*%eta_est)
tau_bar = (1/sqrt(1+alpha_sq))*as.numeric(Z_est%*%ini_Gamma%*%eta_est)
v_1 = (1/sqrt(1+alpha_sq))*(tau_bar+(dn_by_pn(tau_bar)))
v_2 = (1/(1+alpha_sq))*(1 + tau_bar^2 + tau_bar*(dn_by_pn(tau_bar)))

temp = apply(apply(v_1*Z_est,1, FUN = f_1),1,sum)
fin_Psi = (1/N)*(t(Z_est)%*%(Z_est) + sum(v_2)*eta_est%*%t(eta_est)) - (matrix(temp,byrow = TRUE,nrow = p)+matrix(temp,byrow = FALSE,nrow = p))/N

h_stop = abs((llh_SN(fin_Psi)/llh_SN(ini_Psi))-1)

while(h_stop>10^-9){
  ini_Psi = fin_Psi
  ini_Gamma = solve(ini_Psi)
  
  alpha_sq = as.numeric(t(eta_est)%*%ini_Gamma%*%eta_est)
  tau_bar = (1/sqrt(1+alpha_sq))*as.numeric(Z_est%*%ini_Gamma%*%eta_est)
  v_1 = (1/sqrt(1+alpha_sq))*(tau_bar+(dn_by_pn(tau_bar)))
  v_2 = (1/(1+alpha_sq))*(1 + tau_bar^2 + tau_bar*(dn_by_pn(tau_bar)))
  
  temp = apply(apply(v_1*Z_est,1, FUN = f_1),1,sum)
  fin_Psi = (1/N)*(t(Z_est)%*%(Z_est) + sum(v_2)*eta_est%*%t(eta_est)) - (matrix(temp,byrow = TRUE,nrow = p)+matrix(temp,byrow = FALSE,nrow = p))/N
  
  h_stop = abs((llh_SN(fin_Psi)/llh_SN(ini_Psi))-1)
  print(h_stop)
}
diag_sd_inv = diag((sqrt(diag(fin_Psi)))^(-1))
fin_Psi_bar = diag_sd_inv %*% fin_Psi %*% diag_sd_inv

## EM is done, estimate of Psi_bar is in fin_Psi_bar

## MLE estimation ##

full_llh = function(xi_0,Psi_bar_0,omega_0,eta_0,h_0){
  y = t(solve(omega_0)%*%(t(Y) - xi_0))
  g = matrix(nrow = N,ncol = p)
  temp_2 = matrix(nrow = N,ncol = p)
  for(j in 1:p){
    g[,j] = y[,j]*exp(-0.5*W(h_0[j]*y[,j]^2))
    temp_2[,j] =  0.5*W(h_0[j]*y[,j]^2) - log(h_0[j]*y[,j]^2 + exp(W(h_0[j]*y[,j]^2)))
  }
  
  Inv_Psi_0 = solve(Psi_bar_0)
  
  Inv_Omega_0 = solve(Psi_bar_0+eta_0%*%t(eta_0))
  m_1 = g%*%t(chol(Inv_Omega_0))
  temp = 0.5*apply(m_1,1,FUN = sum_of_square)
  temp_1 = -N*sum(log(diag(omega_0)))
  temp_3 = log(2) - 0.5*p*log(2*pi) - 0.5 * log(det(Psi_bar_0+eta_0%*%t(eta_0))) - temp
  temp_4 = pnorm((g%*%Inv_Psi_0%*%eta_0)/sqrt(1+as.numeric(t(eta_0)%*%Inv_Psi_0%*%eta_0)),log.p = TRUE)
  llh = temp_1 + sum(temp_2) + sum(temp_3) + sum(temp_4)
  return(llh)
}

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
snth_fin_time = Sys.time()
x = matrix(rep(1,N),ncol = 1)
st_ini_time = Sys.time()
skew_t = selm.fit(x,as.matrix(Y),family="ST")
st_fin_time = Sys.time()
aic_skew_t = 2*(2*p + 1 + (0.5*p*(p+1))) - 2*skew_t$logL
aic_snth = 2*(4*p + (0.5*p*(p-1))) + 2*H$value

snth_time = snth_fin_time - snth_ini_time
st_time = st_fin_time - st_ini_time

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
