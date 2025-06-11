load("data/fitting.rda")

B.psamp2=fit_mcmc_cpp2$B.psamp
A.psamp2=fit_mcmc_cpp2$A.psamp

B.psamp2 = B.psamp2[,,101:2000] 
A.psamp2 = A.psamp2[,,101:2000] 

protein_names = colnames(Y)
p = dim(A.psamp2)[1]
q = dim(B.psamp2)[2]
niter = dim(A.psamp2)[3]
n_time_points = nrow(unique(x2))


## Psi inverse
psi_inv_arr = array(NA, dim = c(p, p, niter))
for(i in 1:niter){
  psi_inv_arr[,,i] = solve(A.psamp2[,,i])
}


## Kappa_x
kappa_x = function(B, Psi, x){
  # x : unique(x2)
  kap_x = matrix(NA, nrow = n_time_points, ncol = niter)
  tmp = array(NA, dim = c(n_time_points, p, p, niter))
  for(i in 1:n_time_points){
    for(j in 1:niter){
      tmp[i,,,j] = B[,,j] %*% x[i,] %*% t(x[i,]) %*% t(B[,,j]) %*% psi_inv_arr[,,j]
    }
  }
  
  for (i in 1:n_time_points) {
    for (j in 1:niter) {
      kap_x[i, j] = sum(diag(tmp[i,,,j]))
    }
  }
  return(list(Kap_x=kap_x,K_tmp=tmp))
}

Kap = kappa_x(B.psamp2, A.psamp2, unique(x2))
kap_01 = Kap$Kap_x
kap_0 = kap_01[1:(n_time_points/2),]; dim(kap_0)

kap_1 = kap_01[(n_time_points/2 +1):n_time_points,]; dim(kap_1)
k_tmp = Kap$K_tmp


## Omega_x
omega_x = function(B, Psi, x){
  # x : unique(x2)
  omega_arr = array(NA, dim = c(n_time_points, p, p, niter))
  for(i in 1:n_time_points){
    for(j in 1:niter){
      omega_arr[i,,,j] = psi_inv_arr[,,j] - (1/(1+kap_01[i,j])) * psi_inv_arr[,,j] %*% k_tmp[i,,,j]
    }
  }
  omega_arr
}

omg_01 = omega_x(B.psamp2, A.psamp2, unique(x2))
omg_0 = omg_01[1:(n_time_points/2),,,]
omg_1 = omg_01[(n_time_points/2 +1):n_time_points,,,]


## R (scaled Omega)
R_x = function(B, Psi, x){
  # x : unique(x2)
  wi =  omg_01
  pr = array(NA, dim = c(n_time_points, p, p, niter))
  R = array(NA, dim = c(n_time_points, p, p, niter))
  for(i in 1:n_time_points){
    for(j in 1:niter){
      pr[i,,,j] = diag(1/sqrt(diag(wi[i,,,j])))
      R[i,,,j] = - pr[i,,,j] %*% wi[i,,,j] %*% pr[i,,,j]
      diag(R[i,,,j]) = 0
    }
  }
  R
}


R_01 = R_x(B.psamp2, A.psamp2, unique(x2))
R_0 = R_01[1:(n_time_points/2),,,]
R_1 = R_01[(n_time_points/2+1):n_time_points,,,]

R_dn = array(NA, dim = c((n_time_points/2), p, p, niter))
for(i in 1:(n_time_points/2)){
  for(j in 1:niter){
    R_dn[i,,,j] = R_1[i,,,j] - R_0[i,,,j]
  }
}

## credible interval
library(coda)
ci_u = array(NA, dim = c((n_time_points/2), p, p))
ci_l = array(NA, dim = c((n_time_points/2), p, p))
for (i in 1:(n_time_points/2)) {
  for (j in 1:p) {
    for (k in 1:p) {
      ci_l[i, j, k] = HPDinterval(as.mcmc(R_dn[i, j, k, ]), prob = 0.995)[,"lower"]
      ci_u[i, j, k] = HPDinterval(as.mcmc(R_dn[i, j, k, ]), prob = 0.995)[,"upper"]
    }
  }
}


dn = array(NA, dim = c((n_time_points/2), p, p))

for (i in 1:(n_time_points/2)) {
  for (j in 1:p) {
    for (k in 1:p) {
      if (ci_u[i,j,k] < 0 | ci_l[i,j,k] > 0){
        dn[i,j,k] = mean(R_dn[i,j,k,]) # posterior mean
      } else {
        dn[i,j,k] = 0
      }
    }
  }
}
dn