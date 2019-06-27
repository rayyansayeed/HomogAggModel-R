setwd("C:/Users/Rayyan Sayeed/Documents/Northwestern Spring 2018/MKTG 551-3/Homework 3")
library(dplyr, tidyr)
library(MASS)
df=read.csv("car_panel.csv",header=T)

# create market share vector, year sum (ownshare summed by year)

df$ownshare=df$q/df$num_hh
yrsum=aggregate(ownshare~year,data=df,FUN=sum)

# create vector to match dimension of yrsum with datafile

dim1=as.vector(t(table(df$year)))

# create outside share vector

df$outsideshare=1-as.numeric(rep(yrsum$ownshare,times=dim1))

# mean utility

df$meanutility=log(df$ownshare)-log(df$outsideshare)

############# create BLP-instruments now...

###### COST-IV vars already in df, steel_index and wages_index

###### OWN-IV

# create vector to match dimensions
dim2=as.vector(t(table(df$year,df$firm_id)))
dim2 = dim2[dim2!=0]

# Air


own_ag_ac=aggregate(air~firm_id+year,data=df,FUN=sum)
own_ac_vec=rep(own_ag_ac$air,times=dim2)
own_ac_vec=own_ac_vec/(rep(dim2,times = dim2)-1)

own_ac_vec[is.infinite(own_ac_vec)]=0
own_ac_vec[is.nan(own_ac_vec)]=0

# Miles/dollar ratio

own_ag_mpd=aggregate(mpd~firm_id+year,data=df,FUN=sum)
own_mpd_vec=rep(own_ag_mpd$mpd,times=dim2)
own_mpd_vec=own_mpd_vec/(rep(dim2,times = dim2)-1)

own_mpd_vec[is.infinite(own_mpd_vec)]=0
own_mpd_vec[is.nan(own_mpd_vec)]=0


# Size

own_ag_size=aggregate(log_car_size~firm_id+year,data=df,FUN=sum)
own_size_vec=rep(own_ag_size$log_car_size,times=dim2)
own_size_vec=own_size_vec/(rep(dim2,times = dim2)-1)

own_size_vec[is.infinite(own_size_vec)]=0
own_size_vec[is.nan(own_size_vec)]=0

# HP to Weight Ratio

own_ag_hp2wt=aggregate(hp2wt~firm_id+year,data=df,FUN=sum)
own_hp2wt_vec=rep(own_ag_hp2wt$hp2wt,times=dim2)
own_hp2wt_vec=own_hp2wt_vec/(rep(dim2,times = dim2)-1)

own_hp2wt_vec[is.infinite(own_hp2wt_vec)]=0
own_hp2wt_vec[is.nan(own_hp2wt_vec)]=0

########## OTHER-IV

dim3=as.vector(t(table(df$year)))
dim3=dim3[dim3 != 0]

# Air

other_ag_ac=aggregate(air~year,data=df,FUN=sum)
other_ac_vec=rep(other_ag_ac$air,times=dim3)
other_ac_vec=other_ac_vec/(rep(dim3,times = dim3)-1)

other_ac_vec[is.infinite(other_ac_vec)]=0
other_ac_vec[is.nan(other_ac_vec)]=0

# Miles/dollar ratio

other_ag_mpd=aggregate(mpd~year,data=df,FUN=sum)
other_mpd_vec=rep(other_ag_mpd$mpd,times=dim3)
other_mpd_vec=other_mpd_vec/(rep(dim3,times = dim3)-1)

other_mpd_vec[is.infinite(other_mpd_vec)]=0
other_mpd_vec[is.nan(other_mpd_vec)]=0

# Size

other_ag_size=aggregate(log_car_size~year,data=df,FUN=sum)
other_size_vec=rep(other_ag_size$log_car_size,times=dim3)
other_size_vec=other_size_vec/(rep(dim3,times = dim3)-1)

other_size_vec[is.infinite(other_size_vec)]=0
other_size_vec[is.nan(other_size_vec)]=0

# HP to Weight Ratio

other_ag_hp2wt=aggregate(hp2wt~year,data=df,FUN=sum)
other_hp2wt_vec=rep(other_ag_hp2wt$hp2wt,times=dim3)
other_hp2wt_vec=other_hp2wt_vec/(rep(dim3,times = dim3)-1)

other_hp2wt_vec[is.infinite(other_hp2wt_vec)]=0
other_hp2wt_vec[is.nan(other_hp2wt_vec)]=0


############ START OF INSTRUMENTS, 2SLS



# Implement homogenous aggregate logit model with utility

hal_model=lm(meanutility~hp2wt+air+mpd+car_size+p_real, data=df)

# Standard Errors
coef(summary(hal_model))[, "Std. Error"] 

# Estimate model using 2SLS and a) COST-IV b) OWN-IV c) OTHER-IV d) All of the above

markup_vec=df$p_real+1/(1-df$ownshare)

# Estimate Price, p_hat, to use in 2nd stage utility estimation

p_hat=lm(p_real~steel_index+wages_index+1/(1-ownshare), data=df)

# Estimate utility using p_hat

utility=lm(meanutility~hp2wt+air+mpd+car_size+p_hat$fitted.values, data=df)

########################## 2SLS implementation with IVs

# Set up X, input matrix

ones = cbind(rep(1,length(df$id)))
X = cbind(ones,df$hp2wt,df$air,df$mpd,df$car_size,df$p_real)
M = nrow(df) # number of product-market observations, used in variance parts below

################# Part a) Cost-IV

# Z is IV matrix

Z_cost=cbind(df$steel_index,df$wages_index)

# Compute projection matrix of Z_cost, P_Z=Z(Z'Z)^-1*Z'

P_Z_cost=Z_cost%*%ginv((t(Z_cost)%*%Z_cost))%*%t(Z_cost)

# Compute fitted values of X using X_hat=P_Z*X

X_hat_Z_cost=P_Z_cost%*%X

# 2nd stage estimates, beta_hat_2sls=(X_hat'*X)^-1*X_hat'*delta, delta is mean utility

sls_estimates_Z_cost=ginv(t(X_hat_Z_cost)%*%(X_hat_Z_cost))%*%t(X_hat_Z_cost)%*%df$meanutility

# Residuals, epsilon_hat=delta-X*beta_hat_2sls

residuals_Z_cost= df$meanutility-X%*%(sls_estimates_Z_cost)

# Variance, V_2sls=1/M(epsilon_hat')(epsilon_hat)(X'P_zX)^-1, M = # of product-market obsdim

resid_term_cost=(t(residuals_Z_cost)%*%residuals_Z_cost)

V_2sls_Z_cost=(1/M)*resid_term_cost[1,1]*ginv(t(X)%*%P_Z_cost%*%X)

# Standard Error

std_err_cost=V_2sls_Z_cost/sqrt(M)


################# Part b) Own-IV

Z_own=cbind(own_ac_vec,own_mpd_vec,own_size_vec,own_hp2wt_vec)

P_Z_own=Z_own%*%ginv((t(Z_own)%*%Z_own))%*%t(Z_own)

X_hat_Z_own=P_Z_own%*%X

sls_estimates_Z_own=ginv(t(X_hat_Z_own)%*%(X_hat_Z_own))%*%t(X_hat_Z_own)%*%df$meanutility

residuals_Z_own= df$meanutility-X%*%sls_estimates_Z_own

resid_term_own=(t(residuals_Z_own)%*%residuals_Z_own)

V_2sls_Z_own=(1/M)*resid_term_own[1,1]*ginv(t(X)%*%P_Z_own%*%X)

std_err_own=V_2sls_Z_own/sqrt(M)


################# Part c) Other-IV

Z_other=cbind(other_ac_vec,other_mpd_vec,other_size_vec,other_hp2wt_vec)

P_Z_other=Z_other%*%ginv((t(Z_other)%*%Z_other))%*%t(Z_other)

X_hat_Z_other=P_Z_other%*%X

sls_estimates_Z_other=ginv(t(X_hat_Z_other)%*%(X_hat_Z_other))%*%t(X_hat_Z_other)%*%df$meanutility

residuals_Z_other= df$meanutility-X%*%sls_estimates_Z_other

resid_term_other=(t(residuals_Z_other)%*%residuals_Z_other)

V_2sls_Z_other=(1/M)*resid_term_other[1,1]*ginv(t(X)%*%P_Z_other%*%X)

std_err_other=V_2sls_Z_other/sqrt(M)


################# Part d) All of the above

Z_all=cbind(df$steel_index,df$wages_index,own_ac_vec,own_mpd_vec,own_size_vec,own_hp2wt_vec,other_ac_vec,other_mpd_vec,other_size_vec,other_hp2wt_vec)

P_Z_all=Z_all%*%ginv((t(Z_all)%*%Z_all))%*%t(Z_all)

X_hat_Z_all=P_Z_all%*%X

sls_estimates_Z_all=ginv(t(X_hat_Z_all)%*%(X_hat_Z_all))%*%t(X_hat_Z_all)%*%df$meanutility

residuals_Z_all= df$meanutility-X%*%sls_estimates_Z_all

resid_term_all=(t(residuals_Z_all)%*%residuals_Z_all)

V_2sls_Z_all=(1/M)*resid_term_all[1,1]*ginv(t(X)%*%P_Z_all%*%X)

std_err_all=V_2sls_Z_all/sqrt(M)

############ END OF INSTRUMENTS, 2SLS


############ GMM estimation

gmm = function(theta,data,y,Z){
  # GMM estimation function
  # Inputs:
  # data: here, X matrix
  # theta: starting pt of vector of parameters, use 2SLS estimates here
  # y: vector of outcomes (here, mean utility)
  # Z: matrix of instruments
  # Output: GMM val
  
  residuals = y - as.matrix(data)%*%as.matrix(theta)
  
  # Projection matrix: Z(Z'Z)^{-1}Z'
  
  proj_matrix = Z%*%ginv(t(Z)%*%Z)%*%t(Z)
  
  weight = ginv(t(Z)%*%Z)
  
  gmm_val = t(residuals)%*%Z%*%weight%*%t(Z)%*%residuals
  
  return(gmm_val)
}

######## a) GMM estimation, Cost-IV

mean_utility_col = cbind(df$meanutility)

gmm_cost = nlm(gmm,sls_estimates_Z_cost,X,mean_utility_col,Z_cost,print.level = 2, hessian=TRUE)
hansen_cost = gmm_cost$minimum

# Standard Error

H_cost=gmm_cost$hessian
V_cost=solve(H_cost)
stderr_cost=sqrt(diag(V_cost))

# V_cost has negative values on diagonal... how to take sqrt for std err?


######## b) GMM estimation, Own-IV

gmm_own = nlm(gmm,sls_estimates_Z_own,X,mean_utility_col,Z_own,print.level = 2, hessian=TRUE)
hansen_own = gmm_own$minimum


# Standard Error

H_own=gmm_own$hessian
V_own=solve(H_own)
stderr_own=sqrt(diag(V_own))


######## c) GMM estimation, Other-IV

gmm_other = nlm(gmm,sls_estimates_Z_other,X,mean_utility_col,Z_other,print.level = 2, hessian=TRUE)
hansen_other = gmm_other$minimum


# Standard Error

H_other=gmm_other$hessian
V_other=solve(H_other)
stderr_other=sqrt(diag(V_other))


######## d) GMM estimation, All-IV

gmm_all = nlm(gmm,sls_estimates_Z_all,X,mean_utility_col,Z_all,print.level = 2, hessian=TRUE)
hansen_all = gmm_all$minimum


# Standard Error

H_all=gmm_all$hessian
V_all=solve(H_all)
stderr_all=sqrt(diag(V_all))


############ Marginal cost, markup for Cost-IV

alpha = gmm_cost$estimate[6]
own_price = alpha*mean(df$p_real%*%(1-df$meanutility))
part=alpha*mean(df$meanutility%*%(1-df$meanutility))
marg_cost=mean(df$p_real+alpha*(1-df$meanutility))
markup=mean((1-df$meanutility)- alpha)

######### END OF CODE