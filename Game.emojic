################ 
# NOTE: Change lines below to tol=.01 or tol=.001 and nboot=100 or nboot=200 
#            if this takes a long time to run - depends on your machine
            
tol = sqrt(.Machine$double.eps)  # determines convergence of optimizer     
nboot = 500                      # number of bootstrap replicates     
################# 

if (!'plyr' %in% installed.packages()[,1]) install.packages(plyr)   # check if plyr installed
library(plyr)
 
y = window(qinfl, c(1953,1), c(1965,2))  # inflation   
z = window(qintr, c(1953,1), c(1965,2))  # interest   
num = length(y) 
A = array(z, dim=c(1,1,num))
input = matrix(1,num,1) 
    
# Function to Calculate Likelihood   
Linn=function(para){
  phi=para[1]; alpha=para[2]; b=para[3]; Ups=(1-phi)*b
  cQ=para[4]; cR=para[5]  
  kf=Kfilter2(num,y,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)
  return(kf$like)    
  }

# Parameter Estimation   
mu0=1; Sigma0=.01; phi=.84; alpha=-.77; b=.85; cQ=.12; cR=1.1
init.par = c(phi,alpha,b,cQ,cR)    # initial parameters   
est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1,reltol=tol))  
SE = sqrt(diag(solve(est$hessian)))                     
phi = est$par[1]; alpha=est$par[2]; b=est$par[3]; Ups=(1-phi)*b         
cQ=est$par[4]; cR=est$par[5] 
cbind(estimate=est$par, SE)  

# BEGIN BOOTSTRAP   
# Likelihood for the bootstrapped data   
Linn2=function(para){
  phi=para[1]; alpha=para[2]; b=para[3]; Ups=(1-phi)*b
  cQ=para[4]; cR=para[5]  
  kf=Kfilter2(num,y.star,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)
  return(kf$like)  }

# Run the filter at the estimates   
kf = Kfilter2(num,y,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)          

# Pull out necessary values from the filter and initialize  
xp=kf$xp; innov=kf$innov; sig=kf$sig; K=kf$K; e=innov/sqrt(sig)
e.star=e; y.star=y; xp.star=xp; k=4:50 
para.star = matrix(0, nboot, 5)  # to store estimates
init.par=c(.84,-.77,.85,.12,1.1)    

pr <- progress_text()            # displays progress
pr$init(nboot)

for (i in 1:nboot){
 # cat("iteration:", i, "\n")  # old counter
 pr$step()                     # new counter
 e.star[k] = sample(e[k], replace=TRUE)   
 for (j in k){
  xp.star[j] = phi*xp.star[j-1]+Ups+K[j]*sqrt(sig[j])*e.star[j] 
  }   
 y.star[k] = z[k]*xp.star[k]+alpha+sqrt(sig[k])*e.star[k]  
 est.star = optim(init.par, Linn2, NULL, method="BFGS", control=list(reltol=tol))     
 para.star[i,] = cbind(est.star$par[1], est.star$par[2], est.star$par[3], abs(est.star$par[4]), 
                      abs(est.star$par[5]))
} 
                      
# Some summary statistics  
rmse = rep(NA,5)  # SEs from the bootstrap
for(i in 1:5){rmse[i]=sqrt(sum((para.star[,i]-est$par[i])^2)/nboot)
              cat(i, rmse[i],"\n") 
}
              
# phi and sigw  
phi = para.star[,1]; sigw = abs(para.star[,4]) 
phi = ifelse(phi<0, NA, phi)  # any φ < 0 not plotted

panel.hist <- function(x, ...){    # scatterplot with histograms
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

u = cbind(phi, sigw)
colnames(u) = c("f","s")
pairs(u, cex=1.5, pch=1, diag.panel=panel.hist, cex.labels=1.5, font.labels=5)
