#--------program 1. generate basic portafolio

#------------- parameters 

#-------------parameters of Vasicek's model
#        alpha*(rbar-x1) + sigma.r dwt

sigma.r = 0.01384326 
alpha = 1.466184 
rbar = 0.06909307

#theta1-theta2*x1

eta.r = 0.1278
eta.s = 0.13
eta = c(eta.r,eta.s)

#------------------volatilities
sigma.s = 0.24
sigma.T = sigma.r*(1-exp(-alpha*T))/alpha 

#----------- portafolio's composition
ws = 0.6
wr = 0.3
# ws = seq(0.2,0.7,0.1)
w = c(wr,ws)

rho = 0.2

Sigmar = matrix(c(0, rho,0,sqrt(1-rho^2)),2,2,byrow=FALSE)

Sigmas = matrix(c(sigma.T, 0,0,sigma.s),2,2,byrow=FALSE)

Sigma = Sigmas%*%Sigmar

(sigma.onda2 = (t(w)%*%Sigma)%*%(t(Sigma)%*%w))
(theta = t(w)%*%Sigmas%*%eta)

#---------------------------1. portfolio's simulacion 
require(yuima)

param <- list(rbar = rbar, 
alpha = alpha, 
theta = theta,
sigma.r = sigma.r,
sigma.T = sigma.T,
sigma.s = sigma.s,
ws = ws,
wr = wr,
rho = rho)

#-----state system for (rt,At) 

sol <- c("x1","x2") # variable for numerical solution

Df <- c("alpha*(rbar-x1)","x2*(x1+theta)") # drift vector

Sf <- matrix(c("sigma.r",
"x2*(rho*sigma.s*ws+sigma.T*wr)","0",
"x2*ws*sigma.s*sqrt(1-rho^2)"),2,2) # diffusion matrix


mod3 <- setModel(drift = Df, 
diffusion = Sf,state.variable=c("x1","x2"), 
solve.variable=c("x1","x2") )
mod3

x = 62; omega = 110; 
m = 12
T = omega-x
n =T*m
tx = seq(0,T-1/m,1/m)
length(tx)

sampling <- setSampling(Initial = 0.0, 
Terminal = T, n =(T-1/m)*m, delta=1/m)

cmod3 <-setYuima(model=mod3, sampling=sampling)
 

#------------calculate At Xt

#--------vectorizar
# muxt = Vectorize(muxt,"t")
# tpx = Vectorize(tpx,"t")
#---se coloca T-1/m para evitar el valor 0 de acx en T
# T = 110-65
# avxm.v = Vectorize(avxm,"x")
# factor.rp = function(tx,x){exp(-delta*tx)*tpx(tx,x,pars)*avxm.v(x+tx,i,m,pars)/avxm(x,i,m,pars)}

