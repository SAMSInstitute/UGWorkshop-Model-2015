require(truncnorm)
require(deSolve)

tpts = seq(0,500,length=4001)
x0 = c(1,1,1)

OBS_NOISE = F
obs_mu = c(0,0,0)
obs_sig2 = c(0.01,0.01,0.01)

#Dilution rate
D = 0.15


####################
#Note that this section does nothing if epsilon=0
#T = 100.0
TT = 24

eps_f <- function(STOC,eps){
  if(STOC){
    eps <- rtruncnorm(1, a=0, b=3, mean = eps, sd = 0.15)
  }else{
    eps
  }
}

omega_f <- function(STOC,T_arg){
  if(STOC){
    T_arg = rtruncnorm(1, a=0, b=150, mean = 50, sd = 10)
    2*pi/D/T_arg
  }else{
    2*pi/D/T_arg
  }
}

#create the random variable you want for T here
#remember to require T>0


####################
#Other parameters

#Si, mg/l
si = 115.0

#Maximum specific growth rate of prey and predator (1/hr)
mu1 = 0.5 #prey
mu2 = 0.2 #predator

#yield of prey per unit mass of substrate (dimensionless)
y1 = 0.4

#biomass yield of predator per unit mass of prey (dimensionsless)
y2 = 0.6

#half-saturation (Michaelis-Menten) constants for prey and predator
k1 = 8.0
k2 = 9.0

#Constants in dimensionless equations

##### ODE function #####
KotODEs <- function(t,x,params){
  with(as.list(c(x, params)),{
    
    omega = (2*pi/D/TT)*(1-STOC_T)+(2*pi/D/rtruncnorm(1, a=0, b=150, mean = 50, sd = 10))*(STOC_T)
    epsilon = epsilon*(1-STOC_EPS)+rtruncnorm(1, a=0, b=3, mean = epsilon, sd = 0.15)*(STOC_EPS)
    dx = 1 + epsilon*sin((omega*t)) - x[1] - A*x[1]*x[2]/(a+x[1]) # substrate
    dy = A*x[1]*x[2]/(a+x[1]) - x[2] - B*x[2]*x[3]/(b+x[2]) # prey (y) feeding on substrate (x) and eaten by predator (z)
    dz = B*x[2]*x[3]/(b+x[2]) - x[3] # predator(z) feeding on prey (y)
    list(c(dx,dy,dz))})
}


pars = c(epsilon = 0.6,
         STOC_EPS = F,
         TT = 24,
         STOC_T = F,
         A = mu1/D,
         a = k1/si,
         B = mu2/D,
         b = k2/y1/si)
out <- ode(y = x0, times = tpts, func = KotODEs, parms = pars,maxsteps = 10000)
matplot(out[ , 1], out[ , 2:4], type = "l", xlab = "time", ylab = "Conc",
        main = "substrate (x), prey (y) and predator (z) vs time", lwd = 2)
legend("topright", c("substrate","prey", "predator"), col = 1:3, lty = 1:3)

plot(out)  

