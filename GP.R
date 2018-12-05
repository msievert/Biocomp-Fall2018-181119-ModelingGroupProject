rm(list=ls())

library(deSolve)
library(ggplot2)

#define lotka-volterra competition model simulation
LVCR=function(t,y,p){
  H=y[1]# herbavor population
  P=y[2]# predator population
  
  b=p[1]# prey birth rate
  a=p[2]# predator attack rate
  e=p[3]# conversion efficeiency of prey to predtors
  s=p[4]# predator death rat3
  
  dHdt=b*H-a*P*H
  dPdt=e*a*P*H-s*P
  
  return(list(c(dHdt,dPdt)))
}
#sim 1
times=0:100#timestep 0.1
y0=c(25,5)
# params=c(b,a,e,s)
params=c(0.5,0.02,0.1,0.2)
sim=ode(y=y0,times=times,func=LVCR,parms=params)
out=data.frame(time=sim[,1],herb=sim[,2],pred=sim[,3])

ggplot(out,aes(x=time,y=herb))+
  geom_line()+theme_classic()+
  geom_line(data=out,mapping=aes(x=time,y=pred),color='green')+
  theme_classic()+ylab('Species Population')

#sim2 double predator population
times=0:100#timestep 0.1
y0=c(25,10)
params=c(0.5,0.02,0.1,0.2)
sim=ode(y=y0,times=times,func=LVCR,parms=params)
out=data.frame(time=sim[,1],herb=sim[,2],pred=sim[,3])

ggplot(out,aes(x=time,y=herb))+
  geom_line()+theme_classic()+
  geom_line(data=out,mapping=aes(x=time,y=pred),color='green')+
  theme_classic()+ylab('species')

#sim3 double prey population
times=0:100#timestep 0.1
y0=c(50,5)
params=c(0.5,0.02,0.1,0.2)
sim=ode(y=y0,times=times,func=LVCR,parms=params)
out=data.frame(time=sim[,1],herb=sim[,2],pred=sim[,3])

ggplot(out,aes(x=time,y=herb))+
  geom_line()+theme_classic()+
  geom_line(data=out,mapping=aes(x=time,y=pred),color='green')+
  theme_classic()+ylab('Species Population')


####number 2######
RMmodel=function(t,y,p){
  #state variables
  H=y[1] #herbivore population
  P=y[2] #predator population 
  
  #parameters
  b=p[1] # prey birth rate
  e=p[2] # conversion efficiency
  s=p[3] # predator death rate
  w=p[4] # saturating function
  d=p[5] # handling time
  alpha=p[6] # prey self limitation
  
  # calculate change in state variables with time, given parameter values and current value of state variables
  dHdt=(b*H*(1-alpha*H))-(w*(H/(d+H))*P)
  dPdt=(e*w*(H/(d+H))*P)-s*P
  # return list containing change in state variables with time
  return(list(c(dHdt, dPdt)))
}

# Scenario 1
### Define parameters, initial values for state variables, and time steps
# params=c(b, e, s, w, d, alpha)
params=c(0.8, 0.07, 0.2, 5, 400, 0.001)
# inits=c(H, P)
inits=c(500,120)
# times=sequence of time steps
times=seq(from=1, to=300, by = 1)
### Simulate the model using ode()
modelSim=ode(y=inits,times=times,func=RMmodel,parms=params)
modelDF=melt(data.frame(time=modelSim[,1], H=modelSim[,2], P=modelSim[,3]), id.vars = "time", variable.name = "population", value.name = "density")

a=ggplot(modelDF, aes(x=time, y=density, group=population))+
  geom_line(aes(color=population))+
  theme_classic()+ylab("population density")
a

