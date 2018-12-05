rm(list=ls())
setwd('/Users/Owner/Documents/1.0Notre Dame/UNDERC')
#data=read.table('data2.txt',header=TRUE,sep=',')

library(deSolve)
library(ggplot2)
#define lotka-volterra competition model simulation
LV=function(t,y,p){
  n1=y[1]
  n2=y[2]
  r1=p[1]
  r2=p[2]
  a12=p[3]
  a11=p[4]
  a21=p[5]
  a22=p[6]
  
  dN1dt=r1*(1-n1*a11-n2*a12)*n1
  dN2dt=r2*(1-n2*a22-n1*a21)*n2
  
  return(list(c(dN1dt,dN2dt)))
}

#case1 a12>a11 and a21<a22, species 1 has a greater effect on species 2, species 1 outcompetes species 2
times=1:100
y0=c(0.5,0.5)
params=c(0.7,0.5,0.07,0.02,0.02,0.07)
sim=ode(y=y0,times=times,func=LV,parms=params)
out=data.frame(time=sim[,1],spec1=sim[,2],spec2=sim[,3])
ggplot(out,aes(x=time,y=spec1))+
  geom_line()+theme_classic()+
  geom_line(data=out,mapping=aes(x=time,y=spec2),color='green')+
  theme_classic()+ylab('species')

LVCR=function(t,y,p){
  H=y[1]#herbavor population
  P=y[2]#predator population
  
  b=p[1]#prey birth rate
  a=p[2]#predator attack rate
  e=p[3]#conversion efficeiency of prey to predtors
  s=p[4]#predator death rat3
  
  dHdt=b*H-a*P*H
  dPdt=e*a*P*H-s*P
  
  return(list(c(dHdt,dPdt)))
}
#sim 1
times=0:100#timestep 0.1
y0=c(25,5)
params=c(0.5,0.02,0.1,0.2)
sim=ode(y=y0,times=times,func=LVCR,parms=params)
out=data.frame(time=sim[,1],herb=sim[,2],pred=sim[,3])

ggplot(out,aes(x=time,y=herb))+
  geom_line()+theme_classic()+
  geom_line(data=out,mapping=aes(x=time,y=pred),color='green')+
  theme_classic()+ylab('species')

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
  theme_classic()+ylab('species')


####number 2######
RMmodel=function(t,y,p){
  #state variables
  H=y[1] #herbivore population
  P=y[2] #predator population 
  
  #parameters
  b=p[1] 
  e=p[2] 
  s=p[3] 
  w=p[4]
  d=p[5]
  alpha=p[6]
  
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

