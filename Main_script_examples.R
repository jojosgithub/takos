
#REQUIRED PACKAGES
require(data.table)
require(MASS)

library(takos)

###########################################################################################################
#Data Simulation
###########################################################################################################

#user can access the generic gAC function for simulating thermograms setting up the  kinetic parameters and chosen models
#supported models are RO1,RO2,RO3,SB,P1,P2,P3,p4,D1,F1,y14,y13,y12,D3,D4,R2,R3,D2,JMA,Ih,F2

#gAC: simulate a thermogram using a kinetic model chosen by the user
RO1_mod <- gAC(A = exp(35),Ea = 120000,q = 50,T0 = -100,T.end = 300,n=2,npoints=100,prec=10^(-5),rmod="RO1")
RO2_mod <- gAC(A = exp(35),Ea = 120000,q = 50,T0 = -100,T.end = 300,n=2,npoints=100,prec=10^(-5),rmod="RO2")
F2_mod <- gAC(A = exp(35),Ea = 120000,q = 50,T0 = -100,T.end = 300,n=2,npoints=100,prec=10^(-5),rmod="F2")


#Due to their common use JMA and SB have separetaded functions that can be called using JMA and  sbAC 

JMA_mod <- JMA(A = exp(35),Ea = 120000,q = 50,T0 = -100,T.end = 300,npoints=100,n=2)
SB_mod <- sbAC(A = exp(35),Ea = 120000,q = 50,T0 = -100,T.end = 300,npoints=100,n=2)


###########################################################################################################
#Example of simulation of a series of thermogram at different heating rate and calculating the kinetic 
#triplet using several methodologies
par(mfrow=c(1,1))
#choose the rates for the simulation of the thermograms
rates=c(0.5,1,2,5,10,20,50)
#first serie of thermograms for all the chosen rate
a<-lapply(rates, function(x) JMA(A=exp(35),Ea=120000,T0=0,T.end=300,q=x,npoints=500,n=2))
#setup column names
a<-lapply(seq(1,length(a)), function(x) data.table(a[[x]]$time.s, a[[x]]$T.C, a[[x]]$dadT, rates[[x]]))
lapply(seq(1,length(a)), function(x) setnames(a[[x]], c("time.seconds","temperature.s","heat.flow","rates") ) )
###########################################################################################################
#create a plot using the function thermo
amaxH <- max(sapply(a, function(x) max(x$heat.flow))) # calculate the max
plot(c(0,300),c(0,amaxH),main="dataset A 120/60 0.66/0.33",
     ylab="ExothermicHeatFlow", xlab="Temperature")
lapply(a, function(x) lines(x$temperature.s,x$heat.flow,lwd=3))
###########################################################################################################

#prepare the matrix for the calculations using tmat
#keep in mind that we are recalculating alfa even if for simulated data this is not necessary
#but we could just simply copy the columns as ri

tmat <- testMat(a,toselect=c(0,1,2,0,0,0,3,4))

#FRI: Friedman method to calculate the activation energy (Ea)
fri<-FRI(tmat)
summaryTableFri(fri)
plot_fri(fri)

# avrami: Avrami method
avr<-avrami(tmat)
summaryTableA(avr)
plot_avrami(avr)

# lavrami:  linearization of avrami method
lavr<-avrami(tmat)
summaryTableA(lavr)


#KAS: Kissinger-Akahira-Sunose (KAS) method
kas<-KAS(tmat)

#Kiss: Kissinger method to calculate the activation energy (Ea)
kiss<-Kiss(tmat)

#MO: Mo method

mo<-MO(tmat)
summaryTableMo(mo)
plot_mo(mo)

#OFW: Ozawa-Flynn and Wall method
ofw<-OFW(tmat)

#OZ: Ozawa method
oz<-OZ(tmat)
summaryTableOz(oz)
plot_ozawa(oz)

#Starink: Starink method
star<-Starink(tmat)

#VY: Vyazovkin isoconversional method to calculate the activation energy (Ea
#we first select the degree of alfa
as<-select_degree(tmat)
vy<-as[, optimize(function(x) VY(temperature.s.K,rate,x), lower=50,upper=250),by=rit]
vy


###########################################################################################################
#Other functions
###########################################################################################################

#simG: create a simulated spectra with gaussian shape
y=(simG(500,35,1,0,w=20))
plot(y)

#cutSelect and cutValue: cut a region of a spectra and substitutes it with a sequence with initial value i.start and end valye i.end
npoints=1000
x=seq(1,npoints)
y=(dnorm(x, mean=npoints/2, sd=npoints/10))
ycut=cutSelect(y,10,40)
plot(y)
lines(ycut,col="red")

ycut=cutValue(y,10,40,0.003,0.001)
plot(y)
lines(ycut,col="red")

#dadx: calculates the ratio of two differential according to the value of d.step (default =2)
npoints=100
seed=42
x1=round(runif(npoints,0,1), 2)
seed=1234
x2=round(runif(npoints,0,1), 2)
xdiff <- dadx(x1,x2)

#runningIntegral: calculates the running integral for customer input
npoints=1000
x=seq(1,npoints)
y=(dnorm(x, mean=npoints/2, sd=npoints/10))
runningIntegral(x,y)

#RI: calculate the running integral for the selected peak

#data simulation
rates=c(0.5,1,2,5,10,20,50)
a<-lapply(rates, function(x) JMA(A=exp(35),Ea=120000,T0=0,T.end=300,q=x,npoints=1000,n=2))
a<-lapply(seq(1,length(a)), function(x) data.table(a[[x]]$time.s,a[[x]]$T.C,
                                                   a[[x]]$dadT, rates[[x]]))
lapply(seq(1,length(a)), function(x) setnames(a[[x]],
                                              c("time.seconds","temperature.s","heat.flow","rates") ) )
a.dt <-lapply(seq(1,length(a)), function(x) data.table(data.frame(a[[x]])))
a<-rbindlist(a.dt)
a$id<-a$rates
require(pracma)
a.peaks <- a[,.(res.list = list(findpeaks(heat.flow,sortstr=TRUE,npeaks=2))),by=id]
a.peaks$rate<-a.peaks$id
ref.peak=1
a.peaks <- data.table(data.table(a.peaks$rate),rbindlist((lapply(a.peaks$res.list,
                                                                 function(x) data.table(t(x[ref.peak,]))))))
colnames(a.peaks)<- c("rate","peak.value","ind.max","left.lim","right.lim")
a.mat<- lapply(unique(a$rate),function(x)
  #ri
  ri(a[a$rate==x]$time.seconds,a[a$rate==x]$heat.flow,a.peaks[rate==x]))

#smooth.loess: a wrapper for the loess function included in the R base system
npoints=500
x=seq(1,npoints)
y=(dnorm(x, mean=npoints/2, sd=npoints/10))
y.smooth=smooth.loess(x,y)
plot(x,y)

# checkmat: returns checked data frame
npoints=500
x=seq(1,npoints)
y=(dnorm(x, mean=npoints/2, sd=npoints/10))
x=seq(1,1000)
x2=seq(200,500,length.out=1000)
dat=data.frame(x,x2,y)
colnames(dat) <- c("time.seconds", "temperature.s","heat.flow")
cmat<- checkmat(dat,selected=c(1,0,2,0,0,0,3,0))

#testMat return data table ready to be used by all the methods for kinetic analysis
dat=data.table(dat)
dat2=dat
dat$rates=20
dat2$rates=50
toTest=list(dat,dat2)
tested=testMat(toTest)

#TAPPA: calculates the background of a thermogram according to Tangent-area-proportional method
npoints=1000
x=seq(1,npoints)
y=(dnorm(seq(1,npoints), mean=npoints/2, sd=npoints/10)) #simulated peak
y2=y+(dnorm(seq(1,npoints), mean=npoints, sd=npoints/10)) #secondary simulated peak
y2[seq(npoints*0.735,npoints)]=y2[763] #flat the curve at the end of first peak
ytap=TAPPA(x,y2)
#example plot
plot(x,y2)
lines(x,ytap,col="red")


