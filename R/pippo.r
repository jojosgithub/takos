#aggiungo descrizione

sbAC <-  function (time.start=0,T0=0,T.end=1200,qqq=10,A=10^(6.3),Ea=80000,m=1,n=2,npoints=10000,prec=10^(-4.30095790876)) 
{
library(deSolve)
R=8.314  #gas constant
Ts=273.15+ T0 #transformation in K
time.e=(T.end-T0)/(qqq/60) #estimated time for the analysis
time.s=seq(time.start,time.e, length.out=npoints) #vector with all the times
Temp=Ts+(time.s*(qqq/60)) #temperatures calculated at each time
tm=time.s

SBf = function(tm, state, parms)  
{
with(as.list(c(tm, state, parms)),
{
   a1 = parms[["a1"]]
   dy1 =   A*  exp((-Ea)/(R*a1(tm)))  * (y1)^m * (1-y1)^n #this should do the trick
   return(list(c(dy1)))
   })
}
#tm = seq(0, 10, len = 100)
state = c(y1 = prec) #starting value
a1 = approxfun(tm,Temp) #function that changes in time
P = list(a1 = a1)
sol = ode(y = state, times = tm, parms = P, func = SBf) #it works but gives a flat result

plot(sol) #correct!

#df <- df[Reduce(`&`, lapply(df, is.finite)),] ###check values
T.C <- Temp-273.15
T.K <- Temp
sol2 <- sol[,2]
sol2[length(sol2)]=sol2[length(sol2)-1]+(abs(sol2[length(sol2)-1]-(sol2[length(sol2)-2]))) #to avoid NA
fi <- c(0,diff(sol2))
#alfa=sol[,2]*100
alfa=sol2*100
dadT=c(0,diff(sol2)/diff(T.C))
timef <- sol[,1]


my.list <- list("T.C" = T.C,"T.K"=T.K,"fi"=fi,"alfa"=alfa, "dadT"=dadT, "time.s"=timef,"sol"=sol)
return(my.list)
}
