#' Title sbAC
#' @description performs Šesták-Berggren (AC) simulations
#' @param
#'        starting time for the simulations "time.start",
#'	      temperature start "T0",
#'		  end temperature "T.end",
#'		  heating rate "qqq",
#'		  Ea, A, m and n parameters in the equation
#'        number of points "npoints",
#'		  starting value for the equation "prec"
#'
#' @return startgin temperature "T","fi",degree of crystallization "alfa",differential alfa in T "dadT",time in seconds "time.s",differential equation solution "sol"
#' @export
#' @import deSolve
#' @references J. Šesták. Thermophysical Properties of Solids, Their Measurements and Theoretical Analysis. Elsevier: Amsterdam, 1984.
#' @examples  \donttest{
#'
#' res <- sbAC(npoints=100rm(list=ls())00,prec=10^(-4.30095790876))
#'
#' }

#library(deSolve)
sbAC <-  function (time.start=0,
                   T0=0, T.end=500,
                   qqq=50,
                   A=10^(6.3),
                   Ea=80000,
                   m=1,n=2,
                   npoints=10000,
                   prec=10^(-4.30095790876),
                   ...)
{
  R  <- 8.314       #gas constant
  Ts <- 273.15+ T0  #transformation in K

  time.e=(T.end - T0)/(qqq/60) #estimated end time for the analysis depending on the rate
  time.s=seq(time.start,time.e, length.out=npoints) #vector with all the times for solving the equations time in seconds    KARLINE: THIS HAS NO IMPACT ON THE SOLUTION...
  tm = time.s
  Temp = Ts+(time.s*(qqq/60)) #temperatures calculated at each time

####
#function to solve
SBf = function(tm, state, parms)
{
with(as.list(c(state, parms)),
{

   a1_tm = Ts+(tm*qqq/60)   # Karline: this can be calculated here based on the true time.

   dy1 =   A* exp(-Ea/(R*a1_tm))  * (y1)^m * (1-y1)^n #this should do the trick
   return(list(dy1,                 # first return value = derivative (used to integrate)
               Temp = a1_tm,        # second value = calculated temperature
               fi = dy1))           # the derivative again - now used for output
   })
}
#end of the function
####

state = c(y1 = prec) #
P <- list(Ts = Ts, qqq = qqq, A = A)
sol <- ode(y = state, times = tm, parms = P, func = SBf, ...)
###
plot(sol) #a visual check of the results obtained
###
T.C <- Temp-273.15
T.K <- Temp
timef <- sol[,1]
alfa <- sol[,2]*100
dadT=c(0,diff(sol[,2])/diff(T.C))
my.list <- list("T.C" = T.C, "T.K"=T.K,"sol"=sol,"fi"=sol[,4],"alfa"=alfa, "dadT"=dadT, "time.s"=timef)
return(my.list)
}



