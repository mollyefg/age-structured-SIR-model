
require(deSolve);

### modeling the effects of age structure

# four classes: susceptible children, infected children
#		    susceptible adults, infected adults

#The function that describes the ODE's
Age<-function(t,y,p){
	Sc = y[1];
	Ic = y[2];
	Sa = y[3]
	Ia = y[4]
	with(as.list(p), {
		dSc.dt = nu*(1-f)-Sc*(beta[1,1]*Ic + beta[2]*Ia)-lc*Sc
		dIc.dt = Sc*(beta[1]*Ic + beta[2]*Ia) - gamma*Ic - lc*Ic
		dSa.dt = lc*Sc - Sa*(beta[3]*Ic + beta[4]*Ia) - muA*Sa
		dIa.dt = lc*Ic + Sa*(beta[3]*Ic + beta[4]*Ia) - gamma*Ia - muA*Ia
		return(list(c(dSc.dt, dIc.dt, dSa.dt, dIa.dt)));
	})
}

#################################################################################

# defining the parameters

beta = matrix(c(100, 10, 10, 20), nrow = 2, ncol = 2)		# infection rate, specific to age class interactions
nC = .2	# number of children (fraction of the population)
nA = 1-nC	# number of adults
lc = .0667	# growth rate from child to adult
gamma = 10	# death rate due to disease
f = .1	# fraction vaccinated
muA = .0167	# adult-specific death rate
nu = lc*nA	# birth rate in absence of disease

p = list(beta=beta, gamma=gamma, lc = lc, f = f, muA = muA, nu = nu)

# setting the initial conditions

Sc0 = .1
Ic0 = .001
Sa0 = .3
Ia0 = .0001
y0 = c(Sc0, Ic0, Sa0, Ia0)
Steps  = 10
t = seq(from=0,to=600,by=1/Steps) # defining the time interval and step size

out = ode(y=y0,times=t,func=Age,parms=p); # doing the integration

par(mfrow = c(2,1))

#graph one: fraction susceptible
plot(out[,1],(out[,2])/nC,type="l",lwd=1,main="SIR model with age structure", xlab="Time",ylab="Fraction Susceptible", ylim = c(min(out[,4]/nA), max(out[,2]/nC)));
lines(out[,1], (out[,4])/nA, lwd = 2)

legend(300, .4, legend = c("Children", "Adults"), lwd = c(1,2), bty = "n")

#graph two: fraction infected
plot(out[,1],((out[,3])/nC),type="l",lwd=1,xlab="Time",ylab="Fraction Infected", ylim = c(min(out[,5]/nA), max(out[,3]/nC)) );
lines(out[,1], ((out[,5])/nA), lwd = 2)
End = nrow(out)

#print relevant output to the consoel
cat("Fraction susc. adults:", out[End,4]/nA, "\n")


