import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plot
from energy_bal1d import *
from energy_bal1dIE import *
from analytic1d import sol_analytic


# Inputs:
#	h: space grid resolution [m]	
#	dt: time grid resolution [years]
#	Tmax: maximum time [years]
#	(A,B,D,cw,S0,S2,a0,a2): fixed parameters
#	ie: boolean value; if True, implicit euler (ie) chosen over odeint soln
#
# Dependencies: sol_ie, sol_odeint, functions defined in energy_bal1d.py and
#				energy_bal1dIE.py respectively.		
def error_benchmark(h,dt,Tmax,A,B,D,cw,S0,S2,a0,a2,Q,ie):

	x = np.arange(0,1+h,h)
	N = x.shape[0]
	steps = int(dt**(-1)*Tmax)
	t  = np.linspace(0,Tmax,steps)
	
	ub = sol_analytic(h,a0,a2,S2,Q,A,B,D)
	if ie:
		u = sol_ie(h,dt,Tmax,A,B,D,cw,S0,S2,a0,a2,Q)
	else:
		u = sol_odeint(h,dt,Tmax,A,B,D,cw,S0,S2,a0,a2,Q) 
	
	fig, (ax1, ax2) = plot.subplots(1, 2)

	ax1.plot(x, ub)
	ax1.set(xlabel="sin(latitude)", ylabel = "temp [C]")
	ax2.plot(t, u) 
	ax2.set( xlabel="time [years]", ylabel="temp [C]")
	fig.tight_layout(pad=5.0)
	plot.show()
	
	return np.max(np.abs(u[-1,:]-ub))


### MAIN

## INITIALIZE PARAMETERS
Q = 340         # solar constant/4 = 1360/4 (W m^-2). I.e. 
				# "mean solar electromagnetic radiation/unit area" or 
				# "total solar irradiance". Average makes sense to consider
				# Q: Why 1/4th taken here?
A = 203.3       # Outgoing Long-wave Radiation (OLR) when T(t,x) = 0 (W m^-2)
B = 2.09        # OLR temperature dependence (W m^-2 K^-1)
a0 = 0.681      # 1. Legendre Poly. albedo cfft. An order 2 expansion.
a2 = -0.202     # 2. Legendre Poly. albedo cfft
S0 = 1			# \int(S0)=1 while int(sm) = int(S0,Sm) = 0 for all other m.
S2 = -0.477     # solar forcing value in NCC81 (W m^-2)
cw = 6.3        # ocean mixed layer heat capacity (W yr m^-2 K^-1)
D = 0.649       # diffusivity for heat transport (W m^-2 K^-1)


P2 = lambda X: 1/2*(3*X**2-1)

# the grid points: 0 1/n ... .... (n-1)/n 1
#n = 10 # number of cells
#h = 1.0/n
#x = np.arange(0,1+h,h)
# N = n+1 , number of grid points. All functions below defined on this grid.
#N = x.shape[0]	

Tmax = 30

n = 2**np.arange(2,5)
err = np.zeros(n.shape[0])
## for each grid scale, calculate error between benchmark solution and 
#	a. solution found with ODEINT
#	b. solution found with implicit euler scheme
for j in np.arange(0,n.shape[0]):	
	h = 1/n[j]	# cell width
	x = np.arange(0,1+h,h)
	
	# dt chosen small enough for invertibility of advance matrix
	dt = (cw/10)*np.max(np.abs(B - diffm(n[j]+1,h,D,x)))**(-1)
	ie = True
	
	err[j] = error_benchmark(h,dt,Tmax,A,B,D,cw,S0,S2,a0,a2,Q,ie)
	

plot.plot(n,err)
plot.show()
	
	

