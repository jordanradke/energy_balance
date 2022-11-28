import numpy as np
from energy_bal1d import *
## Task 2: Solve linear algebraic system to write an advance mapping.

# ADVANCE MAPPING
def advance(T, h, dt, cw, Qsa, A, B, D, x):
	N = T.shape[0]
	
	# stability criterion: invertibility (perturbation of identity 
	# is invertible)
	if dt/cw*np.linalg.norm( np.ones(N)*B - diffm(N,h,D,x), 1 ) < 0.5:
		M1 = np.eye(N) + dt/cw*( np.ones(N)*B - diffm(N,h,D,x) )
		M2 = T + dt/cw*(Qsa - np.ones(N)*A)
		return np.dot(np.linalg.inv(M1),M2)
	else:
		return "choose smaller dt:" + str(dt) + ". Off by " + str(dt/cw*np.linalg.norm( np.ones(N)*B - diffm(N,h,D,x), 1 ))


# Inputs:
#	N: shape of temp array
#	h: spatial resolution
#	D: diffusivity parameter
def diffm(N,h,D,x):
	
	diffm = np.zeros((N,N))
	for i in np.arange(0,N):
		e = np.eye(1,N,i).T	# std basis vector
		diffm[:,i] = diffusion(e, h, D, x)
	
	return diffm



# Inputs:
#	h: spatial grid res
#	dt: time grid res 
#	Tmax: length of simulation
#	(A,B,D,cw,S0,S2,a0,a2,Q): physical/empirical parameters
def sol_ie(h,dt,Tmax,A,B,D,cw,S0,S2,a0,a2,Q):
	# define space and time grid
	x = np.arange(0,1+h,h)
	N = x.shape[0]
	steps = int(dt**(-1)*Tmax)
	t  = np.linspace(0,Tmax,steps)
	
	# a forcing function evaluated on grid points. Note S0=1 here.	
	# co-alebedo and insolation are approximated with Legendre polynomials.
	P2 = lambda X: 1/2*(3*X**2-1)
	Qsa = Q*(S0+S2*P2(x))*(a0 + a2*P2(x))
	
	T0 = 10*np.ones(N) # initial temp. profile.
	
	# solve pde
	sol_ie = np.zeros((steps, N))	
	sol_ie[0,:] = T0
	for k in np.arange(0,steps-1):
		if k % 100 == 0:
			print("current time step:" + str(k) + "of " + str(steps))
		sol_ie[k+1,:] = advance(sol_ie[k,:], h, dt, cw, Qsa, A, B, D, x)
		
	return sol_ie



# plot results	
#plot.plot(t,sol_ie)
#plot.xlabel("time [years]")
#plot.ylabel("temp [C]")

#plot.show()







