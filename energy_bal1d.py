import numpy as np
from scipy.integrate import odeint
## Task 1. simple solution with odeint

# rhs of a forced heat equation with diffusion given below.
# Note: f, T are (N+1)x1 vectors below. Is this problematic for
# def of forced_heat? The point is you need all (or, at least the
# nearest neighbors) of T in order to define the non-local-in-space
# operator diffusion(T(x,t)). 
def forced_heat(T, t, h, A, B, D, cw, Qsa, x):
	return 1/cw*(Qsa - linearized_sb(T,A,B) + diffusion(T,h,D,x))

# In -10 to 40 C range, not much lost by this linearization of
# \begin{align}\label{sb_law}
#	Inputs:
#		T: An (Nx1) vector
#		A: a scalar parameter
#		B: a scalar parameter
#	Output:
#		an (Nx1) vector of forcing terms
def linearized_sb(T,A,B):
	return A*np.ones(T.shape[0])+B*T
	

	
## FUNCTIONS TO CONSTRUCT PDE

# The spatially inhomogenous diffusion operator
#	Inputs:
#		T: (Nx1) temperature vector
#		h: spatial grid-scale
#		D: diffusion coefficient
#	Output:
#		\D_h[T], an (Nx1) vector 
#	NOTE: The boundary conditions are incorporated into diffusion matrix.
def diffusion(T,h,D,x):
	# T is a Nx1 vector, indices 0, 1, ..., N-1. 
	N = T.shape[0]
	diffT = np.zeros(N) # initialize as Nx1 vector.
	
	# At i=0, use boundary condition to rewrite diffusion operator.
	# We invoke a GHOST POINT at i=-1, where the no-flux condition
	# at equator, implemented via a leapfrog difference, is enforced:
	# \begin{align}
	#	\frac{\D T}{\D x}\Big|_{x=0} \sim D_h^{leap}[T][0] 
	# 	\doteq 1/(2*h)*(T[1]-T[-1])	= 0.
	# \end{align}
	# Therefore, this boundary condition introduces the constraint 
	# T[0] = T[2]. We may then use to re-write the only non-degenerate
	# portion of the diffusion operator $\D_h[T]$ at $i=0$, i.e., 
	# the second-order term of the central-difference operator approximating
	# the term $\frac{\D^2 T}{\D x^2}:
	# \begin{align}
	#	\D_h[T][0] = D*(T[1] - 2*T[0] + T[-1])/(h**2)	
	#				= 2*D*(T[1] - T[0])/(h**2)
	# \end{align}
	diffT[0] = 2*D*(T[1]-T[0])/(h**2)
	# Use boundary condition to write endpoint at pole. Since x[N] = 1,
	# at the pole, the only non-degenerate part of the diffusion operator
	# is the first-order term:
	#	\begin{align}
	#		\D[T](1) = (1-1^2)\frac{\D^2 T(1)}{\D x^2} - 2 \frac{\D T(1)}{\D x}
	#				= -2 \frac{\D T(1)}{\D x}
	#	\end{align}
	# which is implemented via a one-sided difference of O(h) accuracy.
	# does updating before or after main loop affect?
	diffT[-1] = -2*D*(T[-1] - T[-2])/h		
	# central differences for second-order term, leapfrog for first-order
	# for all interior points.
	for i in np.arange(1, N-1):
		diffT[i] = D*( (1-x[i]**2)*cdiff2(T,i,h) - 2*x[i]*leap(T,i,h) )
	
	return diffT
	

#	Inputs:
#		u: a (Nx1) vector.
#		i: index at which scheme is cenetered
#		h: grid resolution
def cdiff2(u,i,h):
    return (u[i+1] - 2*u[i] + u[i-1])/(h**2) 


def leap(u,i,h):
	return (u[i+1] - u[i-1])/(2*h)
	

def sol_odeint(h,dt,Tmax,A,B,D,cw,S0,S2,a0,a2,Q):
	# define space and time grid
	x = np.arange(0,1+h,h)
	N = x.shape[0]
	t  = np.linspace(0,Tmax,int(dt**(-1)*Tmax) )
	
	# a forcing function evaluated on grid points. Note S0=1 here.	
	# co-alebedo and insolation are approximated with Legendre polynomials.	
	P2 = lambda X: 1/2*(3*X**2-1)
	Qsa = Q*(S0+S2*P2(x))*(a0 + a2*P2(x))
	
	T0 = 10*np.ones(N) # initial temp. profile.
	
	sol = odeint(forced_heat, T0, t, args=(h,A,B,D,cw,Qsa,x))
	
	return sol 


# plot results: below, the time evolution for each zonal averaged strip k, T[:,k],
# subjecto to T0[k] = 0 for all $k = 0, \ldots, N$.
#plot.plot(t,sol)
#plot.xlabel("time [years]")
#plot.ylabel("temp [C]")

# albedo and solar irradiance spatial distributions (meridional coordinate):
#fig, (ax1, ax2) = plot.subplots(1, 2)

#ax1.plot(x, a0 + a2*P2(x))
#ax1.set(xlabel="sin(latitude)", ylabel = "albedo estimate")
#ax2.plot(x, S0 + S2*P2(x))
#ax2.set( xlabel="sin(latitude)", ylabel="solar irradiance distribution")
#fig.tight_layout(pad=5.0)
#plot.show()














