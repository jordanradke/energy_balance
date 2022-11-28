import numpy as np
## Analytic solution (testing)
def sol_analytic(h,a0,a2,S2,Q,A,B,D):
	x = np.arange(0,1+h,h)
	
	H0 = a0 + a2*S2/5
	H2 = a0*S2 * a2 + a2*S2*2/7
	H4 = 18/35*a2*S2
	T0 = (Q*H0 - A)/B
	T2 = Q*H2/(6*D+B)
	T4 = Q*H4/(20*D+B)
	P2 = lambda X: 1/2*(3*X**2-1)
	P4 = lambda X: 1/8*(3-30*X**2 + 35*X**4)
	T  = lambda X: T0 + T2*P2(X) + T4*P4(X)
	
	return T(x)
	


# from till's solutions.
def odefunc(T,t,h,A,B,D):
	
	N = T.shape[0]
	
	diffT[0] = D*2*(T[1]-T[0])/(h**2)		# equator b.c.
	diffT[-1] = -D*2*x[-1]*(T[-1]-T[-2])/h	# pole b.c.
	
	for i in range(i,N-1):
		diffT[i] = D/(h**2)*(1-x[i]**2)*(T[i+1]-2*T[i]+T[i-1]) -(D*x[i]/h)*(T[i+1]-T[i-1])
		
	return 1/cw*(Qsa - (A + B*T) + diffT)
 
