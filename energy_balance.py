# ENERGY BALANCE MODEL
#
# (incoming radiation from sun) ~ (outgoing energy from Earth)
#
# Q1: Can temperature be accurately predicted/explained over
#   historical timescales? 
# Q2: What do the steady-state solutions of a 1d (2d) energy-balance
#       model look like?
# A2: HW3 below.
# Q3: Can the bifurcation of steady-state solutions be quantitatively
#     described? See Cardy on renormalization and order parameters.
#     Also see Chen-Hale/Verhulst(Nayfeh?) on bifurcations,  
#
# Implementation: The Boltzmann equation evolves the densities
#   \begin{enumerate}
#       \item Mass: $\rho_a(x_\alpha,t_m) = a, \ldots, K$ where
#                $(\alpha,m)$ denotes a space-time grid point. Let
#                $f_{\alpha\beta m}^a \doteq \frac{1}{|Q_{\alpha}||\tau_m|}\int_{Q_\alpha \times \tau_m} f_0^a dx dv dt$ 
#                By testing the equation
#                \begin{align}\label{boltzmann_vlasov}
#                   \D_t F^a + v \cdot \nabla_x F^a - \nabla_x \Phi^a_\gamma \cdot \nabla_v F^a = Q(F^a, \sum_{b=1}^K F^b)
#                \end{align}
#                The left-hand side represents the transport of the densities $F^a$; in its linearized form it has an
#                explicit representation formula in terms of the backwards trajectories; see [Vidav70s], "Nuetron
#                Transport". By using this "double-duhamel" representation of the solution and the detailed 
#                decomposition of the linearized operator $L^a[f] = -\nu^a(x,v) \vec{f}(t,x,v) + K^a[\vec{f}]$ due to
#                Grad [40s-50s] (but see [Glassey96] for more recent treatment), we wish to prove estimates of the form:
#                \begin{align}\label{stability_est}
#                   \sup_{t\leq T_0} \norm{ F(t) - \mu }_{\L^p_w(\Omega \times \R^3)} 
#                   \leq C^1_{d,p,\Omega,H,T_0}\norm{ F_0 }_{\L^p_w(\Omega \times \R^3)}
#                       + C^2_{d,p,\Omega,H,T_0}\fancy H[F_0]
#                \end{align}
#                
#                   
#                In particular, guided by the Lyupanov functional and the weighted $L^\infty$ estimates that we expect to
#                hold for the nonlinear term look for solutions
#                \begin{align}\label{weighted_Lp_spaces}
#                   \L^p_{w_{H,\beta,\beta_1}}(\Omega \times \R^3) 
#                   =  \left( \int w_{H, \beta, beta_1} |f^a(t,x,v)|^p dx dv \right)^{1/p}
#                \end{align}
#                \begin{align}\label{conservation_laws}
#                   \fancy H[f] \doteq
#                \end{align}
#                Q: in the single-species or multi-species case, having $-\int F^a \log F^a$ low => ?
# Assume further a "moderate diurnal cycle" due to
#   \begin{enumerate}
#       \item Large heat capacity of earth:
#       \item Fast rotation of earth:
#   \end{enumerate}
# Note: there is also diurnal convective action transporting
# water vapor and, one assumes, CO2 (passive tracer). 
#
# At the scales at which particles interact with each other,
# intermolecular forces may be the primary organizing factor
# in the steady states (i.e. thermodynamic equilibria)
# you expect to see for a given system of fixed volume, particle
# numbers $N_1, \ldots, N_k$, (and so $\sum_k N_k \doteq N$ their
# total is also conserved), and bulk-conserved mass, momentum and
# energy.
#
# Task 1: Write down Hamiltonian of the interacting particle system:
# H(t,X,V) = \frac{1}{N} \sum_{i=1}^N |V_i|^N + \frac{\alpha}{N^2} \sum_{i,j=1}^N \phi_{ij}^{\gamma_{ij}}(X_i-X_j) 
# where
# 
# Q1: How is this energy conservation a restatement of the 1st law of thermo (energy balance)?
# Dependencies:
import numpy as np
import scipy
import matplotlib as plot
from shapely.geometry.polygon import Polygon
from itertools import permutations

    L = 10 # [km]
    
    # construct ordered sequence of point tuples: the 'shell' 
    # corners
    hull = [ [L,L],[-L,L], [-L,-L],[L,-L] ]
    omega = Polygon(hull)
    
    # test 1:
    p = plot.plot(*omega.exterior.xy)
    plot.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    



def vol(omega)
    return 
