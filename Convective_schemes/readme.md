Problem Description:  

-> For this problem, consider a steady transport of a scalar $\phi$ by convection process inside a square domain of size $1$.  
-> The governing equation to describe the given problem is

$$\nabla \cdot \left( \rho \mathbf{u} \phi \right) = 0 $$

-> Boundary Conditions:
  - Left boundary => $\phi = 0$
  - Bottom boundary => $\phi = 1$
  - Rest of the boundaries are treated as outflow boundaries where an approximation of high Peclet number is made.  

-> The flow field is given as $\mathbf{v} = 1 \hat{i} + 1 \hat{j}$.    
-> FVM is used to discretise the governing equation on a uniform cartesian grid.  
-> To discretise the convection term, the following schemes are used:
  - Central Differencing Scheme (CDS) $\left( \beta = 1.0 \right)$
  - $1^{st}-$ order upwind scheme $\left( \beta = 0.0 \right)$
  - Hybrid scheme with the blending factor $\left( \beta = 0.9 \right)$ 

-> The results are presented through contours of $\phi$ and variation of $\phi$ on the vertical center line.  
