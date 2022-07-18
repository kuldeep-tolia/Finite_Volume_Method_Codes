Problem Description:  

-> For this problem, consider a steady transport of a scalar $\phi$ by advection-diffusion equation inside a square domain of side $L$.  
-> Boundary Conditions:
  - Left side => $\phi = 0$
  - Bottom side => $\phi = 1$
  - Rest of the sides are treated as outflow boundaries where an approximation of high Peclet number is made.  

-> The flow field is given by $u = x^2 + 1, v = y^2 + 1$. Any other suitable velocity field can be used as well.  
-> Similarly, the source term can also be assumed as $S = f(x, y)$.  
-> Suitable values for $\rho$ and $\Gamma$ are assumed.  
-> The governing equation can be written as:

$$\nabla \cdot \left( \rho \vec{v}\phi \right) = \nabla \cdot \left( \Gamma \nabla \phi \right) + S$$

-> FVM is used to discretise the control volume on a structured cartesian grid.  
-> To discretise the convection term in the governing equation, following schemes are used:
  - Central Differencing Scheme (CDS)
  - QUICK scheme

-> The results are presented through contours of $\phi$ and variation of $\phi$ on the horizontal center line.  
