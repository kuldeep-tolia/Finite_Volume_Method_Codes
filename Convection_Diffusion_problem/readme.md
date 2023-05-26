Problem Description:  

-> For this problem, consider a steady transport of a scalar $\phi$ by convection-diffusion equation inside a square domain of size $L$.  
-> The governing equation can be written as:

$$\nabla \cdot \left( \rho \mathbf{v} \phi \right) = \nabla \cdot \left( \Gamma \nabla \phi \right) + S $$

-> Boundary Conditions:
  - Left side => $\phi = 0$
  - Bottom side => $\phi = 1$
  - Rest of the sides are treated as outflow boundaries where an approximation of high Peclet number is made.  

-> The flow field is given as $\mathbf{v} = f(x, y)$. Any suitable form of velocity field can be used, for example $u = x^3 + 5$.    
-> Similarly, the source term can also be assumed as $S = f(x, y)$.  
-> Suitable values for $\rho$ and $\Gamma$ are assumed.  
-> FVM is used to discretise the governing equation on a uniform cartesian grid.  
-> To discretise the convection term in the governing equation, the following schemes are used:
  - Central Differencing Scheme (CDS)
  - QUICK scheme

-> The results are presented through contours of $\phi$ and variation of $\phi$ on the horizontal center line.  
