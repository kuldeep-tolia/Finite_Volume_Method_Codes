Problem Description:  

-> A more general Couette flow includes a constant pressure gradient, $dP/dx$, in a direction parallel to the plates.  
-> The applied pressure gradient can be an adverse pressure gradient or a favourable pressure gradient.  
-> The size of the computational domain, normalised using the plate separation distance $h$, is taken as 2 x 1 in x (streamwise) and y (wall-normal) directions.  
-> No-slip boundary conditions are applied for the velocity on the plates and periodic boundary conditions are applied on the transverse sides of the domain.  
-> The pressure gradient term is treated as a body force term.  
-> The Reynolds number considered for this problem is $Re = U h / \nu = 15$.  
-> For this problem, a non-dimensional parameter is defined:  

$$ P = -\frac{h^2}{2 \mu U} \left( \frac{dP}{dx} \right) $$

where $U \equiv$ top plate velocity, $\mu \equiv$ dynamic viscosity of the fluid, and $h \equiv$ plate separation distance.

-> Various values of $P$ are tested and it corresponds to:
  - $P < 0$ : adverse pressure gradient
  - $P>0$ : favourable pressure gradient
  - $P=0$ : zero pressure gradient
