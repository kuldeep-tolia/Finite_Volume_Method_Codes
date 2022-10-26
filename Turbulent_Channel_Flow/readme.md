Problem Description:

-> Consider a fully-developed, incompressible turbulent channel flow.  
-> The mean flow is 1D.  
-> Assume that the flow is statistically stationary (i.e. $\frac{\partial <>}{\partial t} = 0$) and spanwise homogenous (i.e. $\frac{\partial <>}{\partial z} = 0$), where $<>$ represents an ensemble averaged quantity.  
-> Since the flow is fully developed, we have $V = \frac{\partial U}{\partial x} = \frac{\partial k}{\partial x} = 0$.  
-> Thus the governing equations are reduced to:  
$$0 = -\frac{1}{\rho}\frac{\partial P}{\partial x} + \frac{\partial}{\partial y} \left( \left( \nu + \nu_{t} \right) \frac{\partial U}{\partial y} \right)$$  
$$0 = \frac{\partial}{\partial y} \left( \left( \nu + \frac{\nu_{t}}{\sigma_k} \right) \frac{\partial k}{\partial y} \right) + P_k - \varepsilon$$  
-> The channel has a height of $y_{max} = 2 \delta$. The flow is driven by the pressure gradient, $\partial P/ \partial x$, and is balanced by the wall shear stress since the flow is fully developed.  
-> The governing friction Reynolds number is assumed to be $Re_{\tau} = 395$.  
-> The results are compared with the DNS data available in the literature.  

Reference: R. D. Moser, J. Kim, and N. N. Mansour. Direct numerical simulation of turbulent channel flow up to $Re_{\tau}=590$. Physics of Fluids, 11(4):943–945, 1999.
