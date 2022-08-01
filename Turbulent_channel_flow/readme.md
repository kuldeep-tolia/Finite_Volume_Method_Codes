Problem Description:

-> For this problem, consider a fully-developed, incompressible turbulent channel flow.   
-> The mean flow is one-dimensional (1D).  
-> Assume that the flow is statistically stationary (i.e. $\frac{\partial <>}{\partial t} = 0$) and statistically homogenous in spanwise direction (i.e. $\frac{\partial <>}{\partial z} = 0$), where $<>$ represents ensemble-average.  
-> Since the flow is fully developed, we have $V = \frac{\partial U}{\partial x} = \frac{\partial k}{\partial x} = 0$.  
-> Thus the governing equations are reduced to:

$$0 = -\frac{1}{\rho} \frac{\partial P}{\partial x} + \frac{\partial}{\partial y} \left( \left( \nu + \nu_t \right) \frac{\partial U}{\partial y} \right)$$

$$0 = \frac{\partial}{\partial y} \left( \left( \nu + \frac{\nu_t}{\sigma_k} \right) \frac{\partial k}{\partial y} \right) + P_k - \varepsilon$$

-> The equations are reduced to 1D diffusion equations with complex source terms and discretised using Finite Volume Method.  
-> The flow is driven by the pressure gradient $\frac{\partial P}{\partial x}$ and is balanced by the wall shear stress.  
-> Here, I have assumed a fully-developed channel flow with $Re_{\tau} = 395$.  
-> For the sake of simplicity, I have assumed the channel half height $\delta = 1$ and $\rho = 1$.  
-> The obtained results are compared with the reference data.  

Reference: R. D. Moser, J. Kim, and N. N. Mansour. Direct numerical simulation of turbulent channel flow up to $Re_{\tau}$ =590. Physics of Fluids, 11(4):943–945, 1999.
