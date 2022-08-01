Problem Description:

-> For this problem, consider a fully-developed, incompressible turbulent channel flow.   
-> The mean flow is one-dimensional (1D).  
-> Assume that the flow is statistically stationary (i.e. $\frac{\partial <>}{\partial t} = 0$) and statistically homogenous in spanwise direction (i.e. $\frac{\partial <>}{\partial z} = 0$), where $<>$ represents ensemble-average.  
-> Since the flow is fully developed, we have $V = \frac{\partial U}{\partial x} = \frac{\partial k}{\partial x} = 0$.  
-> Thus the governing equations are reduced to:

$$0 = -\frac{1}{\rho} \frac{\partial P}{\partial x} + \frac{\partial}{\partial y} \left( \left( \nu + \nu_t \right) \frac{\partial U}{\partial y} \right)$$

$$0 = \frac{\partial}{\partial y} \left( \left( \nu + \frac{\nu_t}{\sigma_k} \right) \frac{\partial k}{\partial y} \right) + P_k - \varepsilon$$

-> The equations are reduced to 1D diffusion equations with complex source terms.  
-> The flow is driven by the pressure gradient ($\frac{\partial P}{\partial x}$) and is balanced by the wall shear stress.  
-> s
