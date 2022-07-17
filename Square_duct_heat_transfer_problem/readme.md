Problem Description:  
  
-> For this problem, consider a flow that is hydrodynamically and thermally fully developed inside a square duct of side $L$.  
-> The walls are imposed with a constant heat flux $q$.  
-> For simplicity, the velocity profile inside the pipe is assumed to be uniform, i.e. $u(x, z) = U$.  
-> All thermal and hydrodynamic properties (like $k$, $C_p$, $\rho$) are assumed to be constant.  
-> The diffusion in the $y$-direction is neglected.  
-> The simplified governing equation can be written as:  

$$ \nabla \cdot \left( k \nabla T \right) = \displaystyle \frac{\partial \left( \rho U C_p T\right)}{\partial y}$$  

-> FVM is used to discretise the control volume on a structured cartesian grid.  
-> The results are plotted using non-dimensional temperature and Nusselt number.  
