!Code for Steady Convection-Diffusion equation using QUICK scheme for convective term

module functions

implicit none

	integer, parameter:: nx = 50, ny = 50				!Set grid size here; ny-->refers horizontal direction, nx-->refers vertical direction
	integer:: i, j, iter, t
	double precision, dimension(:, :), allocatable:: aP, aE, aW, aN, aS, b, phi, phi_old
	double precision, dimension(:, :), allocatable:: De, Dw, Dn, Ds, Fe, Fw, Fn, Fs
	double precision, dimension(:, :), allocatable:: del_e, del_w, del_n, del_s
	double precision, dimension(:), allocatable:: x_f, y_f, x_c, y_c, x_center, phi_center
	double precision:: dx, dy, L, H, rho, diff_coeff, phi_a, phi_b, phi_guess
	double precision:: alpha, tolerance, phil2norm, numerator, denominator
	character (len = 256):: nxx, nyy
	character (len = 1024):: fname1, fname2, fname3
	
contains

	subroutine allocate_memory()
	
	implicit none
	
		allocate (aP(1:nx, 1:ny), aE(1:nx, 1:ny), aW(1:nx, 1:ny), aN(1:nx, 1:ny), aS(1:nx, 1:ny), b(1:nx, 1:ny), phi(1:nx, 1:ny))
		allocate (De(1:nx, 1:ny), Dw(1:nx, 1:ny), Dn(1:nx, 1:ny), Ds(1:nx, 1:ny), Fe(1:nx, 1:ny), Fw(1:nx, 1:ny), Fn(1:nx, 1:ny))
		allocate (phi_old(1:nx, 1:ny), Fs(1:nx, 1:ny))
		allocate (del_e(1:nx, 1:ny), del_w(1:nx, 1:ny), del_n(1:nx, 1:ny), del_s(1:nx, 1:ny))
		allocate (x_f(0:ny), y_f(0:nx), x_c(ny), y_c(nx))
		allocate (x_center(0:ny+1), phi_center(0:ny+1))
		
	end subroutine allocate_memory
	
	subroutine deallocate_memory()
	
	implicit none
	
		deallocate (aP, aE, aW, aN, aS, b, phi, phi_old)
		deallocate (De, Dw, Dn, Ds, Fe, Fw, Fn, Fs)
		deallocate (del_e, del_w, del_n, del_s)
		deallocate (x_f, y_f, x_c, y_c, x_center, phi_center)
	
	end subroutine deallocate_memory
	
	subroutine set_parameters()
	
	implicit none
	
		L = 3.0d0			!Length of domain
		H = 3.0d0			!Height of domain
		rho = 2.0d0			!Density of fluid
		diff_coeff = 3.0d0		!Diffusion coefficient
		
		phi_a = 1.0d0			!Bottom boundary phi value
		phi_b = 0.0d0			!Left boundary phi value
		
		phi_guess = 0.0d0		!Initial guess for phi_field
		
		alpha = 7.0d-1			!Relaxation factor
		tolerance = 1.0d-6		!Tolerance criteria
		phil2norm = 1.0d0		!Initial value of error
		iter = 1			!Initial counter value		
	
	end subroutine set_parameters
	
	subroutine calc_grid_size()
	
	implicit none
	
		dx = L / ny
		dy = H / nx
		
	end subroutine calc_grid_size
	
	subroutine calc_faces()
	
	implicit none
	
		x_f(0) = 0.0d0			
		y_f(0) = 3.0d0			
		
		!Populate x at faces
		do i = 1, ny
			x_f(i) = x_f(i-1) + dx
		end do
		
		!Populate y at face
		do i = 1, nx
			y_f(i) = y_f(i-1) - dy
		end do			
	
	end subroutine calc_faces
	
	subroutine calc_centers()
	
	implicit none
	
		x_c(1) = x_f(0) + (5.0d-1 * dx)
		y_c(1) = y_f(0) - (5.0d-1 * dy)
		
		!Populate x at centers
		do i = 2, ny
			x_c(i) = x_c(i-1) + dx
		end do
		
		!Populate y at centers
		do i = 2, nx
			y_c(i) = y_c(i-1) - dy
		end do
	
	end subroutine calc_centers
	
	subroutine initialise_solver()
	
	implicit none
	
		call calc_grid_size()
		call calc_faces()
		call calc_centers()
	
	end subroutine initialise_solver
	
	subroutine calc_interior_coeff()
	
	implicit none
	
		do i = 2, nx-1
			do j = 2, ny-1
				!Calculating diffusion coefficients				
				De(i, j) = diff_coeff
				Dw(i, j) = diff_coeff
				Dn(i, j) = diff_coeff
				Ds(i, j) = diff_coeff
				
				!Calculating convection coefficients
				Fe(i, j) = rho * dy * (1.0d0 + x_f(j)**2)
				Fw(i, j) = rho * dy * (1.0d0 + x_f(j-1)**2)
				Fn(i, j) = rho * dx * (1.0d0 + y_f(i-1)**2)
				Fs(i, j) = rho * dx * (1.0d0 + y_f(i)**2)
				
				!Calculating neighbouring coefficients
				aE(i, j) = De(i, j)
				aW(i, j) = Dw(i, j) + Fw(i, j)
				aN(i, j) = Dn(i, j)
				aS(i, j) = Ds(i, j) + Fs(i, j)
				
				!Calculating aP coefficient
				aP(i, j) = aE(i, j) + aW(i, j) + aN(i, j) + aS(i, j) + Fe(i, j) - Fw(i, j) + Fn(i, j) - Fs(i, j)
				
				!Calculating b term
				if ((j - 2) .le. 0) then
					del_e(i, j) = ((phi(i, j+1) - phi(i, j-1)) / 4.0d0) + ((phi(i, j+1) + phi(i, j-1) - (2.0d0 * &
							& phi(i, j))) / 8.0d0)
					del_w(i, j) = 0.0d0
					del_n(i, j) = ((phi(i-1, j) - phi(i+1, j)) / 4.0d0) + ((phi(i-1, j) + phi(i+1, j) - (2.0d0 * &
							& phi(i, j))) / 8.0d0)
					del_s(i, j) = ((phi(i, j) - phi(i+2, j)) / 4.0d0) + ((phi(i, j) + phi(i+2, j) - (2.0d0 * &
							& phi(i+1, j))) / 8.0d0)
					b(i, j) = (2.0d0 * dx * dy * (x_c(j) + y_c(i))) - (Fe(i, j) * del_e(i, j)) - (Fn(i, j) * &
							& del_n(i, j)) + (Fs(i, j) * del_s(i, j))
			
				elseif ((i + 2) .ge. (nx + 1)) then
					del_e(i, j) = ((phi(i, j+1) - phi(i, j-1)) / 4.0d0) + ((phi(i, j+1) + phi(i, j-1) - (2.0d0 * &
							& phi(i, j))) / 8.0d0)
					del_w(i, j) = ((phi(i, j) - phi(i, j-2)) / 4.0d0) + ((phi(i, j) + phi(i, j-2) - (2.0d0 * &
							& phi(i, j-1))) / 8.0d0)
					del_n(i, j) = ((phi(i-1, j) - phi(i+1, j)) / 4.0d0) + ((phi(i-1, j) + phi(i+1, j) - (2.0d0 * &
							& phi(i, j))) / 8.0d0)
					del_s(i, j) = 0.0d0
					b(i, j) = (2.0d0 * dx * dy * (x_c(j) + y_c(i))) - (Fe(i, j) * del_e(i, j)) - (Fn(i, j) * &
							& del_n(i, j)) + (Fw(i, j) * del_w(i, j))
				
				
				elseif (((j-2) .le. 0) .and. ((i+2) .ge. (nx+1))) then
					del_e(i, j) = ((phi(i, j+1) - phi(i, j-1)) / 4.0d0) + ((phi(i, j+1) + phi(i, j-1) - (2.0d0 * &
							& phi(i, j))) / 8.0d0)
					del_w(i, j) = 0.0d0
					del_n(i, j) = ((phi(i-1, j) - phi(i+1, j)) / 4.0d0) + ((phi(i-1, j) + phi(i+1, j) - (2.0d0 * &
							& phi(i, j))) / 8.0d0)
					del_s(i, j) = 0.0d0
					b(i, j) = (2.0d0 * dx * dy * (x_c(j) + y_c(i))) - (Fe(i, j) * del_e(i, j)) - (Fn(i, j) * &
							& del_n(i, j))
				
				else
					del_e(i, j) = ((phi(i, j+1) - phi(i, j-1)) / 4.0d0) + ((phi(i, j+1) + phi(i, j-1) - (2.0d0 * &
							& phi(i, j))) / 8.0d0)
					del_w(i, j) = ((phi(i, j) - phi(i, j-2)) / 4.0d0) + ((phi(i, j) + phi(i, j-2) - (2.0d0 * &
							& phi(i, j-1))) / 8.0d0)
					del_n(i, j) = ((phi(i-1, j) - phi(i+1, j)) / 4.0d0) + ((phi(i-1, j) + phi(i+1, j) - (2.0d0 * &
							& phi(i, j))) / 8.0d0)
					del_s(i, j) = ((phi(i, j) - phi(i+2, j)) / 4.0d0) + ((phi(i, j) + phi(i+2, j) - (2.0d0 * &
							& phi(i+1, j))) / 8.0d0)
					b(i, j) = (2.0d0 * dx * dy * (x_c(j) + y_c(i))) - (Fe(i, j) * del_e(i, j)) - (Fn(i, j) * &
							& del_n(i, j)) + (Fs(i, j) * del_s(i, j)) + (Fw(i, j) * del_w(i, j))
				end if				
			end do
		end do
					
	end subroutine calc_interior_coeff
	
	subroutine calc_south_boundary_coeff()
	
	implicit none
	
		do j = 2, ny-1
			!Calculating diffusion coefficients				
			De(nx, j) = diff_coeff
			Dw(nx, j) = diff_coeff
			Dn(nx, j) = diff_coeff
			Ds(nx, j) = 2.0d0 * diff_coeff
			
			!Calculating convection coefficients
			Fe(nx, j) = rho * dy * (1.0d0 + x_f(j)**2)
			Fw(nx, j) = rho * dy * (1.0d0 + x_f(j-1)**2)
			Fn(nx, j) = rho * dx * (1.0d0 + y_f(nx-1)**2)
			Fs(nx, j) = rho * dx * (1.0d0 + y_f(nx)**2)
			
			!Calculating neighbouring coefficients
			aE(nx, j) = De(nx, j)
			aW(nx, j) = Dw(nx, j) + Fw(nx, j)
			aN(nx, j) = Dn(nx, j)
			aS(nx, j) = 0.0d0
			
			!Calculating b term
			if ((j - 2) .ge. 1) then 
				del_e(nx, j) = ((phi(nx, j+1) - phi(nx, j-1)) / 4.0d0) + ((phi(nx, j+1) + phi(nx, j-1) - (2.0d0 * &
							& phi(nx, j))) / 8.0d0)
				del_w(nx, j) = ((phi(nx, j) - phi(nx, j-2)) / 4.0d0) + ((phi(nx, j) + phi(nx, j-2) - (2.0d0 * &
							& phi(nx, j-1))) / 8.0d0)
				del_n(nx, j) = 0.0d0
				del_s(nx, j) = 0.0d0
				b(nx, j) = (2.0d0 * dx * dy * (x_c(j) + y_c(nx))) + (Ds(nx, j) * phi_a) + (Fs(nx, j) * phi_a) - &
						& (Fe(nx, j) * del_e(nx, j)) + (Fw(nx, j) * del_w(nx, j))
			
			else
				del_e(nx, j) = ((phi(nx, j+1) - phi(nx, j-1)) / 4.0d0) + ((phi(nx, j+1) + phi(nx, j-1) - (2.0d0 * &
							& phi(nx, j))) / 8.0d0)
				del_w(nx, j) = 0.0d0
				del_n(nx, j) = 0.0d0
				del_s(nx, j) = 0.0d0
				b(nx, j) = (2.0d0 * dx * dy * (x_c(j) + y_c(nx))) + (Ds(nx, j) * phi_a) + (Fs(nx, j) * phi_a) - &
						& (Fe(nx, j) * del_e(nx, j))	
			end if		
				
			!Calculating aP coefficient
			aP(nx, j) = aE(nx, j) + aW(nx, j) + aN(nx, j) + Fe(nx, j) - Fw(nx, j) + Fn(nx, j) + Ds(nx, j)
		end do
	
	end subroutine calc_south_boundary_coeff
	
	subroutine calc_west_boundary_coeff()
	
	implicit none
	
		do i = 2, nx-1
			!Calculating diffusion coefficients				
			De(i, 1) = diff_coeff
			Dw(i, 1) = 2.0d0 * diff_coeff
			Dn(i, 1) = diff_coeff
			Ds(i, 1) = diff_coeff
			
			!Calculating convection coefficients
			Fe(i, 1) = rho * dy * (1.0d0 + x_f(1)**2)
			Fw(i, 1) = rho * dy * (1.0d0 + x_f(0)**2)
			Fn(i, 1) = rho * dx * (1.0d0 + y_f(i-1)**2)
			Fs(i, 1) = rho * dx * (1.0d0 + y_f(i)**2)
			
			!Calculating neighbouring coefficients
			aE(i, 1) = De(i, 1)
			aW(i, 1) = 0.0d0
			aN(i, 1) = Dn(i, 1)
			aS(i, 1) = Ds(i, 1) + Fs(i, 1)
			
			!Calculating b term
			if ((i + 2) .le. nx) then
				del_e(i, 1) = 0.0d0
				del_w(i, 1) = 0.0d0
				del_n(i, 1) = ((phi(i-1, 1) - phi(i+1, 1)) / 4.0d0) + ((phi(i-1, 1) + phi(i+1, 1) - (2.0d0 * &
							& phi(i, 1))) / 8.0d0)
				del_s(i, 1) = ((phi(i, 1) - phi(i+2, 1)) / 4.0d0) + ((phi(i, 1) + phi(i+2, 1) - (2.0d0 * &
							& phi(i+1, 1))) / 8.0d0)
				b(i, 1) = (2.0d0 * dx * dy * (x_c(1) + y_c(i))) + (Dw(i, 1) * phi_b) + (Fw(i, 1) * phi_b) - &
						& (Fn(i, 1) * del_n(i, 1)) + (Fs(i, 1) * del_s(i, 1))
			
			else 
				del_e(i, 1) = 0.0d0
				del_w(i, 1) = 0.0d0
				del_n(i, 1) = ((phi(i-1, 1) - phi(i+1, 1)) / 4.0d0) + ((phi(i-1, 1) + phi(i+1, 1) - (2.0d0 * &
							& phi(i, 1))) / 8.0d0)
				del_s(i, 1) = 0.0d0
				b(i, 1) = (2.0d0 * dx * dy * (x_c(1) + y_c(i))) + (Dw(i, 1) * phi_b) + (Fw(i, 1) * phi_b) - &
						& (Fn(i, 1) * del_n(i, 1))
			end if
				
			!Calculating aP coefficient
			aP(i, 1) = aE(i, 1) + aN(i, 1) + aS(i, 1) + Fe(i, 1) + Fn(i, 1) - Fs(i, 1) + Dw(i, 1)
		end do
	
	end subroutine calc_west_boundary_coeff
	
	subroutine calc_north_boundary_coeff()
	
	implicit none
	
		do j = 2, ny-1
			!Calculating diffusion coefficients				
			De(1, j) = diff_coeff
			Dw(1, j) = diff_coeff
			Dn(1, j) = 0.0d0
			Ds(1, j) = diff_coeff
			
			!Calculating convection coefficients
			Fe(1, j) = rho * dy * (1.0d0 + x_f(j)**2)
			Fw(1, j) = rho * dy * (1.0d0 + x_f(j-1)**2)
			Fn(1, j) = rho * dx * (1.0d0 + y_f(0)**2)
			Fs(1, j) = rho * dx * (1.0d0 + y_f(1)**2)
			
			!Calculating neighbouring coefficients
			aE(1, j) = De(1, j) 
			aW(1, j) = Dw(1, j) + Fw(1, j)
			aN(1, j) = 0.0d0
			aS(1, j) = Ds(1, j) + Fs(1, j)
			
			!Calculating b term
			if ((j - 2) .ge. 1) then 
				del_e(1, j) = ((phi(1, j+1) - phi(1, j-1)) / 4.0d0) + ((phi(1, j+1) + phi(1, j-1) - (2.0d0 * &
							& phi(1, j))) / 8.0d0)
				del_w(1, j) = ((phi(1, j) - phi(1, j-2)) / 4.0d0) + ((phi(1, j) + phi(1, j-2) - (2.0d0 * &
							& phi(1, j-1))) / 8.0d0)
				del_n(1, j) = 0.0d0
				del_s(1, j) = ((phi(1, j) - phi(1+2, j)) / 4.0d0) + ((phi(1, j) + phi(1+2, j) - (2.0d0 * &
							& phi(1+1, j))) / 8.0d0)
				b(1, j) = (2.0d0 * dx * dy * (x_c(j) + y_c(1))) - (Fe(1, j) * del_e(1, j)) + (Fw(1, j) * &
						del_w(1, j)) + (Fs(1, j) * del_s(1, j))
			
			else
				del_e(1, j) = ((phi(1, j+1) - phi(1, j-1)) / 4.0d0) + ((phi(1, j+1) + phi(1, j-1) - (2.0d0 * &
							& phi(1, j))) / 8.0d0)
				del_w(1, j) = 0.0d0
				del_n(1, j) = 0.0d0
				del_s(1, j) = ((phi(1, j) - phi(1+2, j)) / 4.0d0) + ((phi(1, j) + phi(1+2, j) - (2.0d0 * &
							& phi(1+1, j))) / 8.0d0)
				b(1, j) = (2.0d0 * dx * dy * (x_c(j) + y_c(1))) - (Fe(1, j) * del_e(1, j)) + (Fs(1, j) * del_s(1, j)) 
			end if				
				
			!Calculating aP coefficient
			aP(1, j) = aE(1, j) + aW(1, j) + aS(1, j) + Fe(1, j) - Fw(1, j) + Fn(1, j) - Fs(1, j)
		end do
	
	end subroutine calc_north_boundary_coeff
	
	subroutine calc_east_boundary_coeff()
	
	implicit none
	
		do i = 2, nx-1
			!Calculating diffusion coefficients				
			De(i, ny) = 0.0d0
			Dw(i, ny) = diff_coeff
			Dn(i, ny) = diff_coeff
			Ds(i, ny) = diff_coeff
			
			!Calculating convection coefficients
			Fe(i, ny) = rho * dy * (1.0d0 + x_f(ny)**2)
			Fw(i, ny) = rho * dy * (1.0d0 + x_f(ny-1)**2)
			Fn(i, ny) = rho * dx * (1.0d0 + y_f(i-1)**2)
			Fs(i, ny) = rho * dx * (1.0d0 + y_f(i)**2)
			
			!Calculating neighbouring coefficients
			aE(i, ny) = 0.0d0
			aW(i, ny) = Dw(i, ny) + Fw(i, ny) 
			aN(i, ny) = Dn(i, ny) 
			aS(i, ny) = Ds(i, ny) + Fs(i, ny) 
			
			!Calculating b term
			if ((i + 2) .le. nx) then
				del_e(i, ny) = 0.0d0
				del_w(i, ny) = ((phi(i, ny) - phi(i, ny-2)) / 4.0d0) + ((phi(i, ny) + phi(i, ny-2) - (2.0d0 * &
							& phi(i, ny-1))) / 8.0d0)
				del_n(i, ny) = ((phi(i-1, ny) - phi(i+1, ny)) / 4.0d0) + ((phi(i-1, ny) + phi(i+1, ny) - (2.0d0 * &
							& phi(i, ny))) / 8.0d0)
				del_s(i, ny) = ((phi(i, ny) - phi(i+2, ny)) / 4.0d0) + ((phi(i, ny) + phi(i+2, ny) - (2.0d0 * &
							& phi(i+1, ny))) / 8.0d0)
				b(i, ny) = (2.0d0 * dx * dy * (x_c(ny) + y_c(i))) + (Fw(i, ny) * del_w(i, ny)) - (Fn(i, ny) * &
						& del_n(i, ny)) + (Fs(i, ny) * del_s(i, ny))
						
			else 
				del_e(i, ny) = 0.0d0
				del_w(i, ny) = ((phi(i, ny) - phi(i, ny-2)) / 4.0d0) + ((phi(i, ny) + phi(i, ny-2) - (2.0d0 * &
							& phi(i, ny-1))) / 8.0d0)
				del_n(i, ny) = ((phi(i-1, ny) - phi(i+1, ny)) / 4.0d0) + ((phi(i-1, ny) + phi(i+1, ny) - (2.0d0 * &
							& phi(i, ny))) / 8.0d0)
				del_s(i, ny) = 0.0d0
				b(i, ny) = (2.0d0 * dx * dy * (x_c(ny) + y_c(i))) + (Fw(i, ny) * del_w(i, ny)) - (Fn(i, ny) * &
						& del_n(i, ny))
			end if
				
			!Calculating aP coefficient
			aP(i, ny) = aW(i, ny) + aN(i, ny) + aS(i, ny) + Fe(i, ny) - Fw(i, ny) + Fn(i, ny) - Fs(i, ny)
		end do
	
	end subroutine calc_east_boundary_coeff
	
	subroutine calc_boundary_coeff()
	
	implicit none
	
		call calc_south_boundary_coeff()
		call calc_west_boundary_coeff()
		call calc_north_boundary_coeff()
		call calc_east_boundary_coeff()
	
	end subroutine calc_boundary_coeff
	
	subroutine calc_bottom_left_cell()
	
	implicit none
	
		!Calculating diffusion coefficients				
		De(nx, 1) = diff_coeff
		Dw(nx, 1) = 2.0d0 * diff_coeff
		Dn(nx, 1) = diff_coeff
		Ds(nx, 1) = 2.0d0 * diff_coeff
		
		!Calculating convection coefficients
		Fe(nx, 1) = rho * dy * (1.0d0 + x_f(1)**2)
		Fw(nx, 1) = rho * dy * (1.0d0 + x_f(0)**2)
		Fn(nx, 1) = rho * dx * (1.0d0 + y_f(nx-1)**2)
		Fs(nx, 1) = rho * dx * (1.0d0 + y_f(nx)**2)
		
		!Calculating neighbouring coefficients
		aE(nx, 1) = De(nx, 1) 
		aW(nx, 1) = 0.0d0
		aN(nx, 1) = Dn(nx, 1) 
		aS(nx, 1) = 0.0d0
			
		!Calculating b term
		del_e(nx, 1) = 0.0d0
		del_w(nx, 1) = 0.0d0
		del_n(nx, 1) = 0.0d0
		del_s(nx, 1) = 0.0d0
		b(nx, 1) = (2.0d0 * dx * dy * (x_c(1) + y_c(nx))) + (Dw(nx, 1) * phi_b) + (Ds(nx, 1) * phi_a) + (Fw(nx, 1) * phi_b) + &
			    & (Fs(nx, 1) * phi_a)
				
		!Calculating aP coefficient
		aP(nx, 1) = aE(nx, 1) + aN(nx, 1) + Fe(nx, 1) + Fn(nx, 1) + Ds(nx, 1) + Dw(nx, 1)	
	
	end subroutine calc_bottom_left_cell
	
	subroutine calc_bottom_right_cell()
	
	implicit none
	
		!Calculating diffusion coefficients				
		De(nx, ny) = 0.0d0
		Dw(nx, ny) = diff_coeff
		Dn(nx, ny) = diff_coeff
		Ds(nx, ny) = 2.0d0 * diff_coeff
		
		!Calculating convection coefficients
		Fe(nx, ny) = rho * dy * (1.0d0 + x_f(ny)**2)
		Fw(nx, ny) = rho * dy * (1.0d0 + x_f(ny-1)**2)
		Fn(nx, ny) = rho * dx * (1.0d0 + y_f(nx-1)**2)
		Fs(nx, ny) = rho * dx * (1.0d0 + y_f(nx)**2)
		
		!Calculating neighbouring coefficients
		aE(nx, ny) = 0.0d0
		aW(nx, ny) = Dw(nx, ny) + Fw(nx, ny) 
		aN(nx, ny) = Dn(nx, ny)
		aS(nx, ny) = 0.0d0
			
		!Calculating b term
		del_e(nx, ny) = 0.0d0
		del_w(nx, ny) = ((phi(nx, ny) - phi(nx, ny-2)) / 4.0d0) + ((phi(nx, ny) + phi(nx, ny-2) - (2.0d0 * &
					& phi(nx, ny-1))) / 8.0d0)
		del_n(nx, ny) = 0.0d0
		del_s(nx, ny) = 0.0d0
		b(nx, ny) = (2.0d0 * dx * dy * (x_c(ny) + y_c(nx))) + (Ds(nx, ny) * phi_a) + (Fs(nx, ny) * phi_a) + &
				& (Fw(nx, ny) * del_w(nx, ny)) 
				
		!Calculating aP coefficient
		aP(nx, ny) = aW(nx, ny) + aN(nx, ny) + Fe(nx, ny) - Fw(nx, ny) + Fn(nx, ny) + Ds(nx, ny)	
	
	end subroutine calc_bottom_right_cell
	
	subroutine calc_top_left_cell()
	
	implicit none
	
		!Calculating diffusion coefficients				
		De(1, 1) = diff_coeff
		Dw(1, 1) = 2.0d0 * diff_coeff
		Dn(1, 1) = 0.0d0
		Ds(1, 1) = diff_coeff
		
		!Calculating convection coefficients
		Fe(1, 1) = rho * dy * (1.0d0 + x_f(1)**2)
		Fw(1, 1) = rho * dy * (1.0d0 + x_f(0)**2)
		Fn(1, 1) = rho * dx * (1.0d0 + y_f(0)**2)
		Fs(1, 1) = rho * dx * (1.0d0 + y_f(1)**2)
		
		!Calculating neighbouring coefficients
		aE(1, 1) = De(1, 1) 
		aW(1, 1) = 0.0d0
		aN(1, 1) = 0.0d0
		aS(1, 1) = Ds(1, 1) + Fs(1, 1)
			
		!Calculating b term
		del_e(1, 1) = 0.0d0
		del_w(1, 1) = 0.0d0
		del_n(1, 1) = 0.0d0
		del_s(1, 1) = ((phi(1, 1) - phi(1+2, 1)) / 4.0d0) + ((phi(1, 1) + phi(1+2, 1) - (2.0d0 * &
				& phi(1+1, 1))) / 8.0d0)
		b(1, 1) = (2.0d0 * dx * dy * (x_c(1) + y_c(1))) + (Dw(1, 1) * phi_b) + (Fw(1, 1) * phi_b) + &
				& (Fs(1, 1) * del_s(1, 1)) 
				
		!Calculating aP coefficient
		aP(1, 1) = aE(1, 1) + aS(1, 1) + Fe(1, 1) + Fn(1, 1) - Fs(1, 1) + Dw(1, 1)	
	
	end subroutine calc_top_left_cell
	
	subroutine calc_top_right_cell()
	
	implicit none
	
		!Calculating diffusion coefficients				
		De(1, ny) = 0.0d0
		Dw(1, ny) = diff_coeff
		Dn(1, ny) = 0.0d0
		Ds(1, ny) = diff_coeff
		
		!Calculating convection coefficients
		Fe(1, ny) = rho * dy * (1.0d0 + x_f(ny)**2)
		Fw(1, ny) = rho * dy * (1.0d0 + x_f(ny-1)**2)
		Fn(1, ny) = rho * dx * (1.0d0 + y_f(0)**2)
		Fs(1, ny) = rho * dx * (1.0d0 + y_f(1)**2)
		
		!Calculating neighbouring coefficients
		aE(1, ny) = 0.0d0
		aW(1, ny) = Dw(1, ny) + Fw(1, ny) 
		aN(1, ny) = 0.0d0
		aS(1, ny) = Ds(1, ny) + Fs(1, ny)
			
		!Calculating b term
		del_e(1, ny) = 0.0d0
		del_w(1, ny) = ((phi(1, ny) - phi(1, ny-2)) / 4.0d0) + ((phi(1, ny) + phi(1, ny-2) - (2.0d0 * &
					& phi(1, ny-1))) / 8.0d0)
		del_n(1, ny) = 0.0d0
		del_s(1, ny) = ((phi(1, ny) - phi(1+2, ny)) / 4.0d0) + ((phi(1, ny) + phi(1+2, ny) - (2.0d0 * &
				& phi(1+1, ny))) / 8.0d0)
		b(1, ny) = (2.0d0 * dx * dy * (x_c(ny) + y_c(1))) + (Fw(1, ny) * del_w(1, ny)) + (Fs(1, ny) * &
				& del_s(1, ny))  
				
		!Calculating aP coefficient
		aP(1, ny) = aW(1, ny) + aS(1, ny) + Fe(1, ny) - Fw(1, ny) + Fn(1, ny) - Fs(1, ny)	
	
	end subroutine calc_top_right_cell
	
	subroutine calc_corner_coeff()
	
	implicit none
	
		call calc_bottom_left_cell()
		call calc_bottom_right_cell()
		call calc_top_left_cell()
		call calc_top_right_cell()
		
	end subroutine calc_corner_coeff
	
	subroutine initialise_phi_field()
	
	implicit none
	
		phi(:, :) = phi_guess
		phi_old(:, :) = phi_guess
	
	end subroutine initialise_phi_field
	
	subroutine gauss_seidel()
	
	implicit none
	
		do i = 1, nx
			do j = 1, ny
				if ((i == 1) .and. (j == 1)) then
					phi(i, j) = ((aE(i, j) * phi(i, j+1)) + (aS(i, j) * phi(i+1, j)) + b(i, j)) / aP(i, j)
				end if
				
				if ((i == 1) .and. ((j > 1) .and. (j < ny))) then
					phi(i, j) = ((aE(i, j) * phi(i, j+1)) + (aW(i, j) * phi(i, j-1)) + (aS(i, j) * phi(i+1, j)) + &
							& b(i, j)) / aP(i, j)
				end if
				
				if ((i == 1) .and. (j == ny)) then
					phi(i, j) = ((aW(i, j) * phi(i, j-1)) + (aS(i, j) * phi(i+1, j)) + b(i, j)) / aP(i, j)
				end if
				
				if (((i > 1) .and. (i < nx)) .and. (j == 1)) then
					phi(i, j) = ((aE(i, j) * phi(i, j+1)) + (aN(i, j) * phi(i-1, j)) + (aS(i, j) * phi(i+1, j)) + &
							& b(i, j)) / aP(i, j)
				end if
				
				if (((i > 1) .and. (i < nx)) .and. (j == ny)) then
					phi(i, j) = ((aW(i, j) * phi(i, j-1)) + (aN(i, j) * phi(i-1, j)) + (aS(i, j) * phi(i+1, j)) + &
							& b(i, j)) / aP(i, j)
				end if
				
				if ((i == nx) .and. ((j > 1) .and. (j < ny))) then
					phi(i, j) = ((aE(i, j) * phi(i, j+1)) + (aW(i, j) * phi(i, j-1)) + (aN(i, j) * phi(i-1, j)) + &
							& b(i, j)) / aP(i, j)
				end if
				
				if ((i == nx) .and. (j == 1)) then
					phi(i, j) = ((aE(i, j) * phi(i, j+1)) + (aN(i, j) * phi(i-1, j)) + b(i, j)) / aP(i, j)
				end if
				
				if ((i == nx) .and. (j == ny)) then
					phi(i, j) = ((aW(i, j) * phi(i, j-1)) + (aN(i, j) * phi(i-1, j)) + b(i, j)) / aP(i, j)
				end if
				
				if (((i > 1) .and. (i < nx)) .and. ((j > 1) .and. (j < ny))) then
					phi(i, j) = ((aE(i, j) * phi(i, j+1)) + (aW(i, j) * phi(i, j-1)) + (aN(i, j) * phi(i-1, j)) + &
					    		& (aS(i, j) * phi(i+1, j)) + b(i, j)) / aP(i, j)
				end if
			end do
		end do
	
	end subroutine gauss_seidel
	
	subroutine calc_l2norm()
	
	implicit none
	
		numerator = 0.0d0
		denominator = 0.0d0
		do i = 1, nx
			do j = 1, ny
				numerator = numerator + (phi(i, j) - phi_old(i, j))**2
				denominator = denominator + phi(i, j)**2
			end do
		end do
		phil2norm = sqrt(numerator / denominator)
			
	end subroutine calc_l2norm
	
	subroutine solver()
	
	implicit none
	
		do while (phil2norm > tolerance)
			call gauss_seidel()
			
			!Apply relaxation
			phi(:, :) = (alpha * phi(:, :)) + ((1.0d0 - alpha) * phi_old(:, :))
			
			call calc_l2norm()
			phi_old(:, :) = phi(:, :)
			iter = iter + 1
		end do
	
	end subroutine solver
	
	subroutine post_process_data()
	
	implicit none
	
		print*, ''
		print*, 'The solution is converged to:'
		print*, ''
		
		do i = 1, nx
			write(*, 1) (phi(i, j), j = 1, ny)
			1 format(50f8.4)
		end do
		
		!Creating file name
		write(nxx, '(I5)') nx
		write(nyy, '(I5)') ny
		nxx = trim(adjustl(nxx))
		nyy = trim(adjustl(nyy))
		fname1 = 'quick_converged_solution_'//trim(nxx)//'_'//trim(nyy)//'.txt'
		fname2 = 'quick_horizontal_center_line_data_'//trim(nxx)//'_'//trim(nyy)//'.txt'
		fname3 = 'quick_solver_runtime_details_'//trim(nxx)//'_'//trim(nyy)//'.txt'
		
		open(11, file = fname1)
		open(12, file = fname2)
		open(13, file = fname3)
		
		!Writing solutions in the file
		do j = 1, ny
			do i = nx, 1, -1
				write(11, 2) x_c(j) / L, y_c(i) / H, phi(i, j)
				2 format(50f15.7)
			end do
			write(11, *) ''
		end do
		
		!Getting x data at horizontal center line
		x_center(0) = 0.0d0 / L
		x_center(ny+1) = L / L
		do i = 1, ny
			x_center(i) = x_c(i) / L
		end do
		
		!Getting phi data at horizontal center line
		phi_center(0) = 0.0d0
		phi_center(ny+1) = (phi(nx/2, ny) + phi((nx/2)+1, ny)) / 2.0d0
		do i = 1, ny
			phi_center(i) = (phi(nx/2, i) + phi((nx/2) + 1, i)) / 2.0d0
		end do
		
		!Writing horizontal center line data to the file
		do i = 0, ny+1
			write(12, 1) x_center(i), phi_center(i)
		end do		
	
	end subroutine post_process_data 

end module functions

!Main Program
program quick

use functions

implicit none

	real:: start_time, end_time
	
	call cpu_time(start_time)
	call allocate_memory()
	call set_parameters()
	call initialise_solver()
	call calc_interior_coeff()
	call calc_boundary_coeff()
	call calc_corner_coeff()
	call initialise_phi_field()
	call solver()
	call post_process_data()
	call deallocate_memory()
	call cpu_time(end_time)
	
	print*, ''
	print*, 'Program iterations = ', iter
	print*, 'Program running time = ', end_time-start_time, 'seconds' 
	write(13, *) 'Program iterations = ', iter
	write(13, *) 'Program running time = ', end_time-start_time, 'seconds'

end program quick
