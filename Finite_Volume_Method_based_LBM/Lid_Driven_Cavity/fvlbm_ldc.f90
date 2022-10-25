!Code for Finite Volume LBM lid-driven cavity case using structured Cartesian grid using D2Q9 lattice type

!Module containing required subroutines
module functions

implicit none

	integer, parameter:: n = 128, m = 128, q = 9				!Domain with Nx = n, Ny = m cells
	double precision:: L = 1.0d0, H = 1.0d0				!L = length of domain, H = height of domain		
	double precision, dimension(:, :, :), allocatable:: f_node, feq_node, RHS	
	double precision, dimension(:, :, :), allocatable:: fn_face, fe_face, fw_face, fs_face
	double precision, dimension(:, :, :), allocatable:: dfdx, dfdy
	double precision, dimension(:, :), allocatable:: u, v, rho, u_prev, v_prev
	double precision, dimension(:, :), allocatable:: f_gc_face_top, f_gc_face_bottom, f_gc_face_right, f_gc_face_left
	double precision, dimension(:), allocatable:: w, ex, ey
	double precision, dimension(:), allocatable:: x_f, y_f, x_c, y_c, u_numerical, y_vertical, x_horizontal, v_numerical
	double precision:: t1, t2, num, denom
	double precision:: tau, nu, dx, dy, dt, tolerance, l2norm, Uplate, Re
	integer:: i, j, k, steps, tt
	character (len = 256):: nx, ny
	character (len = 1024):: fname1, fname2, fname3
	
contains
	
	subroutine allocate_memory()
	
	implicit none
	
		allocate (f_node(1:q, 1:n, 1:m), feq_node(1:q, 1:n, 1:m), RHS(1:q, 1:n, 1:m))
		allocate (fn_face(1:q, 1:n, 1:m), fe_face(1:q, 1:n, 1:m), fw_face(1:q, 1:n, 1:m), fs_face(1:q, 1:n, 1:m))
		allocate (dfdx(1:q, 1:n, 1:m), dfdy(1:q, 1:n, 1:m))
		allocate (u(1:n, 1:m), v(1:n, 1:m), rho(1:n, 1:m), u_prev(1:n, 1:m), v_prev(1:n, 1:m))
		allocate (f_gc_face_top(1:q, 1:n), f_gc_face_bottom(1:q, 1:n), f_gc_face_right(1:q, 1:m), f_gc_face_left(1:q, 1:m))
		allocate (x_f(0:n), y_f(0:m), x_c(1:n), y_c(1:m), u_numerical(0:m+1), y_vertical(0:m+1), x_horizontal(0:n+1), v_numerical(0:n+1))
		allocate (w(1:q), ex(1:q), ey(1:q))
	
	end subroutine allocate_memory
	
	subroutine deallocate_memory()
	
	implicit none
	
		deallocate (f_node, feq_node, RHS)
		deallocate (fn_face, fe_face, fw_face, fs_face)
		deallocate (dfdx, dfdy)
		deallocate (u, v, rho, u_prev, v_prev)
		deallocate (f_gc_face_top, f_gc_face_bottom, f_gc_face_right, f_gc_face_left)
		deallocate (x_f, y_f, x_c, y_c, u_numerical, y_vertical, x_horizontal, v_numerical)
		deallocate (w, ex, ey)
		
	end subroutine deallocate_memory

	subroutine fill_wk()
	
	implicit none
	
		w(:) = (/ 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 /)
		
	end subroutine fill_wk
	
	subroutine fill_ek()
	
	implicit none
	
		ex(:) = (/ 0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0, -1.0d0, -1.0d0, 1.0d0 /)
		ey(:) = (/ 0.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0, 1.0d0, 1.0d0, -1.0d0, -1.0d0 /)
		
	end subroutine fill_ek
	
	subroutine set_parameters()
	
	implicit none
	
		u(:, :) = 0.0d0				!Set initial guess value for u-velocity
		v(:, :) = 0.0d0				!Set initial guess value for v-velocity
		u_prev(:, :) = u(:, :)
		v_prev(:, :) = v(:, :)
		rho(:, :) = 1.0d0				!Set initial guess value for density
		dx = L / n
		dy = H / m
		steps = 1
		tolerance = 1.0d-8				!Set tolerance value 
		dt = 5.0d-5
		Uplate = 1.0d-1				!Set plate velocity
		Re = 400
		nu = Uplate * L / Re
		tau = 3.0d0 * nu
		
		print*, 'Re = ', Re
		print*, 'viscosity = ', nu
		
	end subroutine set_parameters
	
	subroutine calculate_faces()
	
	implicit none
	
		x_f(0) = 0.0d0
		y_f(0) = 0.0d0
		
		do i = 1, n
			x_f(i) = x_f(i-1) + dx
		end do
		
		do j = 1, m
			y_f(j) = y_f(j-1) + dy
		end do
	
	end subroutine calculate_faces
	
	subroutine calculate_centers()
	
	implicit none
	
		x_c(1) = x_f(0) + (5.0d-1 * dx)
		y_c(1) = y_f(0) + (5.0d-1 * dy)
		
		!Populate x at centers
		do i = 2, n
			x_c(i) = x_c(i-1) + dx
		end do
		
		!Populate y at centers
		do j = 2, m
			y_c(j) = y_c(j-1) + dy
		end do
	
	end subroutine calculate_centers
	
	subroutine initialise_f()
	
	implicit none
		
		do i = 1, n
			do j = 1, m
				t2 = u(i, j) * u(i, j) + v(i, j) * v(i, j)
				do k = 1, q
					t1 = u(i, j) * ex(k) + v(i, j) * ey(k)
					feq_node(k, i, j) = rho(i, j) * w(k) * (1.0d0 + (3.0d0 * t1) + (4.50d0 * t1 * t1) - (1.50d0 * t2))
					f_node(k, i, j) = feq_node(k, i, j)
				end do
			end do
		end do
		
	end subroutine initialise_f
	
	subroutine initialise_fvlbm()
	
	implicit none
		
		call fill_wk()
		call fill_ek()
		call set_parameters()
		call calculate_faces()
		call calculate_centers()
		call initialise_f()
		
	end subroutine initialise_fvlbm
	
	subroutine update_ghost_faces()
	
	implicit none
	
		!Top boundary -- moving plate, Bottom boundary -- stationary wall
		do i = 1, n
			do k = 1, q
				f_gc_face_top(k, i) = f_node(k, i, m) - feq_node(k, i, m) + (rho(i, m) * w(k) * (1.0d0 + &
				& (3.0d0 * ex(k) * Uplate) + (4.50d0 * (ex(k) * Uplate) ** 2) - &
				& (1.50d0 * Uplate * Uplate)))
				
				f_gc_face_bottom(k, i) = f_node(k, i, 1) - feq_node(k, i, 1) + (rho(i, 1) * w(k) * 1.0d0)
			end do
		end do
		
		!Right wall -- stationary wall, Left wall -- stationary wall
		do j = 1, m
			do k = 1, q
				f_gc_face_right(k, j) = f_node(k, n, j) - feq_node(k, n, j) + (rho(n, j) * w(k) * 1.0d0)
				f_gc_face_left(k, j) = f_node(k, 1, j) - feq_node(k, 1, j) + (rho(1, j) * w(k) * 1.0d0)
			end do
		end do
		
	end subroutine update_ghost_faces
	
	subroutine calc_feq()
	
	implicit none
	
		do i = 1, n
			do j = 1, m
				t2 = u(i, j) * u(i, j) + v(i, j) * v(i, j)
				do k = 1, q
					t1 = u(i, j) * ex(k) + v(i, j) * ey(k)
					feq_node(k, i, j) = rho(i, j) * w(k) * (1.0d0 + (3.0d0 * t1) + (4.50d0 * t1 * t1) - (1.50d0 * t2))
				end do
			end do
		end do
		
	end subroutine calc_feq
	
	subroutine calc_rhs()
	
	implicit none
	
		do i = 1, n
			do j = 1, m
				do k = 1, q
					RHS(k, i, j) = -((1.0d0 / tau) * (f_node(k, i, j) - feq_node(k, i, j))) - &
						&	(((ex(k) * fe_face(k, i, j) * dy) - (ex(k) * fw_face(k, i, j) * dy) + &
						&	(ey(k) * fn_face(k, i, j) * dx) - (ey(k) * fs_face(k, i, j) * dx)) / (dx * dy))
				end do
			end do
		end do
	
	end subroutine calc_rhs
	
	subroutine calc_macro_variables()
	
	implicit none
	
		do i = 1, n
			do j = 1, m
				t1 = 0.0d0
				do k = 1, q
					t1 = t1 + f_node(k, i, j)
				end do
				rho(i, j) = t1
			end do
		end do
			
		do i = 1, n
			do j = 1, m
				t1 = 0.0d0
				t2 = 0.0d0
				do k = 1, q
					t1 = t1 + f_node(k, i, j) * ex(k)
					t2 = t2 + f_node(k, i, j) * ey(k)
				end do
				u(i, j) = t1 / rho(i, j)
				v(i, j) = t2 / rho(i, j)
			end do
		end do
	
	end subroutine calc_macro_variables
	
	subroutine calc_norm()
	
	implicit none
	
		num = 0.0d0
		denom = 0.0d0
			
		do i = 1, n
			do j = 1, m
				num = num + ((u(i, j) - u_prev(i, j))**2 + (v(i, j) - v_prev(i, j))**2)
				denom = denom + (u(i, j)**2 + v(i, j)**2)
			end do
		end do
		
		l2norm = sqrt(num / denom)
	
	end subroutine calc_norm
	
	subroutine calc_dfdx()
	
	implicit none
		
		do i = 2, n-1
			do j = 1, m
				do k = 1, q
					dfdx(k, i, j) = (f_node(k, i+1, j) - f_node(k, i-1, j)) / (2.0d0 * dx)
				end do
			end do
		end do
		
		do j = 1, m
			do k = 1, q
!				dfdx(k, 1, j) = ((4.0d0 * f_node(k, 2, j)) - f_node(k, 3, j) - (3.0d0 * f_node(k, 1, j))) / (2.0d0 * dx)
				dfdx(k, 1, j) = (f_node(k, 2, j) - f_node(k, 1, j)) / dx
!				dfdx(k, n, j) = ((3.0d0 * f_node(k, n, j)) - (4.0d0 * f_node(k, n-1, j)) + f_node(k, n-2, j)) / (2.0d0 * dx)
				dfdx(k, n, j) = (f_node(k, n, j) - f_node(k, n-1, j)) / dx
			end do
		end do			
	
	end subroutine calc_dfdx
	
	subroutine calc_dfdy()
	
	implicit none
		
		do i = 1, n
			do j = 2, m-1
				do k = 1, q
					dfdy(k, i, j) = (f_node(k, i, j+1) - f_node(k, i, j-1)) / (2.0d0 * dy)
				end do
			end do
		end do
		
		do i = 1, n
			do k = 1, q
!				dfdy(k, i, 1) = ((4.0d0 * f_node(k, i, 2)) - f_node(k, i, 3) - (3.0d0 * f_node(k, i, 1))) / (2.0d0 * dy)
				dfdy(k, i, 1) = (f_node(k, i, 2) - f_node(k, i, 1)) / dy
!				dfdy(k, i, m) = ((3.0d0 * f_node(k, i, m)) - (4.0d0 * f_node(k, i, m-1)) + f_node(k, i, m-2)) / (2.0d0 * dy)
				dfdy(k, i, m) = (f_node(k, i, m) - f_node(k, i, m-1)) / dy
			end do
		end do			
	
	end subroutine calc_dfdy
	
	subroutine calculate_f_faces_with_grad()			!Using 2nd  order UDS
	
	implicit none
	
		call calc_dfdx()
		call calc_dfdy()
	
		do i = 2, n-1
			do j = 2, m-1
				do k = 1, q
					if (ex(k) .ge. 0) then
						fe_face(k, i, j) = f_node(k, i, j) + (dfdx(k, i, j) * dx / 2.0d0)
						fw_face(k, i, j) = f_node(k, i-1, j) + (dfdx(k, i-1, j) * dx / 2.0d0)
					else
						fe_face(k, i, j) = f_node(k, i+1, j) - (dfdx(k, i+1, j) * dx / 2.0d0)
						fw_face(k, i, j) = f_node(k, i, j) - (dfdx(k, i, j) * dx / 2.0d0)
					endif
					
					if (ey(k) .ge. 0) then
						fn_face(k, i, j) = f_node(k, i, j) + (dfdy(k, i, j) * dy / 2.0d0)
						fs_face(k, i, j) = f_node(k, i, j-1) + (dfdy(k, i, j-1) * dy / 2.0d0)
					else
						fn_face(k, i, j) = f_node(k, i, j+1) - (dfdy(k, i, j+1) * dy / 2.0d0)
						fs_face(k, i, j) = f_node(k, i, j) - (dfdy(k, i, j) * dy / 2.0d0)
					endif
				end do
			end do
		end do
		
		!Corner cell treatment
		do k = 1, q
			!i = 1, j = 1
			fw_face(k, 1, 1) = f_gc_face_left(k, 1)
			fs_face(k, 1, 1) = f_gc_face_bottom(k, 1)
			
			if (ex(k) .ge. 0) then
				fe_face(k, 1, 1) = f_node(k, 1, 1) + 	(dfdx(k, 1, 1) * dx / 2.0d0)
			else
				fe_face(k, 1, 1) = f_node(k, 2, 1) - (dfdx(k, 2, 1) * dx / 2.0d0)
			endif
			
			if (ey(k) .ge. 0) then
				fn_face(k, 1, 1) = f_node(k, 1, 1) + (dfdy(k, 1, 1) * dy / 2.0d0)
			else
				fn_face(k, 1, 1) = f_node(k, 1, 2) - (dfdy(k, 1, 2) * dy / 2.0d0)
			endif
			
			!i = n, j = 1
			fe_face(k, n, 1) = f_gc_face_right(k, 1)
			fs_face(k, n, 1) = f_gc_face_bottom(k, n)
			
			if (ex(k) .ge. 0) then
				fw_face(k, n, 1) = f_node(k, n-1, 1) + (dfdx(k, n-1, 1) * dx / 2.0d0)
			else
				fw_face(k, n, 1) = f_node(k, n, 1) - (dfdx(k, n, 1) * dx / 2.0d0)
			endif
			
			if (ey(k) .ge. 0) then
				fn_face(k, n, 1) = f_node(k, n, 1) + (dfdy(k, n, 1) * dy / 2.0d0)
			else
				fn_face(k, n, 1) = f_node(k, n, 2) - (dfdy(k, n, 2) * dy / 2.0d0)
			endif
			
			!i = 1, j = m
			fw_face(k, 1, m) = f_gc_face_left(k, m)
			fn_face(k, 1, m) = f_gc_face_top(k, 1)
			
			if (ex(k) .ge. 0) then
				fe_face(k, 1, m) = f_node(k, 1, m) + (dfdx(k, 1, m) * dx / 2.0d0)
			else
				fe_face(k, 1, m) = f_node(k, 2, m) - (dfdx(k, 2, m) * dx / 2.0d0)
			endif
			
			if (ey(k) .ge. 0) then
				fs_face(k, 1, m) = f_node(k, 1, m-1) + (dfdy(k, 1, m-1) * dy / 2.0d0)
			else
				fs_face(k, 1, m) = f_node(k, 1, m) - (dfdy(k, 1, m) * dy / 2.0d0)
			endif
			
			!i = n, j = m
			fe_face(k, n, m) = f_gc_face_right(k, m) 
			fn_face(k, n, m) = f_gc_face_top(k, n)
			
			if (ex(k) .ge. 0) then
				fw_face(k, n, m) = f_node(k, n-1, m) + (dfdx(k, n-1, m) * dx / 2.0d0)
			else
				fw_face(k, n, m) = f_node(k, n, m) - (dfdx(k, n, m) * dx / 2.0d0)
			endif
			
			if (ey(k) .ge. 0) then
				fs_face(k, n, m) = f_node(k, n, m-1) + (dfdy(k, n, m-1) * dy / 2.0d0)
			else
				fs_face(k, n, m) = f_node(k, n, m) - (dfdy(k, n, m) * dy / 2.0d0)
			endif
		end do

		!Boundary cell treatment
		do i = 2, n-1
			do k = 1, q
				fs_face(k, i, 1) = f_gc_face_bottom(k, i)
				fn_face(k, i, m) = f_gc_face_top(k, i)
				
				if (ex(k) .ge. 0) then
					fe_face(k, i, 1) = f_node(k, i, 1) + (dfdx(k, i, 1) * dx / 2.0d0)
					fe_face(k, i, m) = f_node(k, i, m) + (dfdx(k, i, m) * dx / 2.0d0)
					fw_face(k, i, 1) = f_node(k, i-1, 1) + (dfdx(k, i-1, 1) * dx / 2.0d0)
					fw_face(k, i, m) = f_node(k, i-1, m) + (dfdx(k, i-1, m) * dx / 2.0d0)
				else
					fe_face(k, i, 1) = f_node(k, i+1, 1) - (dfdx(k, i+1, 1) * dx / 2.0d0)
					fe_face(k, i, m) = f_node(k, i+1, m) - (dfdx(k, i+1, m) * dx / 2.0d0)
					fw_face(k, i, 1) = f_node(k, i, 1) - (dfdx(k, i, 1) * dx / 2.0d0)
					fw_face(k, i, m) = f_node(k, i, m) - (dfdx(k, i, m) * dx / 2.0d0)
				endif
						
				if (ey(k) .ge. 0) then
					fn_face(k, i, 1) = f_node(k, i, 1) + (dfdy(k, i, 1) * dy / 2.0d0)
					fs_face(k, i, m) = f_node(k, i, m-1) + (dfdy(k, i, m-1) * dy / 2.0d0)
				else
					fn_face(k, i, 1) = f_node(k, i, 2) - (dfdy(k, i, 2) * dy / 2.0d0)
					fs_face(k, i, m) = f_node(k, i, m) - (dfdy(k, i, m) * dy / 2.0d0)
				endif
			end do
		end do
		
		do j = 2, m-1
			do k = 1, q
				fw_face(k, 1, j) = f_gc_face_left(k, j)
				fe_face(k, n, j) = f_gc_face_right(k, j)
				
				if (ex(k) .ge. 0) then
					fe_face(k, 1, j) = f_node(k, 1, j) + (dfdx(k, 1, j) * dx / 2.0d0)
					fw_face(k, n, j) = f_node(k, n-1, j) + (dfdx(k, n-1, j) * dx / 2.0d0)
				else
					fe_face(k, 1, j) = f_node(k, 2, j) - (dfdx(k, 2, j) * dx / 2.0d0)
					fw_face(k, n, j) = f_node(k, n, j) - (dfdx(k, n, j) * dx / 2.0d0)
				endif
				
				if (ey(k) .ge. 0) then
					fn_face(k, 1, j) = f_node(k, 1, j) + (dfdy(k, 1, j) * dy / 2.0d0)
					fn_face(k, n, j) = f_node(k, n, j) + (dfdy(k, n, j) * dy / 2.0d0)
					fs_face(k, 1, j) = f_node(k, 1, j-1) + (dfdy(k, 1, j-1) * dy / 2.0d0)
					fs_face(k, n, j) = f_node(k, n, j-1) + (dfdy(k, n, j-1) * dy / 2.0d0)
				else
					fn_face(k, 1, j) = f_node(k, 1, j+1) - (dfdy(k, 1, j+1) * dy / 2.0d0)
					fn_face(k, n, j) = f_node(k, n, j+1) - (dfdy(k, n, j+1) * dy / 2.0d0)
					fs_face(k, 1, j) = f_node(k, 1, j) - (dfdy(k, 1, j) * dy / 2.0d0)
					fs_face(k, n, j) = f_node(k, n, j) - (dfdy(k, n, j) * dy / 2.0d0)
				endif
			end do
		end do		
						
	end subroutine calculate_f_faces_with_grad
	
	subroutine solve_fvlbm_euler_explicit()
	
	implicit none
	
		call initialise_fvlbm()
	
		l2norm = 1.0d0
		
		print*, ''
		print*, 'Time          Error          u_velocity at different locations'
		print*, ''
		
		do while (l2norm .gt. tolerance .or. steps .lt. 100)	
!		do tt = 1, 1000000
			
			call calc_feq()
			call update_ghost_faces()
			call calculate_f_faces_with_grad()			!using 2nd order UDS
			call calc_rhs()
			
			!Temporal scheme implementation
			do i = 1, n
				do j = 1, m
					do k = 1, q
						f_node(k, i, j) = f_node(k, i, j) + (dt * RHS(k, i, j))
					end do
				end do
			end do
			
			call calc_macro_variables()
			call calc_norm()
			
			u_prev(:, :) = u(:, :)
			v_prev(:, :) = v(:, :)
			
			steps = steps + 1
			
			!Print on terminal
			if (mod(steps, 1000) .eq. 0) then
				write(*, 1) steps * dt, l2norm, u(n/2, 1)/Uplate, u(n/2, m/4)/Uplate, u(n/2, m/2)/Uplate, u(n/2, 3*m/4)/Uplate, u(n/2, m)/Uplate
				1 format(7f15.10)
			endif
		end do			
	
	end subroutine solve_fvlbm_euler_explicit
	
	subroutine post_process()
	
	implicit none
	
		!Creating file name
		write(nx, '(I5)') n
		write(ny, '(I5)') m
		nx = trim(adjustl(nx))
		ny = trim(adjustl(ny))
		fname1 = 'uv_field_c2_'//trim(nx)//'_'//trim(ny)//'.txt'
		fname2 = 'u_prof_c2_'//trim(nx)//'_'//trim(ny)//'.txt'
		fname3 = 'v_prof_c2_'//trim(nx)//'_'//trim(ny)//'.txt'
		
		open(11, file = fname1)
		open(12, file = fname2)
		open(13, file = fname3)
		
		!Writing solutions in the file
		do i = 1, n
			do j = 1, m
				write(11, *) x_c(i), y_c(j), u(i, j), v(i, j), rho(i, j)
			end do
		end do 
		
		!Calculate x locations on horizontal line
		x_horizontal(0) = 0.0d0
		x_horizontal(n+1) = L / L
		do i = 1, n
			x_horizontal(i) = x_c(i) / L
		end do
		
		!Calculate y locations on vertical line
		y_vertical(0) = 0.0d0
		y_vertical(m+1) = H / H
		do j = 1, m
			y_vertical(j) = y_c(j) / H
		end do
		
		!Extracting numerical data at mid vertical line
		u_numerical(0) = 0.0d0
		u_numerical(m+1) = Uplate
		do j = 1, m
			u_numerical(j) = (u(n/2, j) + u((n/2) + 1, j)) / 2.0d0
		end do
		
		!Extracting numerical data at mid horizontal line
		v_numerical(0) = 0.0d0
		v_numerical(n+1) = 0.0d0
		do i = 1, n
			v_numerical(i) = (v(i, m/2) + v(i, (m/2) + 1)) / 2.0d0
		end do
		
		!Writing data to file
		do j = 0, m+1
			write(12, *) y_vertical(j), u_numerical(j) / Uplate
		end do
		
		!Writing data to file
		do i = 0, n+1
			write(13, *) x_horizontal(i), v_numerical(i) / Uplate
		end do

	end subroutine post_process 	
	
end module functions

!Main Program
program FVLBM

use functions

implicit none

	real:: start_time, end_time
		
	call cpu_time(start_time)
	call allocate_memory()
	call solve_fvlbm_euler_explicit()
	call post_process()
	call deallocate_memory()
	call cpu_time(end_time)
	
	print*, 'Program running time = ', end_time - start_time, 'seconds'

end program FVLBM
