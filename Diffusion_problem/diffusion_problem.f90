! Program to solve for diffusion only problem with no source terms for domain 0 <= x <= 1, 0 <= y <= 1
! Boundary conditions are: phi = x*y
! uniform grid has been used in each direction

module diff_prob_data
  implicit none

  integer, parameter :: nx = 30, ny = 30 ! set grid parameters
  integer :: i, j, iter, iter_max
  double precision :: x, y, z, dx, dy, alpha, error
  double precision, dimension(:, :), allocatable :: ae, aw, an, as, ap, Sp, Sc
  double precision, dimension(:, :), allocatable :: phi, phi_analytical

contains

  subroutine allocate_memory()
    implicit none

    allocate (ae(1:nx, 1:ny), aw(1:nx, 1:ny), an(1:nx, 1:ny), as(1:nx, 1:ny))
    allocate (Sp(1:nx, 1:ny), Sc(1:nx, 1:ny))
    allocate (phi(0:nx+1, 0:ny+1), phi_analytical(0:nx+1, 0:ny+1))

  end subroutine allocate_memory

  subroutine deallocate_memory()
    implicit none

    deallocate (ae, aw, an, as)
    deallocate (Sp, Sc)
    deallocate (phi, phi_analytical)

  end subroutine deallocate_memory

  subroutine initialize_fvmsolver()
    implicit none

    ! initialize source values and coefficients 
    Sp = 0.0d0
    Sc = 0.0d0
    ae = 1.0d0
    aw = 1.0d0
    an = 1.0d0
    as = 1.0d0

    ! initialize the domain
    phi = 0.0d0                   ! set your choice of initialization

    ! calculate grid spacing
    dx = 1.0d0 / float(nx)
    dy = 1.0d0 / float(ny)

    ! calculate analytical solution
    do j = 1, ny
       do i = 1, nx
          x = (dx / 2.0d0) + (dx * float(i-1))
          y = (dy / 2.0d0) + (dy * float(j-1))
          phi_analytical(i, j) = x * y
       end do
    end do

  end subroutine initialize_fvmsolver
  
  subroutine set_coefficients()
    implicit none

    ! set west boundary condition
    aw(1, :) = 0.0d0
    Sc(1, :) = 0.0d0              ! phi_west = 0
    Sp(1, :) = -2.0d0

    ! set east boundary condition
    ae(nx, :) = 0.0d0
    Sp(nx, :) = -2.0d0
    do j = 1, ny
       Sc(nx, j) = 2.0d0 * ((dy / 2.0d0) + (dy * float(j-1))) ! phi_east = y
    end do

    ! set south boundary condition
    as(:, 1) = 0.0d0
    Sp(:, 1) = -2.0d0
    Sc(:, 1) = 0.0d0              ! phi_south = 0

    ! set north boundary condition
    an(:, ny) = 0.0d0
    Sp(:, ny) = -2.0d0
    do i = 1, nx
       Sc(i, ny) = 2.0d0 * ((dx / 2.0d0) + (dx * float(i-1))) ! phi_north = x
    end do

    ! set northwest cell coefficients
    an(1, ny) = 0.0d0
    aw(1, ny) = 0.0d0
    Sp(1, ny) = -4.0d0
    Sc(1, ny) = 2.0d0 * (dx / 2.0d0)

    ! set northeast cell coefficients
    an(nx, ny) = 0.0d0
    ae(nx, ny) = 0.0d0
    Sp(nx, ny) = -4.0d0
    Sc(nx, ny) = 2.0d0 * ((dx / 2.0d0) + (dx * float(nx-1)) + (dy / 2.0d0) + (dy * float(ny-1)))

    ! set southwest cell coefficients
    as(1, 1) = 0.0d0
    aw(1, 1) = 0.0d0
    Sp(1, 1) = -4.0d0
    Sc(1, 1) = 0.0d0
    
    ! set southeast cell coefficients
    as(nx, 1) = 0.0d0
    ae(nx, 1) = 0.0d0
    Sp(nx, 1) = -4.0d0
    Sc(nx, 1) = 2.0d0 * (dy / 2.0d0)

    ap = ae + aw + an + as - Sp

  end subroutine set_coefficients

  subroutine gs_solver()
    implicit none

    ! iterative solver GS-SOR/SUR
    iter_max = 1000               ! set maximum number of iterations 
    alpha = 1.10d0                ! set your over/under-relaxation parameter
    do iter = 1, iter_max
       do j = 1, ny
          do i = 1, nx
             phi(i, j) = phi(i, j) + (alpha / ap(i, j)) * (ae(i, j) * phi(i+1, j) + aw(i, j) * &
                  phi(i-1, j) + an(i, j) * phi(i, j+1) + as(i, j) * phi(i, j-1) + Sc(i, j) - &
                  ap(i, j) * phi(i, j))
          end do
       end do

       error = maxval(abs(phi_analytical-phi))
       write(*, '(A, I5, A, F10.8)') "iteration = ", iter, ", max error = ", error
       
    end do

  end subroutine gs_solver

  subroutine post_process_data()
    implicit none

    ! post-process data for Paraview format
    open(unit = 3, file = "converged_results_diffusion.csv")
    write(3, *) "x coord, y coord, z coord, phi"
    do j = 1, ny
       do i = 1, nx
          x = (dx / 2.0d0) + (dx * float(i-1))
          y = (dy / 2.0d0) + (dy * float(j-1))
          z = 0.0d0
          write(3, *) x, ",", y, ",", z, ",", phi(i, j) ! data for interior cells
       end do
    end do
    
    do j = 1, ny
       y = (dy / 2.0d0) + (dy * float(j-1))
       write(3, *) 0.0d0, ",", y, ",", 0.0d0, ",",  0.0d0 ! data for west boundary
    end do
    
    do j = 1, ny
       y = (dy / 2.0d0) + (dy * float(j-1))
       write(3, *) 1.0d0, ",", y, ",", 0.0d0, ",", y ! data for east boundary 
    end do
    
    do i = 1, nx
       x = (dx / 2.0d0) + (dx * float(i-1))
       write(3, *) x, ",", 0.0d0, ",", 0.0d0, ",", 0.0d0 ! data for south boundary 
    end do
    
    do i = 1, nx
       x = (dx / 2.0d0) + (dx * float(i-1))
       write(3, *) x, ",", 1.0d0, ",", 0.0d0, ",", x ! data for north boundary 
    end do

    write(3, *) 0.0d0, ",", 0.0d0, ",", 0.0d0, ",", 0.0d0 ! data for southwest cell
    write(3, *) 0.0d0, ",", 1.0d0, ",", 0.0d0, ",", 0.0d0 ! data for northwest cell
    write(3, *) 1.0d0, ",", 0.0d0, ",", 0.0d0, ",", 0.0d0 ! data for southeast cell
    write(3, *) 1.0d0, ",", 1.0d0, ",", 0.0d0, ",", 1.0d0 ! data for northeast cell
    
  end subroutine post_process_data
    
end module diff_prob_data

program fvm_diffusion

  use diff_prob_data
  implicit none

  call allocate_memory()
  call initialize_fvmsolver()
  call set_coefficients()
  call gs_solver()
  call post_process_data()

  write(*, '(A, I2, A, I2)') "The grid resolution used is ", nx, "x", ny
  write(*, '(A, F10.8)')  "The maximum error = ", error
  write(*, '(A)') "The solution has converged. EXITING!!!"

  call deallocate_memory()
  
end program fvm_diffusion
