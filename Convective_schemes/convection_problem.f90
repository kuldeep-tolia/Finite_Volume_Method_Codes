! Program to solve for convection only problem with no source terms for domain 0 <= x <= 1, 0 <= y <= 1
! Boundary conditions are: phi_west = 1.0, phi_south = 0.0, zero-gradient at the outflow
! The domain has a velocity of 1.0 m/s
! uniform grid has been used in each direction
! a hybrid convective scheme is used (blending between 1st order upwind and central difference scheme)
! beta = 0.0 => 1st order upwind
! beta = 1.0 => CDS

module conv_prob_data
  implicit none

  integer, parameter :: nx = 80, ny = 80 ! set grid parameters
  integer :: i, j, iter, iter_max
  double precision, dimension(:), allocatable :: x, y
  double precision, dimension(:, :), allocatable :: phi, ap_upw, ae_upw, aw_upw, an_upw, as_upw, Sp_upw, Sc_upw
  double precision, dimension(:, :), allocatable :: ap_cds, ae_cds, aw_cds, an_cds, as_cds, Sp_cds, Sc_cds, phi_old
  double precision :: z, dx, dy, blend_factor, rho, phi_west, phi_south, domain_velocity
  double precision :: Fe, Fw, Fn, Fs

contains

  subroutine allocate_memory()
    implicit none

    allocate (x(0:nx+1), y(0:ny+1))
    allocate (phi(0:nx+1, 0:ny+1), ap_upw(0:nx+1, 0:ny+1), ae_upw(0:nx+1, 0:ny+1), aw_upw(0:nx+1, 0:ny+1))
    allocate (an_upw(0:nx+1, 0:ny+1), as_upw(0:nx+1, 0:ny+1), Sp_upw(0:nx+1, 0:ny+1), Sc_upw(0:nx+1, 0:ny+1))
    allocate (ap_cds(0:nx+1, 0:ny+1), ae_cds(0:nx+1, 0:ny+1), aw_cds(0:nx+1, 0:ny+1), an_cds(0:nx+1, 0:ny+1))
    allocate (as_cds(0:nx+1, 0:ny+1), Sp_cds(0:nx+1, 0:ny+1), Sc_cds(0:nx+1, 0:ny+1), phi_old(0:nx+1, 0:ny+1))

  end subroutine allocate_memory

  subroutine deallocate_memory()
    implicit none

    deallocate (x, y)
    deallocate (phi, ap_upw, ae_upw, aw_upw)
    deallocate (an_upw, as_upw, Sp_upw, Sc_upw)
    deallocate (ap_cds, ae_cds, aw_cds, an_cds)
    deallocate (as_cds, Sp_cds, Sc_cds, phi_old)

  end subroutine deallocate_memory

  subroutine initialize_fvmsolver()
    implicit none

    ! calculate grid spacing
    dx = 1.0d0 / float(nx)
    dy = 1.0d0 / float(ny)

    ! set convection scheme blending factor
    write(*, '(A)') "Enter the blending factor between 0 to 1: "
    read(*, *) blend_factor
    
    ! set domain velocity, assuming both u and v have same magnitudes
    domain_velocity = 1.0d0

    ! set density of fluid
    rho = 1.0d0

    ! set phi_west and phi_south boundary values
    phi_west = 1.0d0
    phi_south = 0.0d0

    ! calculate mass flow rates (for this problem, it will same for all cells)
    Fe = rho * domain_velocity * dy
    Fw = rho * domain_velocity * dy
    Fn = rho * domain_velocity * dx
    Fs = rho * domain_velocity * dx

    ! initialize source values
    Sp_upw = -(Fe - Fw + Fn - Fs)
    Sc_upw = 0.0d0
    Sp_cds = -(Fe - Fw + Fn - Fs)
    Sc_cds = 0.0d0

    ! intialize the domain
    phi = 0.0d0                 ! set your choice of initialization
    phi_old = phi
    
    ! set 1st order upwind coefficients
    ae_upw = max(-Fe, 0.0d0)
    aw_upw = max(Fw, 0.0d0)
    an_upw = max(-Fn, 0.0d0)
    as_upw = max(Fs, 0.0d0)

    ! set central difference scheme coefficients
    ae_cds = -Fe / 2.0d0
    aw_cds = Fw / 2.0d0
    an_cds = -Fn / 2.0d0
    as_cds = Fs / 2.0d0    
    
  end subroutine initialize_fvmsolver

  subroutine set_coefficients()
    implicit none

    ! set west boundary condition
    aw_upw(1, :) = 0.0d0
    Sc_upw(1, :) = Fw * phi_west
    Sp_upw(1, :) = -(Fe + Fn - Fs)
    aw_cds(1, :) = 0.0d0
    Sc_cds(1, :) = Fw * phi_west
    Sp_cds(1, :) = -(Fe + Fn - Fs)

    ! set east boundary condition
    ae_upw(nx, :) = 0.0d0
    Sc_upw(nx, :) = 0.0d0
    Sp_upw(nx, :) = -(Fe - Fw + Fn - Fs)
    ae_cds(nx, :) = 0.0d0
    Sc_upw(nx, :) = 0.0d0
    Sp_upw(nx, :) = -(Fe - Fw + Fn - Fs)

    ! set south boundary condition
    as_upw(:, 1) = 0.0d0
    Sc_upw(:, 1) = Fs * phi_south
    Sp_upw(:, 1) = -(Fe - Fw + Fn)
    as_cds(:, 1) = 0.0d0
    Sc_cds(:, 1) = Fs * phi_south
    Sp_cds(:, 1) = -(Fe - Fw + Fn)

    ! set north boundary condition
    an_upw(:, ny) = 0.0d0
    Sc_upw(:, ny) = 0.0d0
    Sp_upw(:, ny) = -(Fe - Fw + Fn - Fs)
    an_cds(:, ny) = 0.0d0
    Sc_cds(:, ny) = 0.0d0
    Sp_cds(:, ny) = -(Fe - Fw + Fn - Fs)

    ! set southwest cell coefficient
    aw_upw(1, 1) = 0.0d0
    as_upw(1, 1) = 0.0d0
    Sc_upw(1, 1) = Fw * phi_west + Fs * phi_south
    Sp_upw(1, 1) = -(Fe + Fn)
    aw_cds(1, 1) = 0.0d0
    as_cds(1, 1) = 0.0d0
    Sc_upw(1, 1) = Fw * phi_west + Fs * phi_south
    Sp_cds(1, 1) = -(Fe + Fn)

    ! set northwest cell coefficient
    aw_upw(1, ny) = 0.0d0
    an_upw(1, ny) = 0.0d0
    Sc_upw(1, ny) = Fw * phi_west
    Sp_upw(1, ny) = -(Fe + Fn - Fs)
    aw_cds(1, ny) = 0.0d0
    an_cds(1, ny) = 0.0d0
    Sc_cds(1, ny) = Fw * phi_west
    Sp_cds(1, ny) = -(Fe + Fn - Fs)

    ! set southeast cell coefficient
    ae_upw(nx, 1) = 0.0d0
    as_upw(nx, 1) = 0.0d0
    Sc_upw(nx, 1) = Fs * phi_south
    Sp_upw(nx, 1) = -(Fe - Fw + Fn)
    ae_cds(nx, 1) = 0.0d0
    as_cds(nx, 1) = 0.0d0
    Sc_cds(nx, 1) = Fs * phi_south
    Sp_cds(nx, 1) = -(Fe - Fw + Fn)

    ! set northeast cell coefficient
    ae_upw(nx, ny) = 0.0d0
    an_upw(nx, ny) = 0.0d0
    Sc_upw(nx, ny) = 0.0d0
    Sp_upw(nx, ny) = -(Fe - Fw + Fn - Fs)
    ae_cds(nx, ny) = 0.0d0
    an_cds(nx, ny) = 0.0d0
    Sc_cds(nx, ny) = 0.0d0
    Sp_cds(nx, ny) = -(Fe - Fw + Fn - Fs)

    ! compute ap for 1st order upwind and central difference schemes
    ap_upw = ae_upw + aw_upw + an_upw + as_upw - Sp_upw
    ap_cds = ae_cds + aw_cds + an_cds + as_cds - Sp_cds

  end subroutine set_coefficients

  subroutine gs_solver()
    implicit none

    iter_max = 100000           ! maximum number of iterations should be large enough to ensure convergence upto desired tolerance
    ! iterative solver Gauss-Seidel
    do iter = 1, iter_max
       do j = 1, ny
          do i = 1, nx

             phi(i, j) = ((ae_upw(i, j) * phi(i+1, j) + aw_upw(i, j) * phi(i-1, j) + an_upw(i, j) * phi(i, j+1) + &
                  as_upw(i, j) * phi(i, j-1) + Sc_upw(i, j)) / ap_upw(i, j)) - ((blend_factor / ap_upw(i, j)) * &
                  ((ap_cds(i, j) * phi(i, j) - ae_cds(i, j) * phi(i+1, j) - aw_cds(i, j) * phi(i-1, j) - &
                  an_cds(i, j) * phi(i, j+1) - as_cds(i, j) * phi(i, j-1) - Sc_cds(i, j)) - (ap_upw(i, j) * phi(i, j) - &
                  ae_upw(i, j) * phi(i+1, j) - aw_upw(i, j) * phi(i-1, j) - an_upw(i, j) * phi(i, j+1) - as_upw(i, j) * &
                  phi(i, j-1) - Sc_upw(i, j))))

          end do
       end do
    end do
    
  end subroutine gs_solver

  subroutine post_process_data()
    implicit none

    open(unit = 3, file = "converged_results_convection.csv")

    ! compute x and y coords of cell centers
    y(1) = 0.5 * dy
    x(1) = 0.5 * dx

    do j = 2, ny
       y(j) = y(j-1) + dy
    end do

    do i = 2, nx
       x(i) = x(i-1) + dx
    end do

    ! writing interior points only
    write(3, *) "x coord, y coord, z coord, phi"
    do j = 1, ny
       do i = 1, nx
          z = 0.0d0
          write(3, *) x(i), ",", y(j), ",", z, ",", phi(i,j)
       end do
    end do

  end subroutine post_process_data
  
end module conv_prob_data

program fvm_convection

  use conv_prob_data
  implicit none

  call allocate_memory()
  call initialize_fvmsolver()
  call set_coefficients()
  call gs_solver()
  call post_process_data()

  write(*, '(A, I2, A, I2)') "The grid resolution used is ", nx, "x", ny
  write(*, '(A)') "The iterations are computed and solution has converged. EXITING!!!"
  
  call deallocate_memory()

end program fvm_convection
