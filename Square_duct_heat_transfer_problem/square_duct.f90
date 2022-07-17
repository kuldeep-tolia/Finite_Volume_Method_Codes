!This is a code to solve for Temperature profile in a square duct with constant heat flux at boundaries for a given cross section

subroutine l2norm(prev, current, n, Tl2norm)

implicit none

integer:: i, j, n
double precision:: num = 0.0d0, denom = 0.0d0, Tl2norm
double precision:: prev(1:n, 1:n), current(1:n, 1:n)

do i = 1, n
	do j = 1, n
		num = num + (prev(i, j) - current(i, j))**2
		denom = denom + (prev(i, j))**2
	end do
end do

Tl2norm = sqrt(num) / sqrt(denom)

end subroutine l2norm

subroutine TDMA(maind, supd, subd, rhsv, phi, n2, j)

implicit none

integer:: i, n2, j
double precision:: maind(n2), supd(n2), subd(n2), rhsv(n2), phi(1:n2, 1:n2), P(n2), Q(n2)

!Initialise P and Q
P(1) = supd(1) / maind(1)
Q(1) = rhsv(1) / maind(1)

!Use recurrence relations to calculate P(i) and Q(i)
do i = 2, n2
	P(i) = supd(i) / (maind(i) - P(i-1) * subd(i))
	Q(i) = (rhsv(i) + Q(i-1) * subd(i)) / (maind(i) - P(i-1) * subd(i))
end do

phi(n2, j) = Q(n2)

!Solve for phi using back substitution
do i = n2-1, 1, -1
	phi(i, j) = Q(i) + P(i) * phi(i+1, j)
end do

end subroutine TDMA

!Main Program
program duct

implicit none

!Declaring variable
integer:: i, j, n1, iter, m
double precision, dimension(:, :), allocatable:: aP, aE, aW, aN, aS, s, T, Told
double precision:: alpha, L, q_flux, dz, dy, dv, sc, k, tolerance, Tl2norm, Tref, Tb, Nu, Twbar
double precision, dimension(:), allocatable:: a, b, c, d, P, Q, Tef, Twf, Tnf, Tsf, theta, Tcenter, z

!Setting the number of grid points for square duct
n1 = 10

!Allocating memory for variables
allocate (aP(1:n1, 1:n1), aE(1:n1, 1:n1), aW(1:n1, 1:n1), aN(1:n1, 1:n1), aS(1:n1, 1:n1))
allocate (s(1:n1, 1:n1), T(1:n1, 1:n1), Told(1:n1, 1:n1))
allocate (a(n1), b(n1), c(n1), d(n1), P(n1), Q(n1))
allocate (Tef(n1), Twf(n1), Tnf(n1), Tsf(n1), theta(n1+2), z(n1+2), Tcenter(n1+2))

!Initialising problem parameters
L = 1.0d0
q_flux = 1.0d0
k = 1.0d0
dz = L / n1
dy = L / n1


!Relaxation factor
alpha = 6.0d-1

!Initialising parameters
dv = dz * dy
sc = -4.0d0 * (q_flux / L) * dv
tolerance = 1.0d-9
Tl2norm = 1.0d0
iter = 1


!Setting reference temperature for a cell
Tref = 400.0d0

!Initialising temperature field --- Input the initial guess value here
T(:, :) = 100.0d0
Told(:, :) = T(:, :)

!Setting up the coefficients after discretisation
!For interior cells
aE(2:n1-1, 2:n1-1) = k * dy / dz
aW(2:n1-1, 2:n1-1) = k * dy / dz
aN(2:n1-1, 2:n1-1) = k * dz / dy
aS(2:n1-1, 2:n1-1) = k * dz / dy
s(2:n1-1, 2:n1-1) = sc

!For boundary cells excluding corner cells
!North face
aE(1, 2:n1-1) = k * dy / dz
aW(1, 2:n1-1) = k * dy / dz
aS(1, 2:n1-1) = k * dz / dy
aN(1, 2:n1-1) = 0.0d0
s(1, 2:n1-1) = sc + (q_flux * dz)

!West face
aE(2:n1-1, 1) = k * dy / dz
aW(2:n1-1, 1) = 0.0d0
aS(2:n1-1, 1) = k * dz / dy
aN(2:n1-1, 1) = k * dz / dy
s(2:n1-1, 1) = sc + (q_flux * dy)

!South face
aE(n1, 2:n1-1) = k * dy / dz
aW(n1, 2:n1-1) = k * dy / dz
aS(n1, 2:n1-1) = 0.0d0
aN(n1, 2:n1-1) = k * dz / dy
s(n1, 2:n1-1) = sc + (q_flux * dz)

!East face
aE(2:n1-1, n1) = 0.0d0
aW(2:n1-1, n1) = k * dy / dz
aS(2:n1-1, n1) = k * dz / dy
aN(2:n1-1, n1) = k * dz / dy
s(2:n1-1, n1) = sc + (q_flux * dy)

!North-East corner cell
aE(1, n1) = 0.0d0
aW(1, n1) = k * dy / dz
aN(1, n1) = 0.0d0
aS(1, n1) = k * dz/ dy
s(1, n1) = sc + (q_flux * dz) + (q_flux * dy)

!South-East corner cell
aE(n1, n1) = 0.0d0
aW(n1, n1) = k * dy / dz
aN(n1, n1) = k * dz/ dy
aS(n1, n1) = 0.0d0
s(n1, n1) = sc + (q_flux * dz) + (q_flux * dy)

!North-West corner cell
aE(1, 1) = k * dy / dz
aW(1, 1) = 0.0d0
aN(1, 1) = 0.0d0
aS(1, 1) = k * dz/ dy
s(1, 1) = sc + (q_flux * dz) + (q_flux * dy)

!South-West corner cell
aE(n1, 1) = k * dy / dz
aW(n1, 1) = 0.0d0
aN(n1, 1) = k * dz/ dy
aS(n1, 1) = 0.0d0
s(n1, 1) = sc + (q_flux * dz) + (q_flux * dy)

aP(:, :) = aE(:, :) + aW(:, :) + aN(:, :) + aS(:, :)

T(1, 1) = Tref

do while (Tl2norm > tolerance)

	do j = 1, n1
		do i = 1, n1

			a(i) = aP(i, j)
			b(i) = aS(i, j)
			c(i) = aN(i, j)
			d(i) = s(i, j)

			if (j > 1) then
				d(i) = d(i) + aW(i, j) * T(i, j-1)
			end if

			if (j < n1) then
				d(i) = d(i) + aE(i, j) * T(i, j+1)
			end if
		end do

		if (j .ne. 1) then

			call TDMA(a, b, c, d, T, n1, j)

		end if

		if (j == 1) then

			T(1, 1) = Tref
			d(2) = d(2) + aN(2, j) * T(1, 1)

			P(2) = b(2) / a(2)
			Q(2) = d(2) / a(2)

			do i = 3, n1
				P(i) = b(i) / (a(i) - c(i) * P(i-1))
				Q(i) = (d(i) + c(i) * Q(i-1)) / (a(i) - c(i) * P(i-1))
			end do

			T(n1, 1) = Q(n1)

			do i = n1-1, 2, -1
				T(i, 1) = Q(i) + P(i) * T(i+1, 1)
			end do
		end if
	end do

	!Applying relaxation
	T(:, :) = (alpha * T(:, :)) + ((1.0d0 - alpha) * Told(:, :))

	Tl2norm = sqrt((T(n1, n1) - Told(n1, n1))**2) / sqrt(T(n1, n1)**2)

	Told(:, :) = T(:, :)

	iter = iter + 1
end do

write(*, *) 'Domain is uniformly initialised using T = 100'
write(*, *) 'Converged solution took iterations =', iter
print*, ''

do i = 1, n1
	write(*, 5) (T(i, j), j = 1, n1)
	5 format(15f10.4)
end do

!-------------Post-processing results-----------------

!Calculating bulk temperature of the cross section
Tb = sum(T(:, :)) * dz * dy / L**2

print*, ''
print*, 'Bulk temperature =', Tb

!Calculating mean wall temperature
Twf(:) = T(:, 1) + (q_flux * dz / (k * 2.0d0))
Tef(:) = T(:, n1) + (q_flux * dz / (k * 2.0d0))
Tnf(:) = T(1, :) + (q_flux * dy / (k * 2.0d0))
Tsf(:) = T(n1, :) + (q_flux * dy / (k * 2.0d0))
Twbar = (sum(Twf(:)) + sum(Tef(:)) + sum(Tnf(:)) + sum(Tsf(:))) / (4 * n1)

print*, ''
print*, 'Mean Wall temeprature =', Twbar

!Calculating Nusselt number of the square duct
Nu = (q_flux * L / k) / (Twbar - Tb)

print*, ''
print*, 'Nusselt number of the square duct =', Nu

!Calculating temperature in center horizontal line
Tcenter(2:n1+1) = (T(n1/2, :) + T((n1/2) + 1, :)) / 2.0d0
Tcenter(1) = (Twf(n1/2) + Twf((n1/2) + 1)) / 2.0d0
Tcenter(n1+2) =  (Tef(n1/2) + Tef((n1/2) + 1)) / 2.0d0
theta(:) = (Tcenter(:) - Tb) / (q_flux * L / k)
z(1) = 0.0d0
z(n1+2) = L
z(2) = dz / 2.0d0
do i = 3, n1+1
	z(i) = z(i-1) + dz
end do

open(17, file = 'theta_profile.txt')
do i = 1, n1+2
	write(17, 3)  z(i) / L, theta(i)
	3 format(2f10.5)
end do

end program duct
