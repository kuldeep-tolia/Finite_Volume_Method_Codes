!Test program for solving system of linear algebraic equations using TDMA

subroutine tdma(a, b, c, d, phi, n)

!a -- main diagonal
!b -- super diagonal
!c -- sub diagonal
!d -- right hand side known vector
!phi -- solution
!n -- number of equations

implicit none

integer:: n
double precision, dimension(n):: a, b, c, d, phi, P, Q
double precision:: temp
integer:: i
	
!Initialise P and Q
P(1) = b(1) / a(1)
Q(1) = d(1) / a(1)
	
!Use recurrence relations to calculate P(i) and Q(i)
do i = 2, n
	temp = a(i) - P(i-1) * c(i)
	P(i) = b(i) / temp
	Q(i) = (d(i) - Q(i-1) * c(i)) / temp
end do
	
phi(n) = Q(n)

!Solve for phi using back substitution
do i = n-1, 1, -1
	phi(i) = Q(i) - P(i) * phi(i+1) 
end do
	
end subroutine tdma

program test_tdma

implicit none

integer, parameter:: n = 5, m = 4
integer:: i
double precision, dimension(:), allocatable:: a1, b1, c1, d1, phi1
double precision, dimension(:), allocatable:: a2, b2, c2, d2, phi2

allocate (a1(1:n), b1(1:n), c1(1:n), d1(1:n), phi1(1:n), a2(1:m), b2(1:m), c2(1:m), d2(1:m), phi2(1:m))

!----------------------Test case 1:------------------------

!Solve the system of linear equations obtained in Assignment-1 Question-3 and it is as follows:
!3phi(1) = phi(2) + 204
!2phi(2) = phi(1) + phi(3) + 12
!2phi(3) = phi(2) + phi(4) + 20
!2phi(4) = phi(3) + phi(5) + 28
!3phi(5) = phi(4) + 1236

!Assigning coefficient values to the arrays
a1(:) = (/ 3.0d0, 2.0d0, 2.0d0, 2.0d0, 3.0d0 /)
b1(:) = (/ -1.0d0, -1.0d0, -1.0d0, -1.0d0, 0.0d0 /)
c1(:) = (/ 0.0d0, -1.0d0, -1.0d0, -1.0d0, -1.0d0 /)
d1(:) = (/ 204.0d0, 12.0d0, 20.0d0, 28.0d0, 1236.0d0 /)


call tdma(a1, b1, c1, d1, phi1, n)

write(*, *) 'Solution for test case 1:'
write(*, *) ''
do i = 1, n
	write(*, 1) phi1(i)
	1 format(f15.6)
end do

!----------------------Test case 2:------------------------

!Solve the system of linear equations as follows:
!104phi(1) = 28phi(2) + 26898
!83.368phi(2) = 28phi(1) + 35.368phi(3) + 6010
!103.368phi(3) = 35.368phi(2) + 48phi(4) + 6010
!77.06phi(4) = 48phi(3) + 8708.56

!Assigning coefficient values to the arrays
a2(:) = (/ 104.0d0, 83.3680d0, 103.3680d0, 77.060d0 /)
b2(:) = (/ -28.0d0, -35.3680d0, -48.0d0, 0.0d0 /)
c2(:) = (/ 0.0d0, -28.0d0, -35.3680d0, -48.0d0 /)
d2(:) = (/ 26898.0d0, 6010.0d0, 6010.0d0, 8708.560d0 /)

call tdma(a2, b2, c2, d2, phi2, m)

write(*, *) ''
write(*, *) 'Solution for test case 2:'
write(*, *) ''
do i = 1, m
	write(*, 1) phi2(i)
end do

deallocate (a1, a2, b1, b2, c1, c2, d1, d2, phi1, phi2)

end program test_tdma
