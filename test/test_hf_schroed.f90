program test_hf_schroed

use types, only: dp
use mesh, only: meshexp
use hf, only: solve_hf_schroed
use feutils, only: get_parent_quad_pts_wts
implicit none

real(dp), allocatable :: xe(:)              ! element coordinates
real(dp), allocatable :: xiq(:), wtq(:)     ! quadrature points and weights
integer :: p, Ne, Nq, Z, i, DOFs
real(dp) :: rmin, rmax, a, err, Etot
real(dp), allocatable :: energies(:), V(:,:)
real(dp), parameter :: Etot_ref = -14.32329062_dp

Z = 4
rmin = 0
rmax = 50
a = 200
Ne = 4
Nq = 53
p = 26

allocate(xe(Ne+1), xiq(Nq), wtq(Nq), V(Nq, Ne))
xe = meshexp(rmin, rmax, a, Ne)
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)

call solve_hf_schroed(Z, p, xiq, wtq, xe, 1e-8_dp, energies, Etot, V, DOFs)

print *
print *, "Hartree-Fock (Slater Exchange) calculation for Be (Z=4)"
print *
print *, "Final Total energy:"
print "(a16,a16,a10)", "E", "E_ref", "error"
err = abs(Etot - Etot_ref)
print "(f16.8, f16.8, es10.2)", Etot, Etot_ref, err
if ( .not. (err < 1e-8_dp)) then
   error stop 'assert failed'
end if
print *
print *, "Final Eigenvalues:"
do i = 1, size(energies)
    print "(i4, f16.8)", i, energies(i)
end do

end program
