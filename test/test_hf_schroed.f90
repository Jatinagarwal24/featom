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
! Psi4 RHF/6-31G energy for Be: -14.56928462073318364389 Hartrees
! Psi4 RHF/aug-cc-pVQZ energy for Be: -14.57541502891604956460 Hartrees
! Psi4 RHF/heavy-aug-cc-pvdz-dk with x2c energy for Be: -14.57489107478905410176 Hartrees

! Reference values from Psi4 RHF/aug-cc-pVQZ
real(dp), parameter :: Etot_ref = -14.5754150289_dp
real(dp), parameter :: energies_ref(*) = [ &
    -4.7329454091_dp, &   ! 1s eigenvalue
    -0.3093151745_dp  ]   ! 2s eigenvalue

Z = 4
rmin = 0
rmax = 180
a = 12
Ne = 14
Nq = 64
p = 30

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
if ( .not. (err < 1e-2_dp)) then
   error stop 'assert failed'
end if
print *
print *, "Final Eigenvalues:"
do i = 1, size(energies)
    print "(i4, f16.8)", i, energies(i)
end do

if (size(energies) /= size(energies_ref)) then
    print *, "ERROR: Mismatched number of eigenvalues."
    error stop 'Eigenvalue count assert failed'
end if

print "(a4,a16,a16,a10)", "n", "E", "E_ref", "error"
do i = 1, size(energies)
    err = abs(energies(i) - energies_ref(i))
    print "(i4, f16.8, f16.8, es10.2)", i, energies(i), energies_ref(i), err
    ! Assert that the methodical error for eigenvalues is also within range
    if (err > 5.0e-0_dp) then
        error stop 'Eigenvalue assert failed'
    end if
end do

end program
