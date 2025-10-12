module hf

use types, only: dp
use mesh, only: meshexp
use feutils, only: define_connect, get_quad_pts, get_parent_quad_pts_wts, &
        get_parent_nodes, phih, dphih, c2fullc2, fe2quad_core, get_nodes, &
        integrate, proj_fn, phih_array, dphih_array
use linalg, only: eigh
use fe, only: assemble_radial_H, assemble_radial_S, assemble_radial_H_setup, &
        assemble_radial_H_complete
use constants, only: pi
use hartree_screening, only: assemble_poisson_A, hartree_potential3
use mixings, only: mixing_anderson
use states, only: get_atomic_states_nonrel_focc, nlf2focc, get_atomic_states_nonrel
use energies, only: thomas_fermi_potential
use iso_c_binding, only: c_double, c_int
use solvers, only: solve_eig_irange, solve_sym_setup
use lapack, only: dsytrf, dsytrs
implicit none
private
public solve_hf_schroed

contains

subroutine solve_hf_schroed(Z, p, xiq, wtq, xe, eps, energies, Etot, V, DOFs)

    integer, intent(in) :: Z, p
    real(dp), intent(in) :: xe(:)        ! element coordinates
    real(dp), intent(in) :: xiq(:)       ! quadrature points
    real(dp), intent(in) :: wtq(:)       ! quadrature weights
    real(dp), intent(in) :: eps
    real(dp), allocatable, intent(out) :: energies(:)
    real(dp), intent(out) :: Etot
    real(dp), intent(out) :: V(:,:)       ! SCF potential
    integer, intent(out) :: DOFs
    integer :: n, Nq
    real(dp), allocatable :: H(:,:), S(:), D(:,:), lam(:)

    real(dp), allocatable :: xin(:)       ! parent basis nodes
    integer, allocatable :: ib(:, :)       ! basis connectivity
    integer, allocatable :: in(:, :)
    real(dp), allocatable :: xq(:, :), fullc(:), uq(:,:), rho(:,:), Vee(:,:), &
        Vx(:,:), Vin(:,:), Vout(:,:), xn(:), xq1(:,:), Am_p(:,:), &
        bv_p(:), Hl(:,:,:)
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: phihq(:,:)   ! parent basis at quadrature points
    real(dp), allocatable :: dphihq(:,:)  ! parent basis derivative at quadrature points

    real(dp), allocatable :: xin_p(:)       ! parent basis nodes for poisson grid
    real(dp), allocatable :: xn_p(:)
    real(dp), allocatable :: phihq_p(:,:)   ! parent basis at quadrature points
    real(dp), allocatable :: dphihq_p(:,:)  ! parent basis derivative at quadrature points
    integer, allocatable :: in_p(:, :), ib_p(:, :)

    integer :: Ne, Nb, Nn, Nb_p, Nn_p, pp
    integer :: l, i, j, Lmax, al, eimin, eimax
    real(dp), allocatable :: focc(:,:), tmp(:)
    real(dp) :: scf_alpha, scf_L2_eps, scf_eig_eps
    integer, allocatable :: no(:), lo(:), focc_idx(:,:)
    real(dp), allocatable :: fo(:)
    real(dp) :: T_s, E_ee, E_en, E_x
    real(dp) :: xiq_lob(p+1), wtq_lob(p+1)


    integer :: nband, scf_max_iter, iter

    Nq = size(xiq)
    Ne = size(xe)-1

    call get_parent_quad_pts_wts(2, p+1, xiq_lob, wtq_lob)

    call get_atomic_states_nonrel_focc(Z, focc)
    call get_atomic_states_nonrel(Z, no, lo, fo)

    Lmax = ubound(focc,2)

    Nn = Ne*p+1

    allocate(xin(p+1))
    call get_parent_nodes(2, p, xin)
    allocate(in(p+1, Ne), ib(p+1, Ne))
    call define_connect(1, 1, Ne, p, in, ib)
    Nb = maxval(ib)
    if ( .not. (Nn == maxval(in)) ) then
       error stop 'Wrong size for Nn'
    end if
    DOFs = Nb
    allocate(xq(Nq, Ne), xn(Nn))
    call get_quad_pts(xe, xiq, xq)
    call get_nodes(xe, xin, xn)

    pp = 2*p
    Nn_p = Ne*pp+1
    allocate(xin_p(pp+1))
    call get_parent_nodes(2, pp, xin_p)
    allocate(in_p(pp+1, Ne), ib_p(pp+1, Ne))
    call define_connect(2, 1, Ne, pp, in_p, ib_p)
    Nb_p = maxval(ib_p)
    if ( .not. (Nn_p == maxval(in_p)) ) then
       error stop 'Wrong size for Nn_p'
    end if
    allocate(xn_p(Nn_p))
    call get_nodes(xe, xin_p, xn_p)


    allocate(phihq(Nq, p+1))
    allocate(dphihq(Nq, p+1))
    allocate(phihq_p(Nq, pp+1))
    allocate(dphihq_p(Nq, pp+1))
    ! tabulate parent basis at quadrature points
    do al = 1, p+1
        call phih_array(xin, al, xiq, phihq(:, al))
        call dphih_array(xin, al, xiq, dphihq(:, al))
    end do
    do al = 1, pp+1
        call phih_array(xin_p, al, xiq, phihq_p(:, al))
        call dphih_array(xin_p, al, xiq, dphihq_p(:, al))
    end do

    n = Nb
    allocate(H(n, n), S(n))
    allocate(D(n, n), lam(n), fullc(Nn), uq(Nq,Ne), rho(Nq,Ne), Vee(Nq,Ne), &
        Vx(Nq,Ne), Vin(Nq,Ne), Vout(Nq,Ne))

    nband = count(focc > 0)
    scf_max_iter = 100
    scf_alpha = 0.25_dp ! A conservative alpha is good for Anderson
    scf_L2_eps = 1e-7_dp
    scf_eig_eps = eps

    allocate(tmp(Nq*Ne))
    allocate(energies(nband))
    Vin = reshape(thomas_fermi_potential(reshape(xq, [Nq*Ne]), Z), [Nq, Ne]) + &
        Z / xq
    iter = 0
    Etot = 0.0_dp
    print *, "Starting SCF calculation..."
    print *, "Iter |   Total Energy   |  dEnergy  |  dV (L2)"

    call assemble_radial_S(xin, xe, ib, wtq_lob, S)
    do i = 1, size(S)
        S(i) = 1/sqrt(S(i))
    end do

    allocate(Am_p(Nb_p, Nb_p), bv_p(Nb_p), ipiv(Nb_p))
    Am_p = 0
    call assemble_poisson_A(xin_p, xe, ib_p, xiq, wtq, dphihq_p, Am_p)
    call solve_sym_setup(Am_p, ipiv)

    allocate(Hl(n,n,0:Lmax))
    call assemble_radial_H_setup(0, Lmax, xin, xe, ib, xiq, wtq, phihq, dphihq, Hl)

    allocate(focc_idx(size(focc,1), 0:Lmax))
    focc_idx = 0
    j=1
    do l=0, Lmax
        do i=1,size(focc,1)
            if (focc(i,l) > 0) then
                focc_idx(i,l) = j
                j = j + 1
            end if
        end do
    end do

    call mixing_anderson &
        (Ffunc, integral, reshape(Vin, [Nq*Ne]), &
        nband, scf_max_iter, scf_alpha, scf_L2_eps, scf_eig_eps, tmp)

contains

    subroutine Ffunc(x, y, eng)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:), eng(:)
    integer :: idx
    real(dp) :: Etot_old
    Etot_old = Etot
    iter = iter + 1
    Vin = reshape(x, shape(Vin))
    rho = 0
    idx = 0
    V = Vin - Z/xq

    do l = 0, Lmax
        call assemble_radial_H_complete(V, xin, xe, ib, xiq, wtq, phihq, Hl(:,:,l), H)

        do concurrent (i = 1:size(S), j = 1:size(S), i>=j)
            H(i, j) = H(i, j)*S(i)*S(j)
        end do

        eimax = 0
        do i=1, size(focc,1)
            if (focc(i,l) > 0) eimax = i
        end do
        if (eimax == 0) cycle

        call solve_eig_irange(H, 1, eimax, lam, D)

        do i = 1, size(S)
            D(i, 1:eimax) = D(i, 1:eimax)*S(i)
        end do

        do i = 1, eimax
            if (focc(i,l) < tiny(1._dp)) cycle

            call c2fullc2(in, ib, D(:Nb,i), fullc)
            call fe2quad_core(xe, xin, in, fullc, phihq, uq)
            rho = rho - focc(i,l)*(uq/xq)**2 / (4*pi)

            idx = focc_idx(i,l)
            eng(idx) = lam(i)
        end do
    end do
    energies = eng

    Vee = hartree_potential3(real(Z, dp), pp, xe, xin_p, in_p, ib_p, xiq, wtq, phihq_p, Am_p, ipiv, -4*pi*rho)
    call slater_exchange_potential(rho, Vx)
    call hf_total_energy_slater(xe, wtq, fo, energies, Vin-Z/xq, Vee, -Z/xq, xq, -rho, T_s, E_ee, E_en, E_x, Etot)

    Vout = Vee + Vx

    if (iter > 1) then
        print '(i4, " | ", f16.8, " | ", es9.2, " | ", es9.2)', iter, Etot, Etot-Etot_old, sqrt(integrate(xe, wtq, (Vout-Vin)**2))
    else
        print '(i4, " | ", f16.8, " | ", a9, " | ", a9)', iter, Etot, '---', '---'
    end if

    y = reshape(Vout, shape(y))
    end subroutine

    real(dp) function integral(x)
    real(dp), intent(in) :: x(:)
    integral = integrate(xe, wtq, reshape(x, shape(uq)))
    end function
end subroutine


subroutine slater_exchange_potential(rho, Vx)
    ! Calculates the Slater local exchange potential (alpha=1)
    ! V_x = - (3/2) * (3/pi)^(1/3) * rho^(1/3)
    real(dp), intent(in) :: rho(:,:)
    real(dp), intent(out) :: Vx(:,:)
    real(dp), parameter :: C_x = -1.5_dp * (3.0_dp / pi)**(1.0_dp/3.0_dp)

    Vx = C_x * sign(1.0_dp, rho) * abs(rho)**(1.0_dp/3.0_dp)
end subroutine


subroutine hf_total_energy_slater(xe, wtq, fo, ks_energies, V_in, V_h, V_coulomb, &
    R, n, T_s, E_ee, E_en, E_x, Etot)
    real(dp), intent(in) :: xe(:), wtq(:)
    real(dp), intent(in) :: R(:,:)
    real(dp), intent(in) :: fo(:), ks_energies(:)
    real(dp), intent(in) :: V_in(:,:), V_h(:,:), V_coulomb(:,:)
    real(dp), intent(in) :: n(:,:)
    real(dp), intent(out) :: Etot
    real(dp), intent(out) :: T_s, E_ee, E_en, E_x

    real(dp) :: rho(size(n,1), size(n,2))
    real(dp) :: E_band
    real(dp) :: V_x(size(n,1), size(n,2))

    rho = -n

    call slater_exchange_potential(rho, V_x)
    E_x = 0.75_dp * 4*pi * integrate(xe, wtq, rho * V_x * R**2)

    E_band = sum(fo * ks_energies)
    T_s = E_band + 4*pi * integrate(xe, wtq, (V_in) * rho * R**2)

    E_ee = -2*pi * integrate(xe, wtq, V_h * rho * R**2)
    E_en =  4*pi * integrate(xe, wtq, (-V_coulomb) * rho * R**2)

    Etot = T_s + E_ee + E_en + E_x
end subroutine

end module
