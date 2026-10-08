! Fixture module supplies only variables used by the unmodified native routines.
module prec_const
    implicit none
    integer, parameter :: rkind = kind(1.0d0)
end module
module globals
    use prec_const
    implicit none
    integer, parameter :: npt = 64, nsp1 = 17
    integer :: ns, nt, ns1, nt1, nupdwn(nsp1*npt), index_out = 1
    real(rkind) :: cpsi(4*nsp1*npt), cpsicl(4*nsp1*npt)
    real(rkind) :: csig(nsp1), ct(npt+1), rc2pi, rc1m14 = 1.0d-14
    real(rkind) :: a, b
    type codeparam_type
        character(len=128) :: output_diag(2)
        integer :: output_flag = 0
    end type
    type equilibrium_type
        type(codeparam_type) :: codeparam
    end type
    type(equilibrium_type) :: eqchease_out(1)
end module
