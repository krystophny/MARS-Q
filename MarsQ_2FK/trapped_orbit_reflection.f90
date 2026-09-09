! Reflect one native trapped half-bounce sample into an ordered full bounce.
! This routine has no shared state so its map can be tested directly. The
! upper turn is owned once; the final lower-turn row is an unowned periodic
! closure of the initial lower-turn row.
subroutine ktrappedfullpoint(nhalf, iout, ell, tauend, omegab, &
                             chihalf, phihalf, tauhalf, hchihalf, vpabshalf, &
                             jhalf, leg, chi, phi, tau, taufraction, &
                             bounceangle, bouncephase, hchifactor, vpstate, &
                             dtaumeasure, orientstate, orientvpar, endpoint, &
                             endpointowner, periodicclosure, kstatus)
    implicit none

    integer, intent(in) :: nhalf, iout
    real(8), intent(in) :: ell, tauend, omegab
    real(8), intent(in) :: chihalf(*), phihalf(*), tauhalf(*)
    real(8), intent(in) :: hchihalf(*), vpabshalf(*)
    integer, intent(out) :: jhalf, leg, endpoint, endpointowner
    integer, intent(out) :: periodicclosure, kstatus
    real(8), intent(out) :: chi, phi, tau, taufraction, bounceangle
    real(8), intent(out) :: bouncephase, hchifactor, vpstate, dtaumeasure
    real(8), intent(out) :: orientstate, orientvpar

    kstatus = 1
    if (nhalf < 2) return
    if (iout < 1 .or. iout > 2*nhalf - 1) return
    if (tauend <= 0.0_8 .or. omegab <= 0.0_8) return

    if (iout <= nhalf) then
        jhalf = iout
        leg = 1
        tau = tauhalf(jhalf)
        if (iout == 1) then
            dtaumeasure = 0.0_8
        else
            dtaumeasure = tauhalf(jhalf) - tauhalf(jhalf - 1)
        end if
    else
        jhalf = 2*nhalf - iout
        leg = -1
        tau = 2.0_8*tauend - tauhalf(jhalf)
        dtaumeasure = tauhalf(jhalf + 1) - tauhalf(jhalf)
    end if
    if (dtaumeasure < 0.0_8) return

    chi = chihalf(jhalf)
    phi = phihalf(jhalf)
    taufraction = tau/(2.0_8*tauend)
    bounceangle = omegab*tau
    bouncephase = ell*bounceangle
    hchifactor = hchihalf(jhalf)
    vpstate = real(leg, 8)*vpabshalf(jhalf)

    endpoint = 0
    if (jhalf == 1) endpoint = -1
    if (jhalf == nhalf) endpoint = 1
    endpointowner = 0
    if (iout == 1 .or. iout == nhalf) endpointowner = 1
    periodicclosure = 0
    if (iout == 2*nhalf - 1) periodicclosure = 1

    if (abs(vpstate) <= tiny(vpstate) .or. &
        abs(hchifactor) <= tiny(hchifactor)) then
        orientstate = 0.0_8
        orientvpar = 0.0_8
    else
        orientstate = sign(1.0_8, vpstate)
        orientvpar = sign(1.0_8, vpstate*hchifactor)
    end if
    kstatus = 0
end subroutine ktrappedfullpoint
