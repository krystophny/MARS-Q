program test_trapped_orbit_reflection
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    integer, parameter :: n_half_test = 4
    integer, parameter :: n_full_test = 2*n_half_test - 1
    integer :: iout, jhalf, leg, endpoint, owner, closure, status
    integer :: upper_owner_count
    real(real64), parameter :: tol = 1.0e-13_real64
    real(real64) :: chi_half(n_half_test), phi_half(n_half_test)
    real(real64) :: tau_half(n_half_test), hchi_half(n_half_test)
    real(real64) :: vpabs_half(n_half_test)
    real(real64) :: expected_tau(n_full_test), expected_measure(n_full_test)
    integer :: expected_j(n_full_test), expected_leg(n_full_test)
    integer :: expected_endpoint(n_full_test)
    real(real64) :: chi, phi, tau, tau_fraction, bounce_angle, bounce_phase
    real(real64) :: hchi, vpstate, measure, orient_state, orient_vpar
    real(real64) :: measure_sum

    chi_half = [0.2_real64, 0.7_real64, 1.4_real64, 2.0_real64]
    phi_half = [0.0_real64, 0.3_real64, 0.9_real64, 1.5_real64]
    tau_half = [0.0_real64, 0.2_real64, 0.8_real64, 1.0_real64]
    hchi_half = [-2.0_real64, -3.0_real64, -4.0_real64, -5.0_real64]
    vpabs_half = [0.0_real64, 0.4_real64, 0.6_real64, 0.0_real64]
    expected_j = [1, 2, 3, 4, 3, 2, 1]
    expected_leg = [1, 1, 1, 1, -1, -1, -1]
    expected_endpoint = [-1, 0, 0, 1, 0, 0, -1]
    expected_tau = [0.0_real64, 0.2_real64, 0.8_real64, 1.0_real64, &
                    1.2_real64, 1.8_real64, 2.0_real64]
    expected_measure = [0.0_real64, 0.2_real64, 0.6_real64, 0.2_real64, &
                        0.2_real64, 0.6_real64, 0.2_real64]

    measure_sum = 0.0_real64
    upper_owner_count = 0
    do iout = 1, n_full_test
        call ktrappedfullpoint(n_half_test, iout, -1.0_real64, 1.0_real64, &
                               acos(-1.0_real64), chi_half, phi_half, tau_half, &
                               hchi_half, vpabs_half, jhalf, leg, chi, phi, tau, &
                               tau_fraction, bounce_angle, bounce_phase, hchi, &
                               vpstate, measure, orient_state, orient_vpar, &
                               endpoint, owner, closure, status)
        if (status /= 0) error stop 1
        if (jhalf /= expected_j(iout)) error stop 2
        if (leg /= expected_leg(iout)) error stop 3
        if (endpoint /= expected_endpoint(iout)) error stop 4
        if (abs(chi - chi_half(jhalf)) > tol) error stop 5
        if (abs(phi - phi_half(jhalf)) > tol) error stop 6
        if (abs(tau - expected_tau(iout)) > tol) error stop 7
        if (abs(measure - expected_measure(iout)) > tol) error stop 8
        if (abs(tau_fraction - 0.5_real64*tau) > tol) error stop 9
        if (abs(bounce_angle - acos(-1.0_real64)*tau) > tol) error stop 10
        if (abs(bounce_phase + bounce_angle) > tol) error stop 11
        if (iout == 1 .or. iout == n_half_test .or. iout == n_full_test) then
            if (abs(vpstate) > tol) error stop 12
            if (abs(orient_state) > tol) error stop 13
            if (abs(orient_vpar) > tol) error stop 14
        else
            if (abs(orient_state - real(leg, real64)) > tol) error stop 15
            if (abs(orient_vpar + real(leg, real64)) > tol) error stop 16
        end if
        if (endpoint == 1 .and. owner == 1) upper_owner_count = upper_owner_count + 1
        if (iout == 1 .and. owner /= 1) error stop 17
        if (iout == n_full_test .and. owner /= 0) error stop 18
        if (closure /= merge(1, 0, iout == n_full_test)) error stop 19
        measure_sum = measure_sum + measure
    end do
    if (upper_owner_count /= 1) error stop 20
    if (abs(measure_sum - 2.0_real64) > tol) error stop 21
end program test_trapped_orbit_reflection

include '../MarsQ_2FK/trapped_orbit_reflection.f90'
