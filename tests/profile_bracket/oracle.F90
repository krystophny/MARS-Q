program profile_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    implicit none
    integer, parameter :: capacity = 128
    real(dp) :: knots(capacity), values(capacity), second(capacity)
    real(dp) :: flux(capacity), output(capacity), deriv(capacity)
    real(dp) :: reference(capacity), refderiv(capacity), scale, rho, expected, d_expected
    real(dp) :: t0, t1, checksum
    integer :: k, sign_case, n, intervals, repeat, repeats, pass, mode
    character(16) :: argument
    call get_command_argument(1, argument)
    n = 12
    intervals = 6
    repeats = 1
    if (argument == 'timing') then
        n = 32
        intervals = 64
        repeats = 200000
    end if
    do k = 1, intervals + 1
        knots(k) = (real(k - 1, dp) / intervals)**1.3_dp
        values(k) = polynomial(knots(k))
        second(k) = 0.4_dp - 0.6_dp * knots(k)
    end do
    do sign_case = 1, 2
        scale = -0.8_dp
        if (sign_case == 2) scale = 1.3_dp
        call set_profile_scale(scale)
        do k = 1, n
            rho = 1.4_dp * real(k - 1, dp) / (n - 1)
            flux(k) = scale * (1.0_dp - rho**2)
        end do
        flux(2) = 1.01_dp * scale ! Clamp beyond the axis, independent of scale sign.
        flux(3) = scale * (1.0_dp - knots(3)**2) ! Exact knot neighbourhood.
        do mode = 0, 1
#ifdef MARS
            call ppspln(n, flux, intervals, knots, values, second, output, mode)
            call baseline_ppspln(n, flux, intervals, knots, values, second, reference, mode)
#else
            call ppspln(n, flux, intervals, knots, values, second, output, deriv, 2)
            call baseline_ppspln(n, flux, intervals, knots, values, second, reference, refderiv, 2)
            if (mode == 1) then
                output(1:n) = deriv(1:n)
                reference(1:n) = refderiv(1:n)
            end if
#endif
            do k = 1, n
                rho = sqrt(max(0.0_dp, 1.0_dp - flux(k) / scale))
                expected = polynomial(rho)
                if (mode == 1) expected = 0.3_dp + 0.4_dp * rho - 0.3_dp * rho**2
                if (abs(output(k) - expected) > 3.0e-13_dp) error stop 'exact cubic oracle'
                if (transfer(output(k), 0_int64) /= transfer(reference(k), 0_int64)) &
                    error stop 'native parent byte equality'
            end do
        end do
    end do
    if (argument /= 'timing') stop
    do pass = 1, 6
        checksum = 0.0_dp
        call cpu_time(t0)
        do repeat = 1, repeats
#ifdef MARS
            if (mod(pass, 2) == 1) then
                call baseline_ppspln(n, flux, intervals, knots, values, second, output, 0)
            else
                call ppspln(n, flux, intervals, knots, values, second, output, 0)
            end if
#else
            if (mod(pass, 2) == 1) then
                call baseline_ppspln(n, flux, intervals, knots, values, second, output, deriv, 0)
            else
                call ppspln(n, flux, intervals, knots, values, second, output, deriv, 0)
            end if
#endif
            checksum = checksum + output(1)
        end do
        call cpu_time(t1)
        write (*, '(i2,2es24.15)') pass, t1 - t0, checksum
    end do
contains
    pure function polynomial(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = 0.8_dp + 0.3_dp*x + 0.2_dp*x**2 - 0.1_dp*x**3
    end function polynomial
end program profile_oracle
