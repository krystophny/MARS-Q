program axis_field_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), external :: native_axis_field
    real(dp) :: x(4), y(4), expected, actual
    integer :: grid, profile, failures
    failures = 0
    do grid = 1, 2
        if (grid == 1) then
            x = [0.11_dp, 0.24_dp, 0.39_dp, 0.62_dp]
        else
            x = [0.0078273_dp, 0.023482_dp, 0.039137_dp, 0.054791_dp]
        end if
        do profile = 1, 4
            expected = 1.0_dp
            select case (profile)
            case (1)
                y = 1.0_dp
            case (2)
                y = 1.0_dp + 0.3_dp * x**2
            case (3)
                y = 1.0_dp + 0.3_dp * x**2 - 0.1_dp * x**3
            case (4)
                y = -1.0_dp - 0.3_dp * x**2
                expected = -1.0_dp
            end select
            actual = native_axis_field(y, x)
            if (abs(actual - expected) > 2.0e-13_dp) failures = failures + 1
            write (*, '(2i3,2es25.16)') grid, profile, expected, actual
        end do
    end do
    if (failures /= 0) error stop 'native axis-field polynomial oracle failed'
end program axis_field_oracle
