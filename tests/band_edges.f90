program band_edge_oracle
    use prec_const, only: rkind
    implicit none
    real(rkind) :: scalar(1), empty(1), singular(2, 2)
    real(rkind) :: band(2, 3), rhs(3), exact(3), dense(3, 3)
    integer :: info, polarity

    empty = 123
    info = 97
    call aldlt(empty, 1.0e-14_rkind, 0, 1, 1, info)
    if (info /= 0 .or. empty(1) /= 123) error stop 'Empty factorization'
    scalar = 2
    info = 97
    call aldlt(scalar, 1.0e-14_rkind, 1, 1, 1, info)
    if (info /= 0 .or. scalar(1) /= 2) error stop 'Singleton factorization'
    scalar = 0
    call aldlt(scalar, 1.0e-14_rkind, 1, 1, 1, info)
    if (info /= -1) error stop 'Zero singleton accepted'
    singular = 0
    singular(1, 2) = 1
    call aldlt(singular, 1.0e-14_rkind, 2, 2, 2, info)
    if (info /= -1) error stop 'Zero first pivot accepted'
    call aldlt(scalar, 1.0e-14_rkind, -1, 1, 1, info)
    if (info /= -1) error stop 'Negative dimension accepted'
    call aldlt(scalar, 1.0e-14_rkind, 1, 2, 1, info)
    if (info /= -1) error stop 'Insufficient storage accepted'

    ! Exact-zero and small Schur complements still indicate singularity.
    do polarity = 0, 1
        singular = reshape([1._rkind, 1._rkind, &
                            1._rkind + polarity*1.0e-12_rkind, 0._rkind], [2, 2])
        call aldlt(singular, 1.0e-10_rkind, 2, 2, 2, info)
        if (info /= -1) error stop 'Singular final Schur complement accepted'
    end do

    do polarity = -1, 1, 2
        dense = reshape([4., 1., 0., 1., 5., 2., 0., 2., 6.], [3, 3])
        dense(2, 2) = polarity*5._rkind
        band = reshape([4._rkind, 1._rkind, dense(2, 2), &
                        2._rkind, 6._rkind, 0._rkind], [2, 3])
        exact = [1._rkind, -2._rkind, 3._rkind]
        rhs = matmul(dense, exact)
        info = 97
        call aldlt(band, 1.0e-14_rkind, 3, 2, 2, info)
        if (info /= 0) error stop 'Regular matrix status'
        call lyv(band, rhs, 3, 3, 2, 2)
        call dwy(band, rhs, 3, 3, 2, 2, 1)
        call ltxw(band, rhs, 3, 3, 2, 2)
        if (.not. all(abs(rhs-exact) <= 1.0e-12_rkind)) error stop 'Dense oracle'
    end do
    ! Changing one unknown's units must not make an independent diagonal singular.
    dense = 0
    dense(1, 1) = 1
    dense(2, 2) = 1
    dense(3, 3) = 1.0e-12_rkind
    band = 0
    band(1, :) = [1._rkind, 1._rkind, 1.0e-12_rkind]
    rhs = matmul(dense, exact)
    call aldlt(band, 1.0e-10_rkind, 3, 2, 2, info)
    if (info /= 0) error stop 'Scaled nonsingular matrix rejected'
    call lyv(band, rhs, 3, 3, 2, 2)
    call dwy(band, rhs, 3, 3, 2, 2, 1)
    call ltxw(band, rhs, 3, 3, 2, 2)
    if (.not. all(abs(rhs-exact) <= 1.0e-12_rkind)) error stop 'Scaled solution'
    if (maxval(abs(matmul(dense, rhs)-matmul(dense, exact)) / &
               abs(matmul(dense, exact))) > 1.0e-12_rkind) error stop 'Scaled residual'
    print *, 'PASS empty, singleton, singular, dimensions, SPD and indefinite'
end program band_edge_oracle
