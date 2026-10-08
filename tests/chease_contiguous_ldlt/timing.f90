program band_timing
    use prec_const
    implicit none
    real(rkind), allocatable :: original(:,:), a(:,:)
    real(rkind) :: t0, t1, old_time, new_time, error
    integer :: n, m, i, j, info
    character(len=32) :: arg
    call get_command_argument(1,arg)
    read(arg,*) n
    call get_command_argument(2,arg)
    read(arg,*) m
    allocate(original(m,n),a(m,n))
    original = 0
    do i = 1, n
        original(1,i) = 2
        do j = 2, min(m,n-i+1)
            original(j,i) = 0.001_rkind*sin(real(i+j,rkind))/m
        end do
    end do
    a = original
    call cpu_time(t0)
    info = 0
    call baseline_aldlt(a,1.0d-14,n,m,m,info)
    call cpu_time(t1)
    if (info /= 0) error stop 'Baseline failed'
    old_time = t1-t0
    original = a
    ! Reconstruct the same independent diagonally dominant input.
    a = 0
    do i = 1, n
        a(1,i) = 2
        do j = 2, min(m,n-i+1)
            a(j,i) = 0.001_rkind*sin(real(i+j,rkind))/m
        end do
    end do
    call cpu_time(t0)
    call aldlt(a,1.0d-14,n,m,m,info)
    call cpu_time(t1)
    if (info /= 0) error stop 'Contiguous factorization failed'
    new_time = t1-t0
    error = maxval(abs(a-original))
    if (error /= 0) error stop 'Timing factors differ'
    print '(a,i0,a,i0,a,es12.5,a,es12.5,a,f9.3,a)', &
        '{"n":',n,',"m":',m,',"baseline_cpu_s":',old_time, &
        ',"contiguous_cpu_s":',new_time,',"speedup":',old_time/new_time,',"factor_max_difference":0}'
end program
