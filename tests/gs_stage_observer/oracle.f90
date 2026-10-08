program observer_oracle
    use gs_stage_observer
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    integer, parameter :: ns=2, nt=4, n=48, width=3, lda=7
    real(real64) :: a(lda,n), original(lda,n), b(n), x(n), context(8)
    real(real64) :: read_a(width,n), read_b(n), read_x(n)
    real(real64) :: s(ns+1), t(nt+1), theta(3), r(3), z(3), second(3)
    integer :: up(n/4), read_up(n/4), dims4(4), dims3(3), i, unit
    character(1024) :: directory
    call get_command_argument(1, directory)
    a = -1.0e200_real64 ! Padding must not enter the serialized active band.
    do i=1,n
        a(1,i)=4.0_real64
        a(2,i)=-0.25_real64
        a(3,i)=0.125_real64
    end do
    original=a
    do i=1,n
        x(i)=real(i,real64)/7.0_real64
        b(i)=4.0_real64*x(i)
        if (i<n) b(i)=b(i)-0.25_real64*real(i+1,real64)/7.0_real64
        if (i>1) b(i)=b(i)-0.25_real64*real(i-1,real64)/7.0_real64
    end do
    do i=1,n/4
        up(i)=n/4+1-i
    end do
    context=[1.0_real64,0.0_real64,1.1_real64,0.01_real64,-0.5_real64, &
             6.2_real64,5.3_real64,1.0_real64]
    call record_system(ns,nt,n,width,lda,a,context)
    s=[0.0_real64,0.5_real64,1.0_real64]
    t=[0.0_real64,1.0_real64,2.0_real64,3.0_real64,4.0_real64]
    theta=[0.0_real64,3.0_real64,6.0_real64]
    r=[1.1_real64,0.9_real64,1.1_real64]
    z=0.0_real64
    second=0.0_real64
    call record_chart(ns,nt,3,s,t,theta,r,z,second,second)
    call record_rhs(ns,nt,n,b,context)
    call record_vector(ns,nt,n,up,x,'reduced',context)
    call record_check(ns,nt,n,width,lda,a,b,context)
    if (any(a/=original)) error stop 'observer mutated original matrix'
    if (directory=='disabled') stop
    open(newunit=unit,file=trim(directory)//'/g0001_matrix.dat',form='unformatted')
    read(unit) dims4
    if (any(dims4/=[ns,nt,n,width])) error stop 'original matrix dimensions'
    read(unit) read_a
    close(unit)
    if (any(read_a/=a(1:width,:))) error stop 'original matrix bytes or padding'
    open(newunit=unit,file=trim(directory)//'/g0001_s000001_rhs.dat',form='unformatted')
    read(unit) dims3
    if (any(dims3/=[ns,nt,n])) error stop 'original RHS dimensions'
    read(unit) read_b
    close(unit)
    if (any(read_b/=b)) error stop 'original RHS bytes'
    open(newunit=unit,file=trim(directory)//'/g0001_s000001_reduced.dat',form='unformatted')
    read(unit) dims3
    read(unit) read_up
    read(unit) read_x
    close(unit)
    if (any(read_up/=up) .or. any(read_x/=x)) error stop 'reduced snapshot bytes'
end program observer_oracle
