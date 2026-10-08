program band_oracle
    use prec_const
    implicit none
    integer :: n, m, mp, i, j, kindcase, info, baseline_info
    real(rkind), allocatable :: dense(:,:), band(:,:), old(:,:), rhs(:), x(:), exact(:)
    real(rkind) :: a1(1), empty(1)
    character(len=16) :: mode
    call get_command_argument(1, mode)
    if (trim(mode) == 'edge') then
        empty = 123
        call aldlt(empty, 1.0d-14, 0, 1, 1, info)
        if (info /= 0 .or. empty(1) /= 123) error stop 'Empty factorization failed'
        a1 = 2
        call aldlt(a1, 1.0d-14, 1, 1, 1, info)
        if (info /= 0 .or. a1(1) /= 2) error stop 'Singleton failed'
        a1 = 0
        call aldlt(a1, 1.0d-14, 1, 1, 1, info)
        if (info /= -1) error stop 'Zero singleton accepted'
        block
            real(rkind) :: singular(2,2)
            singular = 0
            singular(1,2) = 1
            call aldlt(singular, 1.0d-14, 2, 2, 2, info)
            if (info /= -1) error stop 'Zero first pivot accepted'
        end block
        print *, 'PASS empty/singleton/singular'
        stop
    end if
    call panel_threshold_oracle()
    do kindcase = 1, 2
        do n = 3, 75, 8
            do m = 1, n+2, 3
                mp = m+3
                allocate(dense(n,n),band(mp,n),old(mp,n),rhs(n),x(n),exact(n))
                dense = 0
                band = 0
                do i = 1, n
                    do j = i+1, min(n,i+m-1)
                        dense(i,j) = 0.1_rkind*sin(real(3*i+j,rkind))
                        dense(j,i) = dense(i,j)
                    end do
                end do
                do i = 1, n
                    dense(i,i) = 2+sum(abs(dense(i,:)))
                    if (kindcase == 2 .and. mod(i,2) == 0) dense(i,i) = -dense(i,i)
                    band(1:min(m,n-i+1),i) = dense(i,i:min(n,i+m-1))
                    rhs(i) = cos(real(i,rkind))
                end do
                old = band
                baseline_info = 0
            call baseline_aldlt(old,1.0d-14,n,m,mp,baseline_info)
                call aldlt(band,1.0d-14,n,m,mp,info)
                if (info /= 0 .or. info /= baseline_info) error stop 'Factorization rejected regular matrix'
                if (any(transfer(band, [0_8], size(band)) /= transfer(old, [0_8], size(old)))) error stop 'Finite factors changed versus native baseline'
                x = rhs
                call lyv(band,x,n,n,m,mp)
                call dwy(band,x,n,n,m,mp)
                call ltxw(band,x,n,n,m,mp)
                call dense_pivot_solve(dense,rhs,exact)
                if (maxval(abs(x-exact)) > 1.0d-13) error stop 'Independent dense solution differs'
                if (maxval(abs(matmul(dense,x)-rhs)) > 2.0d-13) error stop 'Physical linear residual too large'
                deallocate(dense,band,old,rhs,x,exact)
            end do
        end do
    end do
    print *, 'PASS SPD/indefinite, narrow/wide/padded bands, bitwise factors and dense pivot oracle'
contains
    subroutine panel_threshold_oracle()
        real(rkind) :: a(3,34), old(3,34), delta
        integer :: which, i, status, reference, expected
        do which = 1, 2
            a = 0
            do i = 1, 34
                a(1,i) = 1
            end do
            delta = 0.64_rkind*1.0e-10_rkind
            if (which == 1) then
                delta = delta/2
                expected = -1
            else
                delta = delta*2
                expected = 0
            end if
            ! Pivot 33 receives updates from both sides of the panel boundary.
            ! Its threshold refers to the diagonal before pivot 32 updates it.
            a(3,31) = 0.6_rkind
            a(2,32) = sqrt(0.64_rkind-delta)
            old = a
            reference = 0
            call baseline_aldlt(old,1.0e-10_rkind,34,3,3,reference)
            call aldlt(a,1.0e-10_rkind,34,3,3,status)
            if (status /= expected .or. reference /= expected) &
                error stop 'Panel boundary pivot threshold changed'
            if (status == 0) then
                if (any(transfer(a,[0_8],size(a)) /= &
                        transfer(old,[0_8],size(old)))) &
                    error stop 'Accepted panel threshold factors changed'
            end if
        end do
    end subroutine
    subroutine dense_pivot_solve(a,b,x)
        real(rkind), intent(in) :: a(:,:), b(:)
        real(rkind), intent(out) :: x(:)
        real(rkind) :: work(size(b),size(b)), y(size(b)), row(size(b)), value, factor
        integer :: i, j, k, pivot, nn
        nn = size(b)
        work = a
        y = b
        do k = 1, nn-1
            pivot = k-1+maxloc(abs(work(k:nn,k)),dim=1)
            row = work(k,:)
            work(k,:) = work(pivot,:)
            work(pivot,:) = row
            value = y(k)
            y(k) = y(pivot)
            y(pivot) = value
            do i = k+1, nn
                factor = work(i,k)/work(k,k)
                do j = k+1, nn
                    work(i,j) = work(i,j)-factor*work(k,j)
                end do
                y(i) = y(i)-factor*y(k)
            end do
        end do
        do i = nn, 1, -1
            x(i) = (y(i)-dot_product(work(i,i+1:nn),x(i+1:nn)))/work(i,i)
        end do
    end subroutine
end program
