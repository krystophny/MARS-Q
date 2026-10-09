program band_oracle
    use prec_const
    implicit none
    integer :: n, m, mp, i, j, kindcase, info, reference_info
    real(rkind), allocatable :: dense(:,:), band(:,:), old(:,:), rhs(:), x(:), exact(:)
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
                reference_info = 0
            call dense_factor_reference(old,1.0d-14,n,m,mp,reference_info)
                call aldlt(band,1.0d-14,n,m,mp,info)
                if (info /= 0 .or. info /= reference_info) error stop 'Factorization rejected regular matrix'
                if (any(transfer(band, [0_8], size(band)) /= transfer(old, [0_8], size(old)))) error stop 'Finite factors changed versus dense reference'
                x = rhs
                call lyv(band,x,n,n,m,mp)
                call dwy(band,x,n,n,m,mp)
                call ltxw(band,x,n,n,m,mp)
                call dense_pivot_solve(dense,rhs,exact)
                if (.not. all(abs(x-exact) <= 1.0d-13)) error stop 'Independent dense solution differs'
                if (.not. all(abs(matmul(dense,x)-rhs) <= 2.0d-13)) error stop 'Physical linear residual too large'
                deallocate(dense,band,old,rhs,x,exact)
            end do
        end do
    end do
    print *, 'PASS SPD/indefinite, narrow/wide/padded bands, bitwise factors and dense pivot oracle'
contains
    subroutine dense_factor_reference(a,eps,n,m,mp,info)
        integer, intent(in) :: n,m,mp
        integer, intent(out) :: info
        real(rkind), intent(inout) :: a(mp,n)
        real(rkind), intent(in) :: eps
        real(rkind) :: work(n,n), original(n), threshold, diagonal
        integer :: i,j,k,last
        work=0
        do i=1,n
            do j=i,min(n,i+m-1)
                work(i,j)=a(j-i+1,i)
            end do
        end do
        info=0
        threshold=abs(work(1,1))*eps
        do k=1,n-1
            diagonal=work(k,k)
            if (diagonal==0 .or. abs(diagonal)<threshold) then
                info=-1
                return
            endif
            threshold=abs(work(k+1,k+1))*eps
            last=min(n,k+m-1)
            original=work(k,:)
            work(k,k+1:last)=original(k+1:last)/diagonal
            ! Dense rank-one updates in elimination order, without band addressing.
            do i=k+1,last
                do j=i,last
                    if (original(j)/=0) work(i,j)=work(i,j)+(-original(j))*work(k,i)
                enddo
            enddo
        enddo
        if (work(n,n)==0 .or. abs(work(n,n))<abs(work(n-1,n-1))*eps) info=-1
        do i=1,n
            do j=i,min(n,i+m-1)
                a(j-i+1,i)=work(i,j)
            end do
        end do
    end subroutine
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
            call dense_factor_reference(old,1.0e-10_rkind,34,3,3,reference)
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
