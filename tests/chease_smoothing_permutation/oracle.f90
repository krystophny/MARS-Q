program smoothing_oracle
    use globals
    implicit none
    integer :: i, j, node, mapped, d
    real(rkind) :: s, t, f, df, angle, dangle, error
    ns = 16
    nt = 64
    ns1 = ns + 1
    nt1 = nt + 1
    rc2pi = 2*acos(-1.0_rkind)
    csig(1:ns1) = [(real(i-1, rkind)/ns, i=1,ns1)]
    ct(1:nt1) = [(rc2pi*(i-1)/nt, i=1,nt1)]
    cpsi = huge(1.0_rkind)/4
    do i = 1, ns1
        s = csig(i)
        f = s**3-2*s*s+0.5_rkind*s-1
        df = 3*s*s-4*s+0.5_rkind
        do j = 1, nt
            t = ct(j)
            angle = 1+0.2_rkind*sin(t)
            dangle = 0.2_rkind*cos(t)
            node = (i-1)*nt+j
            ! A nonidentity bijection also rejects an accidental identity copy.
            nupdwn(node) = mod(node+16, ns1*nt)+1
            cpsicl(4*node-3:4*node) = [f*angle, df*angle, f*dangle, df*dangle]
        end do
    end do
    call smooth
    error = 0
    do node = 1, ns1*nt
        mapped = nupdwn(node)
        do d = 1, 4
            error = max(error, abs(cpsi(4*(mapped-1)+d)-cpsicl(4*(node-1)+d)))
        end do
    end do
    if (error /= 0) error stop 'Smoothing lost a physical Hermite jet component'
    ! Cubic radial law and periodic derivative give independent field values.
    do i = 2, ns1
        s = csig(i)
        f = s**3-2*s*s+0.5_rkind*s-1
        df = 3*s*s-4*s+0.5_rkind
        do j = 1, nt
            node = (i-1)*nt+j
            if (abs(cpsicl(4*node-1)-0.2_rkind*f*cos(ct(j))) > 2.0d-6) &
                error stop 'Periodic first-jet oracle failed'
            if (abs(cpsicl(4*node)-0.2_rkind*df*cos(ct(j))) > 2.0d-6) &
                error stop 'Periodic mixed-jet oracle failed'
        end do
    end do
    print *, 'PASS: four native jet components, nonidentity node permutation, cubic/periodic jets'
end program
