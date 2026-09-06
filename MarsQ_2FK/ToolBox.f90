module ToolBox

    real*8, parameter:: epsilon = 1e-6
    integer intMeshPot
    parameter (intMeshPot=10000)
!    parameter (intMeshPot=6000)
    real*8 intMeshX(intMeshPot),intMeshWeight(intMeshPot)
    real*8, allocatable:: x_adp(:),y_adp(:,:),sp_adp(:,:)
    integer, allocatable:: intpwork(:)
    contains
 
    function assignFreeFileUnit () result (unitNum)
		implicit none
		integer::unitNum
		logical::bOpen
		
		bOpen = .TRUE.

		do unitNum = 10, 10000
		    
		    inquire(UNIT=unitNum, OPENED=bOpen)
		    if (.NOT. bOpen) then
			return
		    end if
		    
		end do
		
		unitNum = -1
		
		write (*,*) 'Can not find free Unit.'
		   
	    end function assignFreeFileUnit
	    
	    subroutine createMesh_old (a,b,step,mesh)
		! create a mesh in the range [a,b]
		implicit none
		   
		real*8,intent(in):: a,b,step
		real*8,allocatable,intent(out):: mesh(:)
		integer meshSize,i
		real*8 range,residua
		
		range = b - a
		meshSize = int( range  / step )
		residua = range - meshSize * step
		meshSize = meshSize + 1
		
		if ( residua - epsilon > 0.0 ) then
		    meshSize = meshSize + 1
		else
		    residua = 0.0
		end if
		
		allocate ( mesh( meshSize ) )
		
		mesh(1) = a
		
		do i = 2, meshSize - 1
		    mesh(i) = a + step * (i-1)            
		end do
		
		if ( residua > epsilon ) then
		    mesh(meshSize) = mesh (meshSize - 1) + residua
		else
		    mesh(meshSize) = mesh (meshSize - 1) + step
		end if
	    end subroutine createMesh_old
	    
	    subroutine createMesh (a,b,step,mesh)
		! create a mesh in the range [a,b]
		implicit none
		   
		real*8,intent(in):: a,b,step
		real*8,allocatable,intent(out):: mesh(:)
		integer meshSize,i
		real*8 range
		real*8 meshStep
		
		range = b - a

	    if (b-a .LE. 1e-6) then
		    allocate (mesh(2))
		    mesh(1) = a
		    mesh(2) = b
		else
		    meshSize = int( range  / step ) + 30
		    meshStep = range / (meshSize - 1)

		    allocate ( mesh( meshSize ) )    		
		    mesh(1) = a
		    do i = 2, meshSize - 1
		        mesh(i) = a + meshStep * (i-1)            
		    end do
		    mesh(meshSize) = b
		end if
		
	    end subroutine createMesh
	    
	    
    subroutine getIntegralPointInCell (kCheck,mesh,intpotCell)
        implicit none
        integer kCheck,fileNo
        real*8,allocatable,intent(in):: mesh(:)
        real*8,allocatable,intent(out):: intpotCell(:,:)
        real*8 eta
        integer cellSize, i
!	    real*8,intent(in)::mesh(:)	
	    real*8 midpot,step
        eta = 1.0 / sqrt(3.0)
        if (kCheck == 1) then
           fileNo = assignFreeFileUnit()
           open (unit=fileNo, file='IntPotCellDebug.OUT',STATUS='unknown')
        
        end if    
        cellSize = size(mesh) - 1
        allocate ( intpotCell ( cellSize, 2 ) )
        if (kCheck == 1) then
            write (fileNo,*) 'mesh=',mesh
            write (fileNo,*) 'cellSize=',cellSize
        end if 
        do i = 1, cellSize
            midpot = 0.5 * ( mesh (i) + mesh (i+1) )
            step = midpot - mesh (i)
            intpotCell (i,1) = midpot - eta * step
            intpotCell (i,2) = midpot + eta * step
        
        end do
        

        if (kCheck == 1) then     
            write (fileNo,*)'intpotCell=',intpotCell
            close (fileNo)
          
        end if
    end subroutine getIntegralPointInCell
    
    subroutine interpolateInCell (ynew,xnew,yold,xold)
        implicit none
        real*8,intent(in):: yold(:),xold(:)
        real*8,intent(out):: ynew(:),xnew(:)
        real*8,allocatable:: workArray(:)
        integer newsize,oldsize
        
        newsize = size (xnew)
        oldsize = size (xold)
        allocate( workArray(oldsize) )
!        write (*,*) 'xnew',xnew
!        write (*,*) 'yold',yold
!        write (*,*) 'xold',xold
!        write (*,*) 'newsize',newsize
!        write (*,*) 'oldsize',oldsize
        call spline1d(ynew,xnew,newsize,yold,xold,oldsize,workArray)
        
        deallocate (workArray)
        
    end subroutine interpolateInCell
    
    subroutine splintp(xa,ya,y2a,n,x,yp)
!-----------------------------------------------------------------------
!     cubic spline evaluator--spline must be called first to evaluate
!     y2a
!     ya is a function of xa--both are arrays of length n
!     ya2[1:n] contains spline coefficients calculated in spline
!     x is the argument of y[x] where y is to be evaluated
!     yp=y'[x] is the returned value
!-----------------------------------------------------------------------
      integer n
      real*8 xa(n),ya(n),y2a(n),x,yp
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
!      ap= -1.0/h
!      bp= -ap
!      
!      yp=ap*ya(klo)+bp*ya(khi)+ ( (3 * a**2 * ap - ap)*y2a(klo)+(3 * b**2 * bp - bp)*y2a(khi) )*(h**2)/6.
!     
      yp=(ya(khi)-ya(klo))/h-(3.*a**2-1.)/6.*h*y2a(klo)+(3.*b**2-1.)/6.*h*y2a(khi)
      return
    end subroutine splintp

    subroutine constructSpline(mesh,Func,splinepPara)
        implicit none
        
        real*8,intent(in):: mesh(:),Func(:) 
        real*8,intent(out):: splinepPara(:)
        real*8 yp1,ypn
        integer sizeMesh
        
        sizeMesh = size(mesh)
        yp1=-1.e30
        ypn=-1.e30
        !call spline(xold,yold,nold,yp1,ypn,y2old)
        call spline(mesh,Func,sizeMesh,yp1,ypn,splinepPara)
        
    end subroutine constructSpline
    
    function getSplineValue(mesh,Func,splinePara,x) result (val)
        implicit none
        real*8,intent(in)::mesh(:),Func(:),splinePara(:),x
        integer sizeMesh
        real*8 val
        sizeMesh = size(mesh)
        !call splint(xa,ya,y2a,n,x,y)
        call splint(mesh,Func,splinePara,sizeMesh,x,val)
    end function getSplineValue
    
    ! derivative of spline function    
    function getSplineDerValue(mesh,Func,splinePara,x) result (val)
        implicit none
        real*8,intent(in)::mesh(:),Func(:),splinePara(:),x
        integer sizeMesh
        real*8 val
        sizeMesh = size(mesh)
        !call splintp(xa,ya,y2a,n,x,yp)
        call splintp(mesh,Func,splinePara,sizeMesh,x,val)
        
    end function getSplineDerValue
    
    FUNCTION gammln(xx)
      IMPLICIT NONE
      REAL*8 gammln,xx
      INTEGER j
      REAL*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
    END FUNCTION gammln
    
    SUBROUTINE gaulag(x,w,n,alf)
      IMPLICIT NONE
      INTEGER n,MAXIT
      REAL*8:: alf,w(n),x(n)
      REAL*8 EPS
      PARAMETER (EPS=3.D-14,MAXIT=100)
      INTEGER i,its,j
      REAL*8 ai
      REAL*8 p1,p2,p3,pp,z,z1
      do 13 i=1,n
        if(i.eq.1)then
          z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
        else if(i.eq.2)then
          z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
        else
          ai=i-2
          z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))*(z-x(i-2))/(1.+.3*alf)
        endif
        do 12 its=1,MAXIT
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
11        continue
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        print *,'WARNING: too many iterations in gaulag'
1       x(i)=z
        w(i)=-exp(gammln(alf+n)-gammln(real(n,8)))/(pp*n*p2)
13    continue
      return
    END SUBROUTINE gaulag

    subroutine createEnergyMesh (mesh,a,b,option)

        
        implicit none

        real*8,intent(in):: a,b                             !Define the mesh interval b>a        
        integer meshPots
        real*8,intent(out)::mesh(:)
        
        real*8 RTMP,RLAM
        integer J
        integer,intent(in):: option
        
        meshPots = size (mesh)
        
        if (option .EQ. 1) then
!           Uniform mesh

            RTMP = (b-a)/meshPots
            do J=2,meshPots
                mesh(J) = a + RTMP*(J-1)
            enddo
            mesh(1) = a
            mesh(meshPots) = b
            
        else
        
!           Log10 mesh
            RLAM = 4.0
            RTMP = (LOG10(b-a)+RLAM)/(meshPots-2)
            do J=2,meshPots
                if (option .EQ. 2) then
                    mesh(J) = a + 10.**(-RLAM+(J-2)*RTMP)
                endif
                if (option .EQ. 3) then
                    mesh(meshPots + 1 - J) = b - 10.**(-RLAM+(J-2)*RTMP)
                endif 
            enddo
            mesh(1) = a 
            mesh(meshPots) = b
        endif
    end subroutine createEnergyMesh      

    SUBROUTINE ADAPTIVE_ENERGY_INTEGRAL_LSODE (INT_RESULT,KPARTICLE,INUTYPE,RNTOR,OMEGASN,OMEGAST,OMEGAE,OMEGA,RTMP,RTMP1,NUEFF)
        IMPLICIT NONE
        COMPLEX*16, INTENT(OUT) :: INT_RESULT
        INTEGER, INTENT(IN) :: KPARTICLE,INUTYPE
        REAL*8, INTENT(IN) :: RNTOR,NUEFF,RTMP,RTMP1
        REAL*8,INTENT(IN) :: OMEGASN,OMEGAST,OMEGAE
        COMPLEX*16,INTENT(IN) :: OMEGA
	    INTEGER :: JAC,MF,IOPT,ISTATE,ITASK,ITOL,LIW,LRW
        INTEGER :: NEQ,ISTEP
        PARAMETER (NEQ=2)
        INTEGER,ALLOCATABLE :: IWORK(:)
        REAL*8:: RTOL,ADTOL,E0,E1,TOUT
        REAL*8,ALLOCATABLE :: RWORK(:),ATOL(:)
        COMPLEX*16 :: Y(11)
!       INITIALIZE LSODE SOLOVER
        ALLOCATE (ATOL(1))
        ITOL = 1
        ATOL = 1E-14
        RTOL = 1E-15
        ADTOL = 1E-15
        MF   = 10
        SELECT CASE (MF)
        CASE (10)
            LIW = 20
            LRW = 20 + 16*NEQ
        CASE (22)
            LIW = 20 + NEQ
            LRW = 22 + 9*NEQ + NEQ**2
        END SELECT
        ITASK = 5
        IOPT  = 0
        JAC   = 0
        ALLOCATE (IWORK(LIW),RWORK(LRW))

!       INITIAL THE ENERGY INTEGRAL
        INT_RESULT = 0.0

        E0=0.0
        E1=50.0

        Y = 0.0
        Y(2) = KPARTICLE
        Y(3) = INUTYPE
        Y(4) = RNTOR
        Y(5) = OMEGASN
        Y(6) = OMEGAST
        Y(7) = OMEGAE
        Y(8) = OMEGA
        Y(9) = RTMP        
        Y(10) = RTMP1             ! RNTOR*(DU+DB)
        Y(11) = NUEFF
        DO 
            ISTEP=0
            Y(1)  = 0.0
            IWORK = 0.0
            RWORK = 0.0
            RWORK(1) = E1
            TOUT     = E1
            ISTATE = 1
!       ENERGY INTEGRAL BY LSODE
            DO
                CALL LSODE (THERMAL_PRECESSION_INT,NEQ,Y,E0,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
                IF ( ISTATE .LE. 0) THEN
                    STOP 'LSODE FAIL.'
                 ENDIF                
                 ISTEP = ISTEP +1
                 IF (ABS(E0-E1) .LE. 1E-15) EXIT 
            ENDDO

            INT_RESULT = INT_RESULT + Y(1)

            IF (ABS(Y(1)).LT. ADTOL) EXIT
            E1 = 1.5 * E1
            IF (E1.GT.1E4) STOP 'TOO LARGE ENERGY LIMITE. CHECK LSODE ENERGY INTEGRAL.'
        ENDDO

        IF (KPARTICLE .EQ. 0) INT_RESULT = INT_RESULT * 2.0

        DEALLOCATE (ATOL,IWORK,RWORK)

    END SUBROUTINE ADAPTIVE_ENERGY_INTEGRAL_LSODE

    SUBROUTINE THERMAL_PRECESSION_INT (NEQ,E,Y,YDOT)
        USE GLOBALM, ONLY: NTVNUMEREFAC
        INTEGER, INTENT(IN) :: NEQ
        REAL*8, INTENT(IN) :: E
        COMPLEX*16, INTENT(IN) :: Y(11)
        COMPLEX*16, INTENT(OUT) :: YDOT(1)

        COMPLEX*16 :: KPARTICLE,RNTOR,NUEFF,RTMP,RTMP1,OMEGASN,OMEGAST,OMEGAE,OMEGA
        COMPLEX*16 :: NUMER,DENOM,SQRTE,E15

        KPARTICLE   = Y(2)
        INUTYPE     = Y(3)
        RNTOR       = Y(4)
        OMEGASN     = Y(5)
        OMEGAST     = Y(6)
        OMEGAE      = Y(7)
        OMEGA       = Y(8)
        RTMP        = Y(9)
        RTMP1       = Y(10)
        NUEFF       = Y(11)
        SQRTE  = SQRT(E)
        IF ( ABS(INUTYPE + 1.0) .LT. 1E-15) THEN
           E15 = 1.0
        ELSE
           E15 = E**1.5
        ENDIF

        NUMER = RNTOR*(OMEGASN + (E - 1.5)*OMEGAST + &
            NTVNUMEREFAC*OMEGAE) - OMEGA
        IF (ABS(KPARTICLE) .LT. 1E-15) THEN
           DENOM = RNTOR*OMEGAE + RTMP1*E + RTMP*SQRTE - OMEGA        
           YDOT = NUMER / ( DENOM * E15 - NUEFF * CMPLX(0,1) )
        ELSE
           DENOM = RNTOR*OMEGAE + RTMP*SQRTE - OMEGA        
           YDOT = NUMER / ( DENOM * E15 - NUEFF * CMPLX(0,1) )
           DENOM = RNTOR*OMEGAE - RTMP*SQRTE - OMEGA         
           YDOT = YDOT + NUMER / ( DENOM * E15 - NUEFF * CMPLX(0,1) )
        ENDIF    

        IF ( ABS(INUTYPE + 1.0) .LT. 1E-15) THEN
           YDOT = YDOT * E**2.5 * EXP (-E)
        ELSE
           YDOT = YDOT * E**4 * EXP (-E)
        ENDIF

    END SUBROUTINE THERMAL_PRECESSION_INT

    SUBROUTINE RESIZE_ARRAY (A,ADDSIZE)
    REAL*8, ALLOCATABLE, INTENT(INOUT) :: A(:)
    INTEGER, INTENT(IN):: ADDSIZE
    INTEGER :: SIZEOFA
    REAL*8, ALLOCATABLE :: B(:)
    
    SIZEOFA = SIZE(A)
    ALLOCATE (B(SIZEOFA+ADDSIZE))
    B(1:SIZEOFA)=A
    DEALLOCATE(A)
    CALL MOVE_ALLOC(B,A)
    END SUBROUTINE RESIZE_ARRAY

    FUNCTION LININTERP(OLDX,OLDY,WORK,NEWX) RESULT (NEWY)
        IMPLICIT NONE
        REAL*8,INTENT(IN)::OLDX(:),OLDY(:),NEWX
        INTEGER,INTENT(INOUT):: WORK !ONE LINEAR INTERPRATION ONE WORK
        REAL*8:: NEWY
        INTEGER:: ISTEP,SIZEX

        SIZEX = SIZE(OLDX)

        IF ((WORK.LT.1) .OR. (WORK.GE.SIZEX)) THEN
           WORK=1
        ENDIF
 
        IF (NEWX.LT.OLDX(1)) THEN
           WORK=1
        ELSEIF (NEWX.GT.OLDX(SIZEX)) THEN
           WORK=SIZEX-1
        ELSE
           DO ISTEP=1,SIZEX-1
              IF( (NEWX.GE.OLDX(WORK)).AND.(NEWX .LE.OLDX(WORK+1)) ) THEN
                 EXIT
              ENDIF
              WORK=WORK+1
              IF (WORK .GE. SIZEX) THEN
                 WORK=1
              ENDIF
           ENDDO
           IF (ISTEP.GE.SIZEX) THEN
              STOP 'LININTERP ERROR: NEWX IS NOT IN OLDX'
           ENDIF
        ENDIF

        !INTERP VALUE AT NEWX
        NEWY = OLDX(WORK+1) - OLDX(WORK)
        NEWY = ( (NEWX-OLDX(WORK))*OLDY(WORK+1)+(OLDX(WORK+1)-NEWX)*OLDY(WORK) )/NEWY
        RETURN
    END FUNCTION LININTERP

    SUBROUTINE ADAPTIVE_GRID_LSODE (X,YIN)
        IMPLICIT NONE
        REAL*8, ALLOCATABLE, INTENT(INOUT) :: X(:)
        REAL*8, ALLOCATABLE, INTENT(IN) :: YIN(:,:)

        INTEGER :: SIZEX,SIZEY1,SIZEY2
        INTEGER :: JAC,MF,IOPT,ISTATE,ITASK,ITOL,LIW,LRW
        INTEGER :: NEQ,ISTEP,XLIM
        INTEGER,ALLOCATABLE :: IWORK(:)
        REAL*8:: RTOL,X0,X1,TOUT
        REAL*8,ALLOCATABLE :: RWORK(:),ATOL(:),X_LSODE(:),Y(:),IDX(:)
!       INITIALIZE LSODE SOLOVER
        
        SIZEX=SIZE(X)
        SIZEY1=SIZE(YIN,1)
        SIZEY2=SIZE(YIN,2)
        XLIM=2500

        ALLOCATE(intpwork(SIZEY1))
        ALLOCATE(x_adp(SIZEX),y_adp(SIZEY1,SIZEY2))
        ALLOCATE(sp_adp(SIZEY1,SIZEY2))
        x_adp=X
        y_adp=YIN
        ALLOCATE (Y(SIZEY1),ATOL(1),X_LSODE(XLIM))

        DO ISTEP=1,SIZEY1
            Y(ISTEP)=MAXVAL(ABS(YIN(ISTEP,:)))
            y_adp(ISTEP,:)=y_adp(ISTEP,:)/Y(ISTEP)
            CALL constructSpline(x_adp,y_adp(ISTEP,:),sp_adp(ISTEP,:))
        ENDDO
        intpwork=1


        NEQ=SIZEY1
        ITOL = 1
        ATOL = 1E-10
        RTOL = 1E-12
        MF   = 10
        SELECT CASE (MF)
        CASE (10)
            LIW = 20
            LRW = 20 + 16*NEQ
        CASE (22)
            LIW = 20 + NEQ
            LRW = 22 + 9*NEQ + NEQ**2
        END SELECT
        ITASK = 5
        IOPT  = 0
        JAC   = 0
        ALLOCATE (IWORK(LIW),RWORK(LRW))

        X0=X(1)
        X1=X(SIZEX)

        Y = 0.0
        
        ISTEP=1
        IWORK = 0.0
        RWORK = 0.0
        RWORK(1) = X1
        TOUT     = X1
        ISTATE = 1
        DO
            CALL LSODE (INT_QUANTITY,NEQ,Y,X0,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
            IF ( ISTATE .LE. 0 .OR. ISTEP .GT. 100000) THEN
                STOP 'LSODE FAIL IN ADAPTIVE_GRID_LSODE.'
            ENDIF
            X_LSODE(ISTEP)=X0                
            ISTEP = ISTEP +1
            IF (ISTEP.GT.XLIM) THEN
                CALL RESIZE_ARRAY (X_LSODE,2500)
                XLIM=SIZE(X_LSODE)
            ENDIF
            IF (ABS(X0-X1) .LE. 1E-20) EXIT 
        ENDDO

        ALLOCATE(IDX(ISTEP-1))
        IDX=(/(REAL(XLIM), XLIM=1,ISTEP-1)/)
        IDX=(IDX-1)*(SIZEX-1)/(ISTEP-2)+1
        NEQ=1
        DO XLIM=2,SIZEX-1
            X(XLIM) = LININTERP(IDX,X_LSODE,NEQ,REAL(XLIM))
        ENDDO
        

        DEALLOCATE (ATOL,IWORK,RWORK,Y,X_LSODE,IDX)
        DEALLOCATE (intpwork, x_adp,y_adp,sp_adp)
        RETURN
    END SUBROUTINE ADAPTIVE_GRID_LSODE

    SUBROUTINE INT_QUANTITY (NEQ,X,Y,YDOT)
        INTEGER, INTENT(IN) :: NEQ
        REAL*8, INTENT(IN) :: X
        REAL*8, DIMENSION(NEQ), INTENT(IN) :: Y
        REAL*8, DIMENSION(NEQ), INTENT(OUT) :: YDOT

        INTEGER :: ICOMP

        DO ICOMP=1,NEQ
!           YDOT(ICOMP) = LININTERP(x_adp,y_adp(ICOMP,:),intpwork(ICOMP),X)
!           YDOT(ICOMP) = getSplineDerValue(x_adp,y_adp(ICOMP,:),sp_adp(ICOMP,:),x)
           YDOT(ICOMP) = getSplineValue(x_adp,y_adp(ICOMP,:),sp_adp(ICOMP,:),x)
        ENDDO
    END SUBROUTINE INT_QUANTITY

end module ToolBox
