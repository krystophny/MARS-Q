module KineticExtendAnalysis
    type SurfaceFractionDistribution
        real*8,allocatable:: postionKai(:,:)                         
        real*8,allocatable:: fractionDistribution(:,:)
        real*8 averfraction
    
    end type SurfaceFractionDistribution
        
    contains
       
    
    ! function body
    subroutine createSurfaceFractionDistribution (SFD,arraySize)
        implicit none
        integer,intent(in):: arraySize
        
        type(SurfaceFractionDistribution),pointer:: SFD
        call releaseSurfaceFractionDistribution ( SFD )
        allocate ( SFD )   
        allocate( SFD%postionKai( arraySize,2 ) )
        allocate( SFD%fractionDistribution( arraySize,2 ) )
        
    end subroutine createSurfaceFractionDistribution
    subroutine releaseSurfaceFractionDistribution(SFD)
        implicit none
        
        type(SurfaceFractionDistribution),pointer:: SFD
        if ( associated(SFD) ) then
            if ( allocated(SFD%postionKai) ) deallocate ( SFD%postionKai )
            if ( allocated(SFD%fractionDistribution) ) deallocate ( SFD%fractionDistribution )
            deallocate( SFD )
        end if
        
    end subroutine releaseSurfaceFractionDistribution 

!   calculate SurfacefractionDistribution    
    subroutine calcSurfaceFractionDistribution(kCheck,SFD,energyPartFNum,energyPartFDen,meshCw, Cw, oldsurfJac, minHkNum, oldHkNum, minHkDen, oldHkDen,oldHk)
        
        ! SFD for returning the particle distribtuion fraction 
        ! energyPartFNum energyPartFDen, the parts (numerator and denominator )relate to energy in the integration ( normally we can calculate analitically)
        ! meshCw(:), Cw(:) the mesh and the corresponding value for particle frequency
        ! oldsurfJac(:) Jacobian on the surface
        ! minHkNum, oldHkNum(:), minHkDen, oldHkDen(:): the integratal range for numerator and denominator along Pitch direction
        ! oldHk the Hk value used in the integration
        use ToolBox
        implicit none
        
        integer kCheck,fileNo
        
        real*8, parameter:: pi = 3.14159265358979
        
        type(SurfaceFractionDistribution),pointer,intent(out):: SFD
        
        real*8,intent(in):: energyPartFNum,energyPartFDen
        
       
        ! pitch angle direction
        real*8,intent(in):: meshCw(:),Cw(:)
        real*8,allocatable::spCw(:)  ! spline parameter for Cw
       
        ! kai direction 
        real*8,intent(in):: oldsurfJac(:)
        real*8,intent(in):: oldHkNum(:),oldHkDen(:),oldHk(:)
        real*8,intent(in):: minHkNum,minHkDen
        
        
        real*8,allocatable:: meshKai(:), intpotKaiCell(:,:)
        ! intStepPitch the step length of integration along Pitch direction
        ! intStepKai the step length of integration along Kai direction
        real*8 intStepPitch, intStepKai
        integer sizeKaiMesh,kaiCellSize
        
        integer i
        
        ! the value of Jacobian and Hk on the Kai mesh node
        real*8,allocatable:: meshsurfJac(:),meshHkNum(:),meshHkDen(:),meshHk(:)
        real*8,allocatable:: surfJacCell(:,:),hkNumCell(:,:),hkDenCell(:,:),hkCell(:,:)
              
        ! integration along pitch angle for numerator and denominator
        real*8,allocatable:: intNumPitch(:,:), intDenPitch(:,:)
        
        real*8 hkDen,hkNum,hk
        integer cellNo
        
        real*8 tmp1,tmp2,tmp3
        
        real*8 Cwhkdown,Cwhkup
        
        real*8,allocatable:: meshPitch(:),intpotPitchCell(:,:)
        
        real*8 intpitch1,intpitch2,dpitch,Cwprime1,Cwprime2
        integer pitchCellNo
        real*8 intPitch  ! integration along Pitch direction
        
        real*8 intKaiNum,intKaiDen ! integration along Kai direction
        real*8 dKai
        
        if (kCheck == 1) then
            fileNo = assignFreeFileUnit()
            open (unit=fileNo, file='SFD_DEBUG.OUT',STATUS='unknown')
        end if        
        ! calculate Kai mesh 
        i = size(meshCw)
        intStepPitch = meshCw(i) / (i+1)
        if ( kCheck == 1) then
            write (fileNo,*) 'intStepPitch=',intStepPitch
            write (fileNo,*) 'meshCw(i)=',meshCw(i)
        end if        
        sizeKaiMesh = size(oldsurfJac) + 1
        allocate( meshKai(sizeKaiMesh) )
        intStepKai = 2.0 * pi / ( sizeKaiMesh - 1 )
        
        do i=1,sizeKaiMesh
            meshKai(i) = intStepKai * ( i - 1 )
        end do        
        
        call getIntegralPointInCell (0,meshKai,intpotKaiCell)
        
        allocate( meshsurfJac(sizeKaiMesh),meshHkNum(sizeKaiMesh),meshHkDen(sizeKaiMesh),meshHk(sizeKaiMesh) )
        
        meshsurfJac(1 : sizeKaiMesh - 1) = oldsurfJac
        meshHkNum(1 : sizeKaiMesh - 1) = oldHkNum
        meshHkDen(1 : sizeKaiMesh - 1) = oldHkDen
        meshHk(1 : sizeKaiMesh - 1) = oldHk
        
        meshsurfJac(sizeKaiMesh) = oldsurfJac(1)
        meshHkNum(sizeKaiMesh) = oldHkNum(1)
        meshHkDen(sizeKaiMesh) = oldHkDen(1)
        meshHk(sizeKaiMesh) = oldHk(1)
        
        kaiCellSize = size(intpotKaiCell,dim = 1)
        allocate ( surfJacCell(kaiCellSize,2),hkNumCell(kaiCellSize,2),hkDenCell(kaiCellSize,2),hkCell(kaiCellSize,2) )
        
        do i=1,2
            !interpolateInCell (ynew,xnew,yold,xold)         
            call interpolateInCell (surfJacCell(:,i),intpotKaiCell(:,i),meshsurfJac,meshKai)
            call interpolateInCell (hkNumCell(:,i),intpotKaiCell(:,i),meshHkNum,meshKai)
            call interpolateInCell (hkDenCell(:,i),intpotKaiCell(:,i),meshHkDen,meshKai)
            call interpolateInCell (hkCell(:,i),intpotKaiCell(:,i),meshHk,meshKai)
        end do
        
        if ( kCheck == 1) then
            write (fileNo,*) 'kaiCellSize=',kaiCellSize
            write (fileNo,*) 'hkCell(:,1)=',hkCell(:,1)  
        end if
        
        allocate (intNumPitch(kaiCellSize,2), intDenPitch(kaiCellSize,2))
        ! Denominator calculation before integration along Kai direction
        allocate ( spCw( size(meshCw) ) )
        call constructSpline(meshCw,Cw,spCw)
        if ( kCheck == 1) then
            write (fileNo,*) 'spCw=',spCw
        end if              
        intKaiNum=0.0
        intKaiDen=0.0
        
        ! calculate the numerator and denominator in each cell
        do cellNo = 1, kaiCellSize
            do i=1,2
                
                hk = hkCell(cellNo,i)
                !denominator
                hkDen = hkDenCell(cellNo,i) 
                
                if (minHkDen + epsilon > hkDen) then
                    intDenPitch(cellNo,i) = 0.0
                else

                    tmp1 = 1.0 - hkDen / hk
                    tmp2 = 1.0 - minHkDen / hk
                    if ( tmp1 .le. 0 ) tmp1 = 0.0
                    if ( tmp2 .le. 0 ) tmp2 = 0.0
                    tmp1 =  dsqrt (tmp1)     
                    tmp2 =  dsqrt (tmp2)     
                    
                    intDenPitch(cellNo,i) = -2.0 * ( tmp1 - tmp2 )
                end if
                
                !numerator 
                hkNum = hkNumCell(cellNo,i)
                if ( minHkNum + epsilon > hkNum ) then
                    intNumPitch(cellNo,i) = 0.0
                    goto 1
                end if
                
                Cwhkup = getSplineValue (meshCw,Cw,spCw,hkNum) 
                Cwhkdown = getSplineValue (meshCw,Cw,spCw,minHkNum)
                tmp1 = 1.0 - hkNum / hk
                tmp2 = 1.0 - minHkNum / hk
                if ( tmp1 .le. 0 ) tmp1 = 0.0
                if ( tmp2 .le. 0 ) tmp2 = 0.0
                tmp1 = Cwhkup * dsqrt (tmp1)
                tmp2 = Cwhkdown * dsqrt (tmp2)
                tmp1 = -2.0 * ( tmp1 - tmp2 )
                
                intPitch = 0.0
                if ( kCheck == 1) then
                        write (fileNo,*) 'minHkNum,hkNum,intStepPitch=',minHkNum,hkNum,intStepPitch                        
                end if                  
                call createMesh (minHkNum,hkNum,intStepPitch,meshPitch)
                call getIntegralPointInCell (kCheck,meshPitch,intpotPitchCell)
                
                intPitch = 0.0
                if ( kCheck == 1) then
                    write (fileNo,*) 'meshPitch=',meshPitch
                    write (fileNo,*) 'intpotPitchCell=',intpotPitchCell    
                    write (fileNo,*) 'size(intpotPitchCell,dim = 1)',size(intpotPitchCell,dim = 1)
                end if   
                do pitchCellNo = 1, size(intpotPitchCell,dim = 1)
                    
                    intpitch1 = intpotPitchCell (pitchCellNo,1)
                    intpitch2 = intpotPitchCell (pitchCellNo,2)
                    dpitch = meshPitch (pitchCellNo+1) - meshPitch (pitchCellNo)
                    Cwprime1 = getSplineDerValue(meshCw,Cw,spCw,intpitch1)
                    Cwprime2 = getSplineDerValue(meshCw,Cw,spCw,intpitch2)
                    if ( kCheck == 1) then
                        write (fileNo,*) 'Cwprime1=',Cwprime1
                        write (fileNo,*) 'Cwprime2=',Cwprime2 
                    end if                    
                    tmp2 = dsqrt( 1.0 - intpitch1 / hk ) 
                    tmp3 = dsqrt( 1.0 - intpitch2 / hk )
                    
                    tmp2 = ( Cwprime1 * tmp2 + Cwprime2 * tmp3 ) * dpitch ! 1/2 is canncled by 2 from formula
                    
                    
                    intPitch = intPitch + tmp2
                end do
                    intNumPitch(cellNo,i) = tmp1 + intPitch
                
                deallocate (meshPitch,intpotPitchCell)
               
1            end do
            
            dKai = meshKai (cellNo + 1) - meshKai (cellNo)
            intKaiNum = intKaiNum + 0.5 * ( intNumPitch(cellNo,1) * surfJacCell(cellNo,1) + intNumPitch(cellNo,2) * surfJacCell(cellNo,2) ) * dKai
            intKaiDen = intKaiDen + 0.5 * ( intDenPitch(cellNo,1) * surfJacCell(cellNo,1) + intDenPitch(cellNo,2) * surfJacCell(cellNo,2) ) * dKai
            if ( kCheck == 1) then
                write (fileNo,*) 'intKaiNum=',intKaiNum
                write (fileNo,*) 'intKaiDen=',intKaiDen 
            end if   
             
        end do
        
        call createSurfaceFractionDistribution (SFD,kaiCellSize)
        
        tmp1 = energyPartFDen * intKaiDen
        if ( tmp1 < epsilon ) then
            SFD%averfraction = 0.0
        else 
            SFD%averfraction = energyPartFNum * intKaiNum / tmp1
        end if
        SFD%postionKai = intpotKaiCell
        
        do cellNo = 1, kaiCellSize
            do i=1,2
                tmp1 =  energyPartFDen * intDenPitch(cellNo,i)
                if (tmp1 < epsilon) then
                    SFD%fractionDistribution(cellNo,i)  = 0.0
                else
                    SFD%fractionDistribution(cellNo,i) = energyPartFNum * intNumPitch(cellNo,i) / tmp1
                end if
            end do
        end do 
        
        if ( kCheck == 1) then
            close (fileNo) 
        end if 
        
         
        deallocate ( spCw, meshKai, intpotKaiCell )
        deallocate ( surfJacCell,hkNumCell,hkDenCell,hkCell )
        deallocate ( meshsurfJac,meshHkNum,meshHkDen,meshHk )
        deallocate ( intNumPitch, intDenPitch )
     end subroutine calcSurfaceFractionDistribution


end module KineticExtendAnalysis

