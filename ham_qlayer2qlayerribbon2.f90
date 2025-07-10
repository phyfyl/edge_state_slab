subroutine ham_qlayer2qlayerribbon2(k,H00,H01)
     use para

     implicit none

     ! loop index
     integer :: iR, i1,i2,j1,j2

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ia,inew_ib

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k

     complex(dp) :: ratio

     complex(Dp) :: Hij(-ijmax:ijmax,-ijmax:ijmax,Num_wann,Num_wann)

     complex(Dp), intent(out) :: H00(nslab*Ndim, nslab*Ndim), H01(nslab*Ndim, nslab*Ndim)
  
     complex(Dp) :: hm(nslab*Ndim,nslab*Ndim), he(nslab*Ndim,nslab*Ndim)
     Hij=0.0d0
     do iR=1,Nrpts
        ia=irvec(1,iR)
        ib=irvec(2,iR)
        ic=irvec(3,iR)

        !> new lattice
        call latticetransform(ia, ib, ic, new_ia, new_ib, new_ic)
      !  call rotation(ia, ib, ic, new_ia, new_ib, new_ic)

        inew_ia= int(new_ia)
        inew_ib= int(new_ib)
        if (abs(new_ia).le.ijmax)then
        if (abs(new_ib).le.ijmax)then
           kdotr=k*new_ic
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Hij(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           =Hij(inew_ia, inew_ib, 1:Num_wann, 1:Num_wann )&
           +HmnR(:,:,iR)*ratio/ndegen(iR)
        endif
        endif

     enddo
     H00=0d0
     H01=0d0
     ! i1,j1 row index
     do i1=1,Np
     do j1=1,nslab
     ! i2,j2 column index
     do i2=1,Np
     do j2=1,nslab
       if (abs(i2-i1).le.ijmax.and.abs(j2-j1).le.ijmax)then
         H00((i1-1)*nslab*Num_wann+(j1-1)*Num_wann+1: &
             (i1-1)*nslab*Num_wann+(j1-1)*Num_wann+Num_wann,&
             (i2-1)*nslab*Num_wann+(j2-1)*Num_wann+1: &
             (i2-1)*nslab*Num_wann+(j2-1)*Num_wann+Num_wann )&
         = Hij(i2-i1,j2-j1,1:Num_wann,1:Num_wann)

       endif
     enddo
     enddo
     enddo
     enddo
     
     He = 0d0
     ! i1,j1 row index
     do i1=1,Np
     do j1=1,nslab
        do i2=1, Num_wann
         He((i1-1)*nslab*Num_wann+(j1-1)*Num_wann+i2, &
           (i1-1)*nslab*Num_wann+(j1-1)*Num_wann+i2 )&
         =E*(-nslab*Origin_cell%cell_parameters(E_direction)/2+&
               (j1-1)*Origin_cell%cell_parameters(E_direction)+&
               Origin_cell%wannier_centers_cart(E_direction,i2))*&
               eV2Hartree/Angstrom2atomic
      enddo
      enddo
      enddo
     
     hm=0.0d0
  !   do i1=1, Np
  !   do i2=1, Num_wann/2
  !       hm((i1-1)*nslab*Num_wann+i2,(i1-1)*nslab*Num_wann+i2)&
  !         =mz*eV2Hartree
  !       hm((i1-1)*nslab*Num_wann+i2+Num_wann/2,&
  !          (i1-1)*nslab*Num_wann+i2+Num_wann/2)&
  !         =-mz*eV2Hartree
  !       hm((i1-1)*nslab*Num_wann+(nslab-1)*Num_wann+i2,&
  !          (i1-1)*nslab*Num_wann+(nslab-1)*Num_wann+i2)&
  !         =-mz*eV2Hartree
  !       hm((i1-1)*nslab*Num_wann+(nslab-1)*Num_wann+i2+Num_wann/2, &
  !          (i1-1)*nslab*Num_wann+(nslab-1)*Num_wann+i2+Num_wann/2)&
  !         =mz*eV2Hartree
  !   enddo
   !  enddo

     H00=H00+hm+he

     ! i1,j1 row index
     do i1=1,Np
     do j1=1,nslab
     ! i2,j2 column index
     do i2=Np+1,Np*2
     do j2=1,nslab
       if (abs(i2-i1).le.ijmax.and.abs(j2-j1).le.ijmax)then
         H01((i1-1)*nslab*Num_wann+(j1-1)*Num_wann+1: &
             (i1-1)*nslab*Num_wann+(j1-1)*Num_wann+Num_wann,&
             (i2-1-Np)*nslab*Num_wann+(j2-1)*Num_wann+1: &
             (i2-1-Np)*nslab*Num_wann+(j2-1)*Num_wann+Num_wann )&
         = Hij(i2-i1,j2-j1,1:Num_wann,1:Num_wann)

       endif
     enddo
     enddo
     enddo
     enddo

     return
  end subroutine ham_qlayer2qlayerribbon2
