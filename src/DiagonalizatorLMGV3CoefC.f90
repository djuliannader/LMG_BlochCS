!******************
! compile as:
! gfortran -o LMGcoefC.exe ./src/DiagonalizatorLMGV3CoefC.f90 -llapack -lblas
!************************
!------------------------------------------------------------------------!
! Author Daniel Julian Nader                                             !
! Collaboration Sergio Lerma                                             !
!                                                                        !
!  This code:                                                            !          
!     -  diagonalize LGM Model in SU(2) representation                   !
!     -  Calculate  the profile of the Bloch                             !
!        coherent state |z(theta,phi)>                                   !
!        by projecting the eigenstate ck^2=|<E_k|z>|^2                   !                                           
!                                                                        !
! The output for each parity is printed in the file: - resultsp+ck.dat   !
!                                                    - resultsp-ck.dat   !
!                                                                        !
! Gnuplot visualization:                                                 !
!                     gnulplot>plot 'resultsp+ck.dat' u 2:3 w lp         !
!************************************************************************!
!---------------!
 program LGMmodel
!---------------!

   !
   !********* variables definition  ******************************
   !
   
 implicit none

 integer k,nn

 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(16)

 integer :: i,ii,j,jq,jj,kk,kkk,n,nMC,m,mq
 real(8), allocatable :: m0p(:,:),m1p(:,:),m2p(:,:),eigp(:)
 real(8), allocatable :: m0n(:,:),m1n(:,:),m2n(:,:),eign(:)
 real(8) :: epsilon,uMC, w, v, gammax, gammay  ! Hamiltonian
 real(DP) :: alabs,alere,aleim,theta,phi  ! Husimi
 real(DP) :: Q2temp,QM2,Q2sum,Qpsitemp,Qpsi,Wtemp,Wpsi  ! Husimi

 real(DP) :: Qh,u,vv,delta,bin,pi

 pi=3.14159265359
 


  !**************  Read data  ****************
 
 print*,"Angular Momentum J"
 read*, jj
 print*,"gap energy epsilon (parameter)"
 read*, epsilon
 print*,"gammax (parameter)"
 read*, gammax
 print*,"gammay (parameter)"
 read*, gammay
 print*,'Theta (initial coherent state)'
 read*,theta
 print*,'Phi   (initial coherent state)'
 read*,phi

 v=epsilon*(gammax-gammay)/(2.0*(2.0*jj-1.0))
 w=epsilon*(gammax+gammay)/(2.0*(2.0*jj-1.0))

  ! ************   Diagonalization       ****************

 



 
 allocate (m0p(jj+1,jj+1))
 allocate (m1p(jj+1,jj+1))
 allocate (m2p(jj+1,jj+1))
 allocate (eigp(jj+1))
 allocate (m0n(jj,jj))
 allocate (m1n(jj,jj))
 allocate (m2n(jj,jj))
 allocate (eign(jj))



 
!!$ !********************************************************************
!!$ ! ************   Calculate matrix elements**********************
!!$ !********************************************************************
 ! All elements starting as zero
do i=1,jj
do j=1,jj
   m0p(i,j)=0.0
 enddo
enddo

do i=1,jj-1
do j=1,jj-1
   m0p(i,j)=0.0
 enddo
enddo


 !!**** Elements in the diagonal**********!

 do i=0,jj
    m0p(i+1,i+1)=epsilon*(-jj+2*i)+(w)*(jj*(jj+1)-(-jj+2*i)**2)    ! <------   !!!positive parity
 enddo

 do i=0,jj-1
    m0n(i+1,i+1)=epsilon*(-jj+2*i+1)+(w)*(jj*(jj+1)-(-jj+2*i+1)**2)    ! <------   !!! negative parity
 enddo



 !!**** Elements out of the diagonal**********!

    do i=0,jj-1
       m0p(i+1,i+2)=1.0*(v/2.0)*(jj*(jj+1)-(-jj+2*i+2)*((-jj+2*i+2)-1))**(1.0/2.0)* &
           & (jj*(jj+1)-((-jj+2*i+2)-1)*((-jj+2*i+2)-2))**(1.0/2.0)  ! <------   !!! positive parity            
    enddo

    do i=0,jj-2
       m0n(i+1,i+2)=1.0*(v/2.0)*(jj*(jj+1)-(-jj+2*i+3)*((-jj+2*i+3)-1))**(1.0/2.0)* &
           & (jj*(jj+1)-((-jj+2*i+3)-1)*((-jj+2*i+3)-2))**(1.0/2.0) ! <------   !!! negative parity
    enddo




 do i=0,jj-1
    m0p(i+2,i+1)=m0p(i+1,i+2)         ! <------   !!! positive parity
 enddo

 do i=0,jj-2
    m0n(i+2,i+1)=m0n(i+1,i+2)         ! <------   !!! negative parity
 enddo

  
 !   

 m1p(:,:)=m0p(:,:)
 m1n(:,:)=m0n(:,:)

! Call routini which will perform diagonalization (positive/negative parity separated)
 call diasym(m1p,eigp,jj+1)
 call diasym(m1n,eign,jj)

 




!******************************************************************************
 !                           Printing results
!******************************************************************************


!***parameters*******!
print*,'Diagonalization'
print*, 'Angular momentum= ', jj
print*,'Parameters'
print*,'epsilon= ',epsilon
print*,'gammax= ',gammax
print*,'gammay= ',gammay



!***Energias******!

 print*,'Eigenvalues (positive parity):'
 do i=1,jj+1
    !if (eig(i)< 0.0) then
     print*,i,eigp(i)
!    10 format(I3,'   ',f14.8)
     !endif
  enddo

 print*,'Eigenvalues (negative parity):'
 do i=1,jj
    !if (eig(i)< 0.0) then
     print*,i,eign(i)
!    10 format(I3,'   ',f14.8)
     !endif
  enddo 

   print*,'************************'
   print*,'Husimi Representation   *'
   print*,'see files:              *'
   print*,'        resultsp+ck.dat *'
   print*,'        resultsp-ck.dat *'
   print*,'************************'


!**** eigenvectors***!

!!$    print*,'Eigenvector:'
!!$    do i=1,jj+1
!!$    print*,i,m1p(i,:)
!!$     20 format(i3,'   ',10f14.8)
!!$    enddo
!!$ print*,

!******************************************************************************
!***** Printing results in the output file*********!
!!$!******************************************************************************
  open(1,file='resultsp+ck.dat',status='new')

  write(1,*)'# ********Parameters**********  #'
  write(1,*)'# angular momentum',jj
  write(1,*)'# Epsilon',epsilon
  write(1,*)'# gammax',gammax
  write(1,*)'# gammay',gammay
  write(1,*)'# Epsilon',epsilon
  write(1,*)'# Theta (initial coherent state)',theta
  write(1,*)'# Phi   (initial coherent state)',phi
  write(1,*)'# positive parity '
  write(1,*)'# ****************************  #'
  write(1,*)'#   k                  E_k/J                   c_k^2    #'





  alabs=1.0d0*tan(theta/2.0)
  do ii=1,jj+1
  !***Husimi************!
   alere=0.0d0
   aleim=0.0d0
   do m=0,jj
   call binomial(2*jj,jj+(-jj+2*m),bin)
   alere=alere+m1p(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m)))*cos((jj+(-jj+2*m))*phi)
   aleim=aleim-m1p(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m)))*sin((jj+(-jj+2*m))*phi)
   enddo
   Qh= 1.0d0*(alere**2.0+aleim**2.0)
   write(1,*)ii,' ' ,eigp(ii)/jj, '  ',Qh 

   enddo
   
   close(1)



   !********* Negative parity ***********************!

     open(1,file='resultsp-ck.dat',status='new')

  write(1,*)'# ********Parameters**********  #'
  write(1,*)'# angular momentum',jj
  write(1,*)'# Epsilon',epsilon
  write(1,*)'# gammax',gammax
  write(1,*)'# gammay',gammay
  write(1,*)'# Epsilon',epsilon
  write(1,*)'# Theta (initial coherent state)',theta
  write(1,*)'# Phi   (initial coherent state)',phi
  write(1,*)'# negative parity '
  write(1,*)'# ****************************  #'
  write(1,*)'#   k                  E_k/J                   c_k^2    #'


  alabs=1.0d0*tan(theta/2.0)
  do ii=1,jj+1
  !***Husimi************!
   alere=0.0d0
   aleim=0.0d0
   do m=0,jj-1
   call binomial(2*jj,jj+(-jj+2*m+1),bin)
   alere=alere+m1n(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m+1)))*cos((jj+(-jj+2*m+1))*phi)
   aleim=aleim-m1n(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m+1)))*sin((jj+(-jj+2*m+1))*phi)
   enddo
   Qh= 1.0d0*(alere**2.0+aleim**2.0)
   write(1,*)ii,' ' ,eign(ii)/jj, '  ',Qh 

   enddo
   
   close(1)




  
  m2p=matmul(transpose(m1p),m0p)
  m0p=matmul(m2p,m1p)









  deallocate(m0p); deallocate(m1p); deallocate(m2p);  deallocate(eigp)

 end program LGMmodel


!-------------------!

!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
 subroutine diasym(a,eig,n)
 implicit none

 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

 end subroutine diasym
!---------------------!


!********   Subroutine for the calculation of binomial factors ****!
 
subroutine binomial(x,y,b)
     implicit none

     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(16)
     integer :: x, y, i
     REAL(DP) :: b, binnum, binden1, binden2

     binnum=1.0d0
     do i=0,X-1
     binnum=binnum*(X-i)
     enddo
     binden1=1.0d0
     do i=0,Y-1
     binden1=binden1*(Y-i)
     enddo
     binden2=1.0
     do i=0,X-Y-1
     binden2=binden2*(X-Y-i)
     enddo
     b=1.0d0*(binnum/(binden1*binden2))

     end subroutine binomial
