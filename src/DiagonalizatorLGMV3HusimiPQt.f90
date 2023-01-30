!******************
! compile as:
! gfortran -o LGMHusimiPQt.exe ./DiagonalizatorLMGV3HusimiPQt.f90 -llapack -lblas
!************************
!------------------------------------------------------------------------!
! Author Daniel Julian Nader                                             !
! Collaboration Sergio Lerma                                             !
!                                                                        !
!  This code:                                                            !        
!     -  diagonalize LGM Model in SU(2) representation                   !
!     -  Calculate  the time evolution of the Husimi function            !
!        Q(P,Q,t)=|<z(alpha0,phi0,t=0)|z(P,Q,t)>|                        !
!                                                                        !
! The output is printed in the file:   - HusimiPQt.dat                   !
!                                                                        !
! Gnuplot visualization:                                                 !
!                     gnulplot>splot 'HusimiPQt.dat' u 1:2:3 w l         !
!************************************************************************!

!---------------!
 program LGMmodel
!---------------!

   !
   !********* variables definition ******************************
   !
   
 implicit none

 integer k,nn

 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(16)

 integer :: i,ii,j,jq,jj,kk,n,nMC,m,mq
 real(8), allocatable :: m0p(:,:),m1p(:,:),m2p(:,:),eigp(:)
 real(8), allocatable :: m0n(:,:),m1n(:,:),m2n(:,:),eign(:)
 real(8), allocatable :: Q(:),ThetaL(:),PhiL(:)
 real(8) :: epsilon,uMC, Ps, w, v, gammax, gammay  ! Hamiltonian
 real(DP) :: bkire,bkiim,bkre,bkim,absbksqrt,thetai,phii,theta,phi,alabs  
 real(DP) :: sumre,sumim,t,alabsi
 real(DP) :: Q2temp,QM2,Q2sum,Qpsitemp,Qpsi,Wtemp,Wpsi  ! Husimi function

 real(DP) :: Qh,u,vv,delta,bin,pi,test
 real(DP) :: QQ,PP,jz

 pi=3.14159265359
 


  !**************  Reading data  ****************
 
 print*,"Angular Momentum J"
 read*, jj
 print*,"gap energy epsilon (parameter)"
 read*, epsilon
 print*,"gammax (parameter)"
 read*, gammax
 print*,"gammay (parameter)"
 read*, gammay
 print*,'initial theta (coherent state)'
 read*,thetai
 print*,'initial phi (coherent state)'
 read*,phii
 print*,'time'
 read*,t
 
 v=epsilon*(gammax-gammay)/(2.0*(2.0*jj-1.0))
 w=epsilon*(gammax+gammay)/(2.0*(2.0*jj-1.0))

  ! ************   Diagonalization       ****************

 
! Defining the size of the matrix


 
 allocate (m0p(jj+1,jj+1))
 allocate (m1p(jj+1,jj+1))
 allocate (m2p(jj+1,jj+1))
 allocate (eigp(jj+1))
 allocate (m0n(jj,jj))
 allocate (m1n(jj,jj))
 allocate (m2n(jj,jj))
 allocate (eign(jj))
 allocate (Q(30000))




 
!!$ !********************************************************************
!!$ ! ************  Calculating the matriz elements**********************
!!$ !********************************************************************
 ! starting all elements in zero
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


 !!**** Elements on the diagonal**********!

 do i=0,jj
    m0p(i+1,i+1)=epsilon*(-jj+2*i)+(w)*(jj*(jj+1)-(-jj+2*i)**2)    ! <------   !!!paridad positiva
 enddo

 do i=0,jj-1
    m0n(i+1,i+1)=epsilon*(-jj+2*i+1)+(w)*(jj*(jj+1)-(-jj+2*i+1)**2)    ! <------   !!!paridad negativa
 enddo



 !!**** Elements out of the diagonal **********!

    do i=0,jj-1
       m0p(i+1,i+2)=1.0*(v/2.0)*(jj*(jj+1)-(-jj+2*i+2)*((-jj+2*i+2)-1))**(1.0/2.0)* &
           & (jj*(jj+1)-((-jj+2*i+2)-1)*((-jj+2*i+2)-2))**(1.0/2.0)  ! <------   !!!positive parity
    enddo

    do i=0,jj-2
       m0n(i+1,i+2)=1.0*(v/2.0)*(jj*(jj+1)-(-jj+2*i+3)*((-jj+2*i+3)-1))**(1.0/2.0)* &
           & (jj*(jj+1)-((-jj+2*i+3)-1)*((-jj+2*i+3)-2))**(1.0/2.0) ! <------   !!!negative parity
    enddo




 do i=0,jj-1
    m0p(i+2,i+1)=m0p(i+1,i+2)         ! <------   !!!positive parity
 enddo

 do i=0,jj-2
    m0n(i+2,i+1)=m0n(i+1,i+2)         ! <------   !!!negative parity
 enddo

  

 m1p(:,:)=m0p(:,:)
 m1n(:,:)=m0n(:,:)

! Calling routine which perform diagonalization
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



!***Energies******!

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

   print*,'***************************'
   print*,'Husimi Representation     *'
   print*,'see file:                 *'
   print*,'         HusimitQPt.dat   *'
   print*,'                          *'
   print*,'**************************'



!******************************************************************************
!***** Calculating the time evolution of the Husimi function         *********!
!!$!******************************************************************************
  open(1,file='HusimiPQt.dat',status='new')

!  t=0.0
  write(1,*)'# ********Parameters**********  #'
  write(1,*)'# angular momentum',jj
  write(1,*)'# Epsilon',epsilon
  write(1,*)'# gammax',gammax
  write(1,*)'# gammay',gammay
  write(1,*)'# Epsilon',epsilon
  write(1,*)'# Initial Theta (coherent state)',thetai
  write(1,*)'# Initial Phi   (coherent state)',phii
  write(1,*)'# positive parity '
  write(1,*)'# time ',t
  write(1,*)'# ****************************  #'
  write(1,*)'#     u             v                Q   '       

  nn=50
  delta=pi/(nn*1.0d0)
  theta=0.0d0
  alabsi=1.0d0*tan(thetai/2.0)
  do k=0,nn
    phi=0.0d0
    do kk=0,2*nn
       !***Husimi************!
     sumre = 0.0
     sumim = 0.0
     alabs=1.0d0*tan(theta/2.0)
       do ii=1,jj+1
         bkre  = 0.0d0
         bkim  = 0.0d0
         bkire = 0.0d0
         bkiim = 0.0d0
!      test=0.0
         do m=0,jj
           call binomial(2*jj,jj+(-jj+2*m),bin)
           bkre=bkre+m1p(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
           &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m)))*cos((jj+(-jj+2*m))*phi)
           bkim=bkim-m1p(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
           &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m)))*sin((jj+(-jj+2*m))*phi)
           bkire=bkire+m1p(m+1,ii)*(1.0/(1.+alabsi**2.0)**jj)*  &
           &  bin**(1./2.)*((alabsi)**(jj+(-jj+2*m)))*cos((jj+(-jj+2*m))*phii)
           bkiim=bkiim+m1p(m+1,ii)*(1.0/(1.+alabsi**2.0)**jj)*  &
           &  bin**(1./2.)*((alabsi)**(jj+(-jj+2*m)))*sin((jj+(-jj+2*m))*phii)
           enddo
           sumre = sumre   &
                & + bkire*bkre*cos(eigp(ii)*t) - bkire*bkim*sin(eigp(ii)*t) &
                & + bkiim*bkim*cos(eigp(ii)*t) + bkiim*bkre*sin(eigp(ii)*t)
           sumim = sumim   &
                & - bkire*bkim*cos(eigp(ii)*t) - bkire*bkre*sin(eigp(ii)*t) &
                & + bkiim*bkre*cos(eigp(ii)*t) - bkiim*bkre*sin(eigp(ii)*t)
       enddo

       do ii=1,jj
         bkre  = 0.0d0
         bkim  = 0.0d0
         bkire = 0.0d0
         bkiim = 0.0d0
         do m=0,jj-1
           call binomial(2*jj,jj+(-jj+1+2*m),bin)
           bkre=bkre+m1n(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
           &  bin**(1./2.)*((alabs)**(jj+(-jj+1+2*m)))*cos((jj+(-jj+1+2*m))*phi)
           bkim=bkim-m1n(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
           &  bin**(1./2.)*((alabs)**(jj+(-jj+1+2*m)))*sin((jj+(-jj+1+2*m))*phi)
           bkire=bkire+m1n(m+1,ii)*(1.0/(1.+alabsi**2.0)**jj)*  &
           &  bin**(1./2.)*((alabsi)**(jj+(-jj+1+2*m)))*cos((jj+(-jj+1+2*m))*phii)
           bkiim=bkiim+m1n(m+1,ii)*(1.0/(1.+alabsi**2.0)**jj)*  &
           &  bin**(1./2.)*((alabsi)**(jj+(-jj+1+2*m)))*sin((jj+(-jj+1+2*m))*phii)
         enddo
           sumre = sumre   &
                & + bkire*bkre*cos(eign(ii)*t) - bkire*bkim*sin(eign(ii)*t) &
                & + bkiim*bkim*cos(eign(ii)*t) + bkiim*bkre*sin(eign(ii)*t)
           sumim = sumim   &
                & - bkire*bkim*cos(eign(ii)*t) - bkire*bkre*sin(eign(ii)*t) &
                & + bkiim*bkre*cos(eign(ii)*t) - bkiim*bkre*sin(eign(ii)*t)
       enddo

       u=(1.0d0-cos(theta))*cos(phi)
       vv=(1.0d0-cos(theta))*sin(phi)

       jz=-cos(theta)
       QQ=sqrt(2*(1+jz))*cos(phi)
       PP=-sqrt(2*(1+jz))*sin(phi)
       
       Qh= 1.0d0*(sumre**2.0+sumim**2.0)
 
       write(1,*)QQ,'  ',PP,' ', Qh 
       
       phi=phi+delta
    enddo
    theta=theta+delta
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


!********   Subroutine which calculates binomial factors ****!
 
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
