!******************
! compile as:
! gfortran -o SurProb.exe ./src/LMGV3SurProbBlochCS.f90 -llapack -lblas
!************************
!------------------------------------------------------------------------!
! Author Daniel Julian Nader                                             !
! Collaboration Sergio Lerma                                             !
!                                                                        !
!  This code:                                                            !        
!     -  diagonalize LGM Model in SU(2) representation                   !
!     -  Calculate  the survival probability SP of an initial Bloch      !
!        coherent state |z0(theta0,phi0)>                                !
!                                                                        !
! The output for each parity are printed in the files  - Sp_cs.dat       !
!                                                                        !
! Gnuplot visualization:                                                 !
!                     gnulplot>plot 'SP_cs.dat' u 1:2 w l                !
!************************************************************************!

!---------------!
 program LGMmodel
!---------------!

   !
   !********* variables definition******************************
   !
   
 implicit none

 integer k,nn

 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(16)

 integer :: i,ii,j,jq,jj,kkk,n,nMC,m,mq
 real(8), allocatable :: m0p(:,:),m1p(:,:),m2p(:,:),eigp(:)
 real(8), allocatable :: m0n(:,:),m1n(:,:),m2n(:,:),eign(:)
 real(8), allocatable :: Q(:),ThetaL(:),PhiL(:)
 real(8) :: epsilon,uMC, Ps, w, v, gammax, gammay  ! Hamiltonian
 real(DP) :: bkre,bkim,absbksqrt,theta,phi,alabs  
 real(DP) :: sumre,sumim,dt
 real(DP) :: Q2temp,QM2,Q2sum,Qpsitemp,Qpsi,Wtemp,Wpsi  ! Husimi

 real(DP) :: Qh,u,vv,delta,bin,pi,test

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
 print*,'theta (coherent state)'
 read*,theta
 print*,'phi (coherent state)'
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
 allocate (Q(30000))




 
!!$ !********************************************************************
!!$ ! ************  Calculating matrix elements**********************
!!$ !********************************************************************
 ! all matrix elements starts in zero
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
    m0p(i+1,i+1)=epsilon*(-jj+2*i)+(w)*(jj*(jj+1)-(-jj+2*i)**2)    ! <------   !!!paridad positiva
 enddo

 do i=0,jj-1
    m0n(i+1,i+1)=epsilon*(-jj+2*i+1)+(w)*(jj*(jj+1)-(-jj+2*i+1)**2)    ! <------   !!!paridad negativa
 enddo



 !!**** Elements out of the diagonal**********!

    do i=0,jj-1
       m0p(i+1,i+2)=1.0*(v/2.0)*(jj*(jj+1)-(-jj+2*i+2)*((-jj+2*i+2)-1))**(1.0/2.0)* &
           & (jj*(jj+1)-((-jj+2*i+2)-1)*((-jj+2*i+2)-2))**(1.0/2.0)  ! <------   !!!paridad positiva            
    enddo

    do i=0,jj-2
       m0n(i+1,i+2)=1.0*(v/2.0)*(jj*(jj+1)-(-jj+2*i+3)*((-jj+2*i+3)-1))**(1.0/2.0)* &
           & (jj*(jj+1)-((-jj+2*i+3)-1)*((-jj+2*i+3)-2))**(1.0/2.0) ! <------   !!!paridad negativa
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

   print*,'************************'
   print*,'Survival Probability   *'
   print*,'see file:              *'
   print*,'         SP_cs.dat     *'
   print*,'                       *'
   print*,'************************'



!******************************************************************************
!*****  Printing the survival probability to a file                 *********!
!!$!******************************************************************************
  open(1,file='SP_cs.dat',status='new')

  dt=0.01
  write(1,*)'# ********Parameters**********  #'
  write(1,*)'# angular momentum',jj
  write(1,*)'# Epsilon',epsilon
  write(1,*)'# gammax',gammax
  write(1,*)'# gammay',gammay
  write(1,*)'# Epsilon',epsilon
  write(1,*)'# Theta (coherent state)',theta
  write(1,*)'# Phi   (coherent state)',phi
  write(1,*)'# positive parity '
  write(1,*)'# dt #',dt
  write(1,*)'# ****************************  #'
  write(1,*)'#     t             P   '       

  
  
  do kkk=0,10000
     alabs=1.0d0*tan(theta/2.0)
     sumre = 0.0
     sumim = 0.0
   do ii=1,jj+1
      bkre = 0.0d0
      bkim = 0.0d0
!      test=0.0
   do m=0,jj
   call binomial(2*jj,jj+(-jj+2*m),bin)
   bkre=bkre+m1p(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m)))*cos((jj+(-jj+2*m))*phi)
   bkim=bkim+m1p(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m)))*sin((jj+(-jj+2*m))*phi)
   enddo
   absbksqrt= bkre**2.0 + bkim**2.0
   sumre = sumre + absbksqrt*cos(eigp(ii)*kkk*dt)
   sumim = sumim - absbksqrt*sin(eigp(ii)*kkk*dt)
   enddo

   do ii=1,jj
   bkre = 0.0d0
   bkim = 0.0d0
   do m=0,jj-1
   call binomial(2*jj,jj+(-jj+1+2*m),bin)
   bkre=bkre+m1n(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+1+2*m)))*cos((jj+(-jj+1+2*m))*phi)
   bkim=bkim+m1n(m+1,ii)*(1.0/(1.+alabs**2.0)**jj)*  &
        &  bin**(1./2.)*((alabs)**(jj+(-jj+1+2*m)))*sin((jj+(-jj+1+2*m))*phi)
   enddo
   absbksqrt= bkre**2.0 + bkim**2.0
   sumre = sumre + absbksqrt*cos(eign(ii)*kkk*dt)
   sumim = sumim - absbksqrt*sin(eign(ii)*kkk*dt)
   enddo

   Ps= 1.0d0*(sumre**2.0+sumim**2.0)
 
   write(1,*)kkk*dt,'  ',Ps 

   
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


!********   Subroutine for calculatin binomial factors ****!
 
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
