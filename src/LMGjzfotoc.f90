!******************
! compile as:
! gfortran -o LGMjz.exe ./LMGV3expectationjz.f90 -llapack -lblas
!************************
!------------------------------------------------------------------------!
! Author Daniel Julian Nader                                             !
! Collaboration Sergio Lerma                                             !
!                                                                        !
!  This code:                                                            !        
!     -  diagonalize LGM Model in SU(2) representation                   !
!     -  calculates the expectation value <z0|Jz^2(t)|z0>                !
!     -  calculates fotoc  <z0|Jz^2(t)|z0> - <z0|Jz(t)|z0>^2             !
!                                                                        !
! The output is printed in the file: jz_expectation.dat                  !
!                                                                        !
! Gnuplot visualization:                                                 !
!                     gnulplot>plot 'HusimiPQp+.dat' u 1:2 w l           !
!************************************************************************!


!---------------!
 program LGMmodel
!---------------!

   !
   !********* variables definition******************************
   !
   
 implicit none

 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(16)

 integer :: i,j,jq,jj,k,kp,kkk,n,nMC,m,mm
 real(8), allocatable :: m0p(:,:),m1p(:,:),m2p(:,:),eigp(:)
 real(8), allocatable :: m0n(:,:),m1n(:,:),m2n(:,:),eign(:)
 real(8), allocatable :: Q(:),ThetaL(:),PhiL(:)
 real(8), allocatable :: bpre(:),bnre(:),bpim(:),bnim(:)
 real(8) :: epsilon,uMC, Ps, w, v, gammax, gammay  ! Hamiltonian
 real(DP) :: bkre,bkim,bkpre,bkpim,sumint,sumint2
 real(DP) :: absbksqrt,theta,phi,alabs  
 real(DP) :: sumre,sumim,dt
 real(DP) :: sumre2,sumim2
 real(DP) :: Q2temp,QM2,Q2sum,Qpsitemp,Qpsi,Wtemp,Wpsi  ! Husimi

 real(DP) :: Qh,u,vv,delta,bin,pi,test

 character(len = 50) :: outputfile 

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
 print*,'theta (initial coherent state)'
 read*,theta
 print*,'phi (initial coherent state)'
 read*,phi
 print*,'output file name'
 read*,outputfile
  
 
 v=epsilon*(gammax-gammay)/(2.0*(2.0*jj-1.0))
 w=epsilon*(gammax+gammay)/(2.0*(2.0*jj-1.0))

  ! ************   Diagonalization       ****************

 
! Defining the size of the matrix


 
 allocate (m0p(jj+1,jj+1))
 allocate (m1p(jj+1,jj+1))
 allocate (m2p(jj+1,jj+1))
 allocate (eigp(jj+1))
 allocate (bpre(jj+1),bpim(jj+1))
 allocate (m0n(jj,jj))
 allocate (m1n(jj,jj))
 allocate (m2n(jj,jj))
 allocate (eign(jj))
 allocate (bnre(jj),bnim(jj))
 allocate (Q(30000))



 
!!$ !********************************************************************
!!$ ! ************ Calculating the matriz elements**********************
!!$ !********************************************************************
 ! starting all elements from zero
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


 !!**** elements on the diagonal**********!

 do i=0,jj
    m0p(i+1,i+1)=epsilon*(-jj+2*i)+(w)*(jj*(jj+1)-(-jj+2*i)**2)    ! <------   !!!paridad positiva
 enddo

 do i=0,jj-1
    m0n(i+1,i+1)=epsilon*(-jj+2*i+1)+(w)*(jj*(jj+1)-(-jj+2*i+1)**2)    ! <------   !!!paridad negativa
 enddo



 !!**** Elements out of the diagonal**********!


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
 !                          Printing results
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
   print*,'Husimi Representation  *'
   print*,'see file:              *'
   print*,'   fotoc_jz.dat        *'
   print*,'                       *'
   print*,'************************'



!******************************************************************************
!****  Calculating the expectation value <Jz>                *********!
!!$!******************************************************************************
  open(1,file=outputfile,status='new')

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
  write(1,*)'#     t         <J_z>         fotoc                    '       


   alabs=1.0d0*tan(theta/2.0)
   !****** calculating bk=<E_k|alpha0> positive parity*********************************************!
   do k=1,jj+1
   do m=0,jj
   call binomial(2*jj,jj+(-jj+2*m),bin)
   bpre(k)=bpre(k)+m1p(m+1,k)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m)))*cos((jj+(-jj+2*m))*phi)
   bpim(k)=bpim(k)+m1p(m+1,k)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+2*m)))*sin((jj+(-jj+2*m))*phi)
   enddo
   enddo
   
   !****** calculating bk=<E_k|alpha0> negative parity*********************************************!
   do k=1,jj
   do m=0,jj-1
   call binomial(2*jj,jj+(-jj+1+2*m),bin)
   bnre(k)=bnre(k)+m1n(m+1,k)*(1.0/(1.+alabs**2.0)**jj)*  &
         &  bin**(1./2.)*((alabs)**(jj+(-jj+1+2*m)))*cos((jj+(-jj+1+2*m))*phi)
   bnim(k)=bnim(k)+m1n(m+1,k)*(1.0/(1.+alabs**2.0)**jj)*  &
        &  bin**(1./2.)*((alabs)**(jj+(-jj+1+2*m)))*sin((jj+(-jj+1+2*m))*phi)
   enddo
   enddo




  !****   loop for time ************************************!
  do kkk=0,500                       !  
     sumre = 0.0
     sumim = 0.0
     sumre2 = 0.0
     sumim2 = 0.0
   !****   loop for k ************************************!
   do k=1,jj+1                        !  
   !****   loop for kprime ************************************!
   do kp=1, jj+1
   !****   loop for M ************************************!
   sumint=0.0d0
   sumint2=0.0d0
   !****  calculating sum_M <JM|E_k><JM|E_kp> (positive parity) **********************!
   do mm=0,jj
      sumint=sumint+m1p(mm+1,k)*m1p(mm+1,kp)*(-jj+2*mm)               ! This sum is introduced for J_z|JM>
      sumint2=sumint2+m1p(mm+1,k)*m1p(mm+1,kp)*(-jj+2*mm)*(-jj+2*mm)   ! This sum is introduced for J_z^2|JM>
   enddo
   sumre=sumre+((bpre(k)*bpre(kp)+bpim(k)*bpim(kp))*cos(kkk*dt*(eigp(k)-eigp(kp))))*sumint+ &
        & ((-bpim(k)*bpre(kp)+bpre(k)*bpim(kp))*sin(kkk*dt*(eigp(k)-eigp(kp))))*sumint
   sumre2=sumre2+((bpre(k)*bpre(kp)+bpim(k)*bpim(kp))*cos(kkk*dt*(eigp(k)-eigp(kp))))*sumint2+ &
         & ((-bpim(k)*bpre(kp)+bpre(k)*bpim(kp))*sin(kkk*dt*(eigp(k)-eigp(kp))))*sumint2
   sumim=sumim-((bpre(k)*bpre(kp)+bpim(k)*bpim(kp))*sin(kkk*dt*(eigp(k)-eigp(kp))))*sumint+ &
        & ((-bpim(k)*bpre(kp)+bpre(k)*bpim(kp))*cos(kkk*dt*(eigp(k)-eigp(kp))))*sumint
   sumim2=sumim2-((bpre(k)*bpre(kp)+bpim(k)*bpim(kp))*sin(kkk*dt*(eigp(k)-eigp(kp))))*sumint2+ &
         & ((-bpim(k)*bpre(kp)+bpre(k)*bpim(kp))*cos(kkk*dt*(eigp(k)-eigp(kp))))*sumint2
   enddo
   enddo
  !*****  loop for k negative parity
   do k=1,jj                          
   !****  loop for bkprime negative parity ************************************!
   do kp=1,jj
   !****   loop for M ************************************!
     sumint=0.0d0
     sumint2=0.0d0
   !****  calculating sum_M <JM|E_k><JM|E_kp> (negative parity) **********************!
   do mm=0,jj-1
      sumint=sumint+m1n(mm+1,k)*m1n(mm+1,kp)*(-jj+2*mm+1)                    ! This sum is introduced for J_z|JM>
      sumint2=sumint2+m1n(mm+1,k)*m1n(mm+1,kp)*(-jj+2*mm+1)*(-jj+2*mm+1)     ! This sum is introduced for J_z^2|JM>
   enddo
   sumre=sumre+((bnre(k)*bnre(kp)+bnim(k)*bnim(kp))*cos(kkk*dt*(eign(k)-eign(kp))))*sumint+ &
        & ((-bnim(k)*bnre(kp)+bnre(k)*bnim(kp))*sin(kkk*dt*(eign(k)-eign(kp))))*sumint
   sumre2=sumre2+((bnre(k)*bnre(kp)+bnim(k)*bnim(kp))*cos(kkk*dt*(eign(k)-eign(kp))))*sumint2+ &
        & ((-bnim(k)*bnre(kp)+bnre(k)*bnim(kp))*sin(kkk*dt*(eign(k)-eign(kp))))*sumint2
   sumim=sumim-((bnre(k)*bnre(kp)+bnim(k)*bnim(kp))*sin(kkk*dt*(eigp(k)-eigp(kp))))*sumint+ &
        & ((-bnim(k)*bnre(kp)+bnre(k)*bnim(kp))*cos(kkk*dt*(eigp(k)-eigp(kp))))*sumint
   sumim2=sumim2-((bnre(k)*bnre(k)+bnim(k)*bnim(kp))*sin(kkk*dt*(eigp(k)-eigp(kp))))*sumint2+ &
         & ((-bnim(k)*bnre(kp)+bnre(k)*bnim(kp))*cos(kkk*dt*(eigp(k)-eigp(kp))))*sumint2
   enddo
   enddo


 
   write(1,*)kkk*dt,' ', sumre/(1.0*jj),' ',(sumre2-sumre**2)/(1.0*jj)**2

   
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


!******** subroutine which calculates binomial factors ****!
 
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
