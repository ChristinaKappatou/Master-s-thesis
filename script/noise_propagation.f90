program noise_propagation
implicit none
complex, allocatable :: beta(:)
complex, allocatable :: m1p(:), m1e(:), m1w(:), m2p(:), m2e(:), m2w(:)
complex, allocatable :: psi(:), psinew(:)
complex, allocatable :: source(:)
complex, parameter   :: im=(0,1)  
integer              :: i, j, nr, nz, rcounter 
integer              :: ic, approx, absorption
integer              :: k
double precision     :: c0, fr, brefr, z0, flowres 
double precision     :: lamda, dr, dz, pi, omega, pref, rmax, zmax, zsource, zlayer
double precision     :: r, z, csound, at
double precision     :: p2av, spl, dl
double precision     :: kfr
double precision     :: t01, temp, csat, rhosat, rhor, rh, hum, fro, t20, tr, frn, b1, b2, abscoef, ki
complex              :: pfree
complex              :: alpha, gama, sigma1, sigma2, tau1, tau2, kwav, dk2, k0
complex              :: p, q
complex              :: zeta, reflcoef

external tridag

open(1,file='input.txt')
open(2,file='output.dat')
open(3,file='dr.dat')

read(1,*)c0           !ground speed of sound [m/s]
read(1,*)brefr        !refraction factor [m/s]
read(1,*)z0           !roughness lenght [m]
read(1,*)nr           !r-wise section number
read(1,*)nz           !z-wise section number
read(1,*)flowres      !flow resistivity [Pa*s/m2]
read(1,*)zsource      !source height [m]
read(1,*)ic           !initial conditions
read(1,*)approx       !approximation type
read(1,*)absorption   !with/without atm.abs.
read(1,*)t01          !water triple point temp.
read(1,*)t20          !reference temp.
read(1,*)temp         !absolute temp.
read(1,*)rhor         !relative pressure
read(1,*)rh           !realtive humidity

allocate (beta(nz))
allocate (m1p(nz), m1e(nz), m1w(nz), m2p(nz), m2e(nz), m2w(nz), psi(nz), psinew(nz), source(nz))

!*every equation is with reference to Salomons, "Computational Atmospheric Acoustics"

print*,"start program"

!-----loop for spectre calculations-----------------------------------
do k=1,17
 kfr=k-10.d0-1.d0
 fr=1000.d0*10.d0**(0.1*kfr)      !1/3 octave band frequencies 
! fr=1000.d0
!-----------------------------------------------------------------------term definitions body-----------------------------------------------------------------------------------------------
 lamda=c0/fr                      !wave length [m]
 dr=lamda/10.d0                   !horizontal grid spacing [m]
! dr=1.d0/30.d0
 do rcounter=1,120
!  print*,1.d0/rcounter
  if (fr.lt.400.d0.and.abs(dr-1.d0/rcounter).lt.0.025d0.or.fr.ge.400.d0.and.abs(dr-1.d0/rcounter).lt.0.003d0) then  !grid density to fit each frequency
!   print*,fr,rcounter
   dr=1.d0/rcounter 
   goto 15
  else
   dr=lamda/10.d0
  endif
 enddo
!15 write (3,*) fr, dr
!15 dz=lamda/10.d0                   !vertical grid spacing [m]
15 dz=dr
 zmax=dz*nz                       !upper bound height [m]
 rmax=dr*nr                       !right bound lenght [m] 
 pi=4.d0*atan(1.d0)
 omega=2.d0*pi*fr                 !angular frequency [rad/s]
 pref=2.e-5                       !reference pressure [Pa]
 zlayer=zmax-50.d0*lamda          !start of absorbing layer [m]
! zeta=1.d0+9.08d0*(1000.d0*fr/flowres)**(-0.75d0)+im*(11.9d0*(1000.d0*fr/flowres)**(-0.73d0)) !normalized impedance
 zeta=5.96+im*2.46
 reflcoef=(zeta-1.d0)/(zeta+1.d0) !reflection coefficient (G.77)

!-------atmospheric absorption coefficient (sect.B5)--------------------------------------------------------------------------------------------
 csat=4.6151d0-6.8346d0*(t01/temp)**1.261d0
 rhosat=10.d0**csat
 hum=rh*rhosat/rhor
 fro=rhor*(24.d0+40400.d0*hum*(0.02d0+hum)/(0.391d0+hum))
 tr=temp/t20
 frn=(rhor/sqrt(tr))*(9.d0+280.d0*hum*exp(-4.17d0*(tr**(-1.d0/3.d0)-1.d0)))
 b1=0.1068d0*(exp(-3352.d0/temp))/(frn+(fr**2.d0)/frn)
 b2=0.01275d0*(exp(-2239.1d0/temp))/(fro+(fr**2.d0)/fro)

 if (absorption.eq.0) then
  abscoef=0.d0                                                                          !ignores atm. absorption
 elseif (absorption.eq.1) then
  abscoef=8.686d0*(fr**2.d0)*(sqrt(tr))*(1.84d0/(rhor*10.d0**11.d0)+(b1+b2)/(tr**3.d0)) !calculates atm. absorption
 endif

 ki=abscoef/(20.d0*log10(exp(1.d0)))

!-------at calculation (top absorbing layer coefficient)------------------------------------------------------------------------
!(linear interpolation has been used for calculations)
 if(fr.ge.30.d0.and.fr.lt.125.d0) then
  at=0.2d0+0.2d0*(fr-30.d0)/95.d0
 elseif(fr.ge.125.d0.and.fr.lt.500.d0) then
  at=0.4d0+0.1d0*(fr-125.d0)/375.d0
 elseif(fr.ge.500.d0.and.fr.lt.1000.d0) then
  at=0.5d0+0.5d0*(fr-500.d0)/500.d0
 elseif(fr.ge.1000.d0) then
  at=1.d0
 endif
 
!------------------------------------------------------------------------------------------------------------------------------------ 
 k0=omega/c0+im*ki                           !wave number at ground surface
 alpha=im/(2.d0*k0)                          !matrix coefficients
 gama=alpha/(dz**2.d0)
!-----------------------------------------------------------------------end of term definitions body-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------initial conditions body----------------------------------------------------------------------------------------------
 if (ic.eq.1) then
  goto 1
 elseif (ic.eq.2) then 
  goto 2
 endif 
 
!-------1st order i.c.-----------------------------------------------------------------------------------------------------
1 sigma1=1.d0/(1.d0-im*k0*dz/zeta)
 sigma2=0.d0
 tau1=1.d0/(1.d0+im*k0*dz)
 tau2=0.d0
 goto 3

!-------2nd order i.c.------------------------------------------------------------------------------------------------------
2 sigma1=4.d0/(3.d0-2.d0*im*k0*dz/zeta)
 sigma2= -1.d0/(3.d0-2.d0*im*k0*dz/zeta)
 tau1=4.d0/(3.d0+2.d0*im*k0*dz)
 tau2=-1.d0/(3.d0+2.d0*im*k0*dz)
 goto 3

!----------------------------------------------------------------------------------------------------------------------------
3 do j=1,nz
  z=j*dz
  csound=c0+brefr*log(1.d0+z/z0)                                          !logarithmic sound speed profile (ch.G.9)
   if (z.ge.zlayer.and.z.le.zmax) then
    kwav=omega/csound+im*ki+im*at*((z-zlayer)**2.d0)/((zmax-zlayer)**2.d0)!wave number calculation inside the absorbing layer
   else 
    kwav=omega/csound+im*ki                                               !wave number calculation underneath the absorbing layer
   endif
  dk2=kwav**2.d0-k0**2.d0

  beta(j)=(im*dk2)/(2.d0*k0)

!-------psi for r=0-initial condition-----------------------------------------------------------------------------------------
  if (approx.eq.1) then 
   goto 4
  elseif (approx.eq.0) then
   goto 5
  endif

!-------wide angle solution (G.75) in (G.76)                                         
4 psi(j)=(sqrt(im*k0))*((1.3717d0-0.3701d0*(k0**2.d0)*(z-zsource)**2.d0)*exp(-(1.d0/3.d0)*(k0**2.d0)*(z-zsource)**2.d0)+reflcoef*(1.3717d0-0.3701d0*(k0**2.d0)*(z+zsource)**2.d0)*exp(-(1.d0/3.d0)*(k0**2.d0)*(z+zsource)**2.d0))
  goto 6     

!-------narrow angle solution (G.64) in (G.76)
5 psi(j)=(sqrt(im*k0))*(exp(-(k0**2.d0)*((z-zsource)**2.d0)/2.d0)+reflcoef*exp(-(k0**2.d0)*((z+zsource)**2.d0)/2.d0)) 
  goto 6 
6 enddo
!------------------------------------------------------------------------end of initial conditions body--------------------------------------------------------------------------------------


!------------------------------------------------------------------------calculations body--------------------------------------------------------------------------------------------------
!-------vertical loop---------------------------------------------------------------------------------------------------------------------------------
  do i=1,nr
   r=i*dr
   if (approx.eq.0) then
    goto 7
   elseif (approx.eq.1) then
    goto 8
   endif

!--------narrow angle case matrices-------------------------------------------------------------------------------------------------------------------
!--------lower bound i.c. (G.32)-narrow angle approx.-------------------------------------------------------------------------------------------------
7  m2p(1)=1.d0-0.5d0*dr*(gama*(-2.d0+sigma1)+beta(1))
   m2e(1)= -0.5d0*dr*gama*(1.d0+sigma2)

   m1p(1)=1.d0+0.5d0*dr*(gama*(-2.d0+sigma1)+beta(1))
   m1e(1)=  0.5d0*dr*gama*(1.d0+sigma2)

   source(1)=m1p(1)*psi(1)+m1e(1)*psi(2)

!-------upper bound i.c. (G.32)-narrow angle approx.--------------------------------------------------------------------------------------------------
   m2p(nz)=1.d0-0.5d0*dr*(gama*(-2.d0+tau1)+beta(nz))
   m2w(nz)= -0.5d0*dr*gama*(1.d0+tau2)

   m1p(nz)=1.d0+0.5d0*dr*(gama*(-2.d0+tau1)+beta(nz))
   m1w(nz)=  0.5d0*dr*gama*(1.d0+tau2)

   source(nz)=m1p(nz)*psi(nz)+m1w(nz)*psi(nz-1)

!-------body (horizontal loop) (G.32)-narrowangle approx.---------------------------------------------------------------------------------------------
   do j=2, nz-1
    m2w(j)= -0.5d0*dr*gama
    m2p(j)=1.d0-0.5d0*dr*(-2.d0*gama+beta(j))
    m2e(j)= -0.5d0*dr*gama

    m1w(j)=  0.5d0*dr*gama
    m1p(j)=1.d0+0.5d0*dr*(-2.d0*gama+beta(j))
    m1e(j)=  0.5d0*dr*gama

    source(j)=m1p(j)*psi(j)+m1w(j)*psi(j-1)+m1e(j)*psi(j+1)
   enddo
   goto 9
!---------------------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------------

!--------wide angle case matrices-------------------------------------------------------------------------------------------------------------------
!--------lower bound i.c. (G.33)-wide angle approx.-------------------------------------------------------------------------------------------------
8 m2p(1)=1.d0+(-0.5d0*dr+1.d0/(2.d0*im*k0))*(gama*(-2.d0+sigma1)+beta(1))
  m2e(1)= (-0.5d0*dr+1.d0/(2.d0*im*k0))*gama*(1.d0+sigma2)

  m1p(1)=1.d0+(0.5d0*dr+1.d0/(2.d0*im*k0))*(gama*(-2.d0+sigma1)+beta(1))
  m1e(1)=  (0.5d0*dr+1.d0/(2.d0*im*k0))*gama*(1.d0+sigma2)

  source(1)=m1p(1)*psi(1)+m1e(1)*psi(2)

!-------upper bound i.c. (G.33)-wide angle approx.-------------------------------------------------------------------------------------------------
  m2p(nz)=1.d0+(-0.5d0*dr+1.d0/(2.d0*im*k0))*(gama*(-2.d0+tau1)+beta(nz))
  m2w(nz)=(-0.5d0*dr+1.d0/(2.d0*im*k0))*gama*(1.d0+tau2)

  m1p(nz)=1.d0+(0.5d0*dr+1.d0/(2.d0*im*k0))*(gama*(-2.d0+tau1)+beta(nz))
  m1w(nz)=  (0.5d0*dr+1.d0/(2.d0*im*k0))*gama*(1.d0+tau2)

  source(nz)=m1p(nz)*psi(nz)+m1w(nz)*psi(nz-1)

!-------body (horizontal loop) (G.33)-wide angle approx.------------------------------------------------------------------------------------------
  do j=2, nz-1
   m2w(j)= (-0.5d0*dr+1.d0/(2.d0*im*k0))*gama
   m2p(j)=1.d0+(-0.5d0*dr+1.d0/(2.d0*im*k0))*(-2.d0*gama+beta(j))
   m2e(j)= (-0.5d0*dr+1.d0/(2.d0*im*k0))*gama

   m1w(j)=  (0.5d0*dr+1.d0/(2.d0*im*k0))*gama
   m1p(j)=1.d0+(0.5d0*dr+1.d0/(2.d0*im*k0))*(-2.d0*gama+beta(j))
   m1e(j)=  (0.5d0*dr+1.d0/(2.d0*im*k0))*gama

   source(j)=m1p(j)*psi(j)+m1w(j)*psi(j-1)+m1e(j)*psi(j+1) !m1*psi-matrix multiplication
  enddo 
  goto 9
!-----------------------------------------------------------------------end of calculations body---------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------solution/results body------------------------------------------------------------------------------------------------
9 call tridag(m2w, m2p, m2e, source, psinew, nz)           !thomas algorithm solution for the tridiagonal matrix system
  do j=1,nz
   psi(j)=psinew(j)                                        !vector values update
   z=j*dz   
!-------pressure field----------------------------------------------------------------------------------------------------------------------------
   q=psi(j)*exp(im*k0*r)            !(G.4)
   p=q/sqrt(r)                      !(G.2)
   if ((abs(p)).lt.pref) then       !covers the case of too small pressure values
    p=pref
   endif    
!-------SPL-DL calculation------------------------------------------------------------------------------------------------------------------------
   p2av=0.5d0*(abs(p))**2.d0
   spl=10.d0*log10(p2av/pref**2.d0)                                                  !sound pressure level(2.8)

   pfree=exp(im*k0*sqrt((z-zsource)**2.d0+r**2.d0))/sqrt((z-zsource)**2.d0+r**2.d0)  !free field pressure (3.1) 
!   pfree=exp(im*k0)  !free field pressure (3.1) 
   dl=10.d0*log10(((abs(p))**2.d0)/((abs(pfree))**2.d0))                             !relative sound pressure level(3.6)
!--------------------------------------------------------------------------------------------------------------------------------------------------------
   if (abs(z-2.d0).lt.0.05d0.and.r.eq.30.d0) then       !used when calculating the spectre on a specific grid point
!   if (z.gt.50.00124d0.and.z.lt.50.00126d0) then
!   if (mod(j,15).eq.0) then                       !plots one every 2 lines
!    write(2, '(2f10.2, 1e20.7)')  r, z, dl   
    write (2,'(2f10.2, 1e20.7)') fr, dl 
   endif
  enddo 
!  write(2,*)
 enddo
! write(2,*)
enddo
!-----------------------------------------------------------------------end of solution/results body-----------------------------------------------------------------------------------------


close(1)
close(2)
deallocate (beta)
deallocate (m1p, m1e, m1w, m2p, m2e, m2w, psi, psinew, source)
print*,"end program"
end program

subroutine tridag(a,b,c,r,k,n)
integer :: n,j
complex :: gam(n), a(n), b(n), c(n), r(n), k(n)
complex :: bet
if(b(1).eq.0.d0)pause 
bet=b(1)
k(1)=r(1)/bet
do 11 j=2,n
gam(j)= c(j-1)/bet
bet=b(j)-a(j)*gam(j)
if(bet.eq.0.d0)pause
k(j)=(r(j)-a(j)*k(j-1))/bet
11 continue
do 12 j=n-1,1,-1
k(j)=k(j)-gam(j+1)*k(j+1)
12 continue
return
end subroutine tridag 

