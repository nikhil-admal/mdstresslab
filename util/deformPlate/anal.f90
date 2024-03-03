! Nikhil
! rhole = radius of the hole in the plate
! youngs_mod = Young's Modulus 
! nu = Poisson's ratio
! rig_mod = Rigidity Modulus
! psi = Angle at which the force is being applied
! force = Magnitude of the force being applied

! Added the following Modules, Subroutines and Functions related to
! analytical solution
! Modules:     Plate
!              Plate_constants
!              stress
! Subroutines: Anal_Sol
!              Prep_constants
!              roots
! Functions:   x_disp
!              y_disp

module Plate
    double precision             ::      rhole
    double precision             ::      E,nu,G
    double precision             ::      force
    double precision             ::      psi
end module Plate

module Plate_constants
    double precision,dimension(6,6) :: beta
    complex(kind=8),dimension(2)    :: mu
    complex(kind=8)                 :: a1,b1
end module Plate_constants    

module Continuum_Stress
    double precision :: sigx,sigy,Txy,u,v
    double precision::  sigx0,sigy0,Txy0,u0,v0
end module Continuum_Stress

subroutine Prep_constants
use Plate
use Continuum_Stress
use Plate_constants
implicit none
complex(kind=8)                 :: x0=(1.d0,1.d0)
double precision,dimension(6,6) :: alpha
integer                         :: i,j

alpha=reshape( (/ 1.d0/E ,  -nu/E  ,  -nu/E  ,    0.d0   ,    0.d0  ,    0.d0  , &
                & -nu/E  ,  1.d0/E ,  -nu/E  ,    0.d0   ,    0.d0  ,    0.d0  , &
                & -nu/E  ,  -nu/E  ,  1.d0/E ,    0.d0   ,    0.d0  ,    0.d0  , &
                &   0.d0 ,    0.d0 ,    0.d0 ,    1.d0/G  ,   0.d0  ,    0.d0  , &
                &   0.d0 ,    0.d0 ,    0.d0 ,    0.d0   ,    1.d0/G  ,  0.d0  , &
                &   0.d0 ,    0.d0 ,    0.d0 ,    0.d0   ,    0.d0  ,    1.d0/G  /),(/6,6/) )

beta = reshape( (/((alpha(i,j)-alpha(i,3)*alpha(j,3)/alpha(3,3),i=1,6),j=1,6)/),(/6,6/) )
a1 = (0.d0,0.d0)
b1 = cmplx(0,-force*rhole/2)
call roots(beta,mu,x0)
sigx0 = force*COS(psi)**2
sigy0 = force*SIN(psi)**2
Txy0  = force*SIN(psi)*COS(psi)
end subroutine Prep_constants

! x0 = intial guess

subroutine roots(beta,mu,x0)
implicit none
double precision, dimension(6,6), intent(in) :: beta
complex(kind=8), dimension(2), intent(out)   :: mu
complex(kind=8),intent(inout)                :: x0
! Local variables
complex(kind=8)  :: pol,dpol,x1
double precision :: eps,delta,r,i
eps = 1.d-16
delta = 1.d-16

pol = beta(1,1)*x0**4 - 2*beta(1,6)*x0**3 + (2*beta(1,2)+beta(6,6))*x0**2 - 2*beta(2,6)*x0 + beta(2,2) 
dpol = 4*beta(1,1)*x0**3 - 6*beta(1,6)*x0**2 + 2*(2*beta(1,2)+beta(6,6))*x0 - 2*beta(2,6)
if(abs(pol) .ne. 0.d0 .and. abs(dpol) .ne. 0.d0) then
    do
        x1 = x0          
        x0 = x0 - pol/dpol
        pol = beta(1,1)*x0**4 - 2*beta(1,6)*x0**3 + (2*beta(1,2)+beta(6,6))*x0**2 - 2*beta(2,6)*x0 + beta(2,2) 
        if(abs(x0-x1) < eps .or. abs(pol)<delta) exit
    enddo
endif
    
mu(1) = x0
r = beta(1,6)/beta(1,1) - real(x0)
i = sqrt(beta(2,2)/(beta(1,1)*abs(x0)**2) - r**2)
mu(2) = cmplx(r,i)  
!write(*,*) mu
end subroutine roots


subroutine Anal_Sol(x,y)
use Plate
use Plate_constants
use Continuum_Stress
implicit none
double precision,external    :: x_disp,y_disp
double precision,intent(in)  :: x,y
!local variables
complex(kind=8),dimension(2) :: p,q,xn,yn,z
complex(kind=8)              :: mu_diff
complex(kind=8),dimension(2) :: zeta,dzeta,zetainv,delz,phi,dphi
complex(kind=8),parameter    :: i = (0,1)
integer                      :: k

!write(*,*) mu(1), mu(2)
u = 0.d0
v = 0.d0
sigx = 0.d0
sigy = 0.d0
Txy = 0.d0
if (rhole == 0.d0) then
    u = x_disp(x,y)
    v = y_disp(x,y)
    sigx = sigx0
    sigy = sigy0
    Txy = Txy0
else
        z = x+mu*y
        xn = sqrt(rhole**2 + (mu*rhole)**2)
    yn = i*xn
    delz = sqrt(z**2 - rhole**2 - (mu*rhole)**2)

    do k=1,2
        if( real(z(k)*conjg(xn(k))) * real(delz(k)*conjg(xn(k))) .gt. 0.d0 .and. real(z(k)*conjg(yn(k))) * real(delz(k)*conjg(yn(k))) .gt. 0.d0 ) then
            zeta(k) = (z(k) + delz(k))/(rhole-i*mu(k)*rhole)
            dzeta(k) = (1 + z(k)/delz(k)) / (rhole-i*mu(k)*rhole)
        else
            zeta(k) = (z(k) - delz(k))/(rhole-i*mu(k)*rhole)
            dzeta(k) = (1 - z(k)/delz(k)) / (rhole-i*mu(k)*rhole)
        endif
    enddo


    zetainv = 1/zeta
    mu_diff = mu(1)-mu(2)

    phi(1) =  zetainv(1) * (b1-mu(2)*a1)/mu_diff
    phi(2) = -zetainv(2) * (b1-mu(1)*a1)/mu_diff

    dphi(1) = -(zetainv(1)**2) * dzeta(1) * (b1-mu(2)*a1)/mu_diff
    dphi(2) =  (zetainv(2)**2) * dzeta(2) * (b1-mu(1)*a1)/mu_diff

    p = beta(1,1)*mu**2 + beta(1,2) - beta(1,6)*mu
    q = beta(1,2)*mu + beta(2,2)/mu - beta(2,6)

    sigx = 2*real((mu(1)**2) * dphi(1) +(mu(2)**2) * dphi(2)) + sigx0
    sigy = 2*real(dphi(1) + dphi(2)) + sigy0
    Txy = -2*real(mu(1)*dphi(1) + mu(2)*dphi(2)) + Txy0
    u = 2*real(p(1)*phi(1) + p(2)*phi(2)) + x_disp(x,y)
    v = 2*real(q(1)*phi(1) + q(2)*phi(2)) + y_disp(x,y)
endif
!write(92,'(9E20.10)') x,y,u,v,sigx,sigy,Txy,nu,E

end subroutine Anal_Sol

! The following functions are only valid for psi=0
double precision function x_disp(x,y)
use Plate
implicit none
double precision,intent(in) :: x,y
!local variables
double precision            :: strain
strain = (1-nu**2)*force/E
x_disp = x*strain
end function x_disp

double precision function y_disp(x,y)
use Plate
implicit none
double precision,intent(in) :: x,y
!local variables
double precision            :: strain
strain = -nu*(1+nu)*force/E
y_disp = y*strain
end function y_disp
