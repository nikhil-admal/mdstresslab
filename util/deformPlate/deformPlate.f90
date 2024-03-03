program deformPlate
!
!  This is the main driver for the molecular dynamics program.
!
use Plate
use Continuum_Stress
implicit none

double precision,parameter :: PI = 3.141592653589793D0
integer, parameter :: DIM=3         ! change DIM to switch between 2D and 3D !
integer :: N=0
double precision, dimension(DIM) :: pos, new_pos,Mass_Center
double precision, dimension(:,:), allocatable    :: positions
character*80     :: SampIn        !  Name of file containing input sample
character*80     :: SampOut       !  Name of file to contain output sample

integer :: i,lis

! Read input and output filenames
read(*,'(a)',end=200,err=300) SampIn                       ! "
read(*,'(a)',end=200,err=300) SampOut                      ! "
! Read E, nu, G, rhole, psi
read(*,*,end=200,err=300) E
read(*,*,end=200,err=300) nu
read(*,*,end=200,err=300) G
read(*,*,end=200,err=300) force
read(*,*,end=200,err=300) rhole
read(*,*,end=200,err=300) psi
psi = psi * PI/180.d0                   ! Convert to radians

call Prep_constants

! open the input file
lis = len_trim(SampIn)
open(unit=1,file=SampIn(1:lis),status='old', &
     action='read',err=700)   
! open the output file
lis = len_trim(SampOut)
open(unit=2,file=SampOut(1:lis),status='unknown', &
     action='write',err=600)   

! Read number of atoms
read(1,*,end=800,err=900) N

!  read the coordinates (one line per atom) and write the deformed coordinate
!
allocate(positions(DIM,N))
do i=1,N
   read(1,*,end=800,err=900) pos
   positions(:,i) = pos
enddo

!  compute center of mass coordinates
!
Mass_Center = sum( positions , dim=2 ) / N

do i=1,N
   pos = positions(:,i) - Mass_Center
   call Anal_Sol(pos(1),pos(2))
   new_pos(1) = positions(1,i)+u
   new_pos(2) = positions(2,i)+v
   new_pos(3) = positions(3,i)
   write(2,'(1X,3E23.15)',err=1000) new_pos
enddo

close(unit=1)
close(unit=2)
return
!
200 continue
   print*,'Read_Input: FATAL: premature end-of-file in standard input'
   stop
300 continue
   print*,'Read_Input: FATAL: read error in standard input'
   stop
600 continue
   print*,'Read_Sample: FATAL: ',SampOut(1:lis),' not found.'
   stop
700 continue
   print*,'Read_Sample: FATAL: ',SampIn(1:lis),' not found.'
   stop
800 continue
   print*,'Read_Sample: FATAL: premature end-of-file at atom ',i
   close(unit=1)
   stop   
900 continue
   print*,'Read_Sample: FATAL: read error in ',SampIn(1:lis)
   close(unit=1)
   stop
1000 continue
   print*,'Write_Sample: WARNING: error when writing ',SampOut(1:lis)
   close(unit=2)
   return

end program deformPlate



