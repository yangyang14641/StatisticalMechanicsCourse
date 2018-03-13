!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Program name: ThreeBodyMD.f90
!
! Program Purpose:
!                  Solve three-bodies probelm in molecule simulation.
!  Version message:
!        Date          Programmer            Description of change
!      ========       ============          =========================
!  1.  13/3/2016      Yang Yang :-)               Source Code
! 
!
!
! A brief algorithm introduction: 
! First time step is caculated using velocity-verlet algorithm and other time step is
! caculated using verlet algorithm. Velocity is caculated using centered differential scheme.
! $r(0 + \delta t) = 2r(0) - \delta t*v(0) + \frac{1}{2} \delta \t^{2} a(0)$
!!$v(0 + \delta t) = v(t) + \frac{1}{2} \delta t [a(0) + a(0 + \delta t)]$
! $r(t+\delta t) = 2r(t) - r(t - \delta t)+ \delta \t^{2}a(t)$
! $v(t) = [r(t + \delta t) - r(t - \delta t)] / 2\delta t$ 
! $v(N) = [r(N) - r(N - \delta t)] / \delta t$ 
! external field : $V(r) = \frac{1}{2} \r^{2} 
! internal reaction : Leonard-Jones 6-12 potential  $V(r_{ij}) = 4[(\frac{1}{r_{ij}})^{12} - (\frac{1}{r_{ij}})^{6}]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!------------------------------------------------------------------------------------------------
! Setup computing parameters 
MODULE ParametersBlock
implicit none

integer, parameter:: N_particles = 3;                                      ! Number of particle
integer, parameter:: N_dimension = 2;                                      ! Problem dimension

! real*8, parameter:: Mass = 1.d0                                          ! Particles' mass (this term is not necessary!)

integer, parameter:: N_time = 1000000000;                                  ! Number of time step
real*8, parameter::  delta_t = 1.d-6;                                      ! Time step size  (really important)

integer, parameter:: Output_interval = 100000                              ! Output interval

character(len=20), parameter:: OutputFile = 'result.txt'                   ! Output file name

real*8, parameter:: x1_ini = 0.d0, y1_ini = 1.d0                           ! Particle No.1's position
real*8, parameter:: x2_ini = -sqrt(3.d0) / 2.d0, y2_ini = -1.d0 / 2.d0     ! Particle No.2's position
real*8, parameter:: x3_ini = sqrt(3.d0) / 2.d0, y3_ini = -1.d0 / 2.d0      ! Particle No.3's position

!real*8, parameter:: x1_ini = sqrt(2.d0)/2.d0, y1_ini = sqrt(2.d0)/2.d0      ! Particle No.1's position
!real*8, parameter:: x2_ini = 0.d0, y2_ini = 1.d0                            ! Particle No.2's position
!real*8, parameter:: x3_ini = -1.d0, y3_ini = 0.d0                           ! Particle No.3's position

END MODULE ParametersBlock
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------------
! The main program
PROGRAM Main
use ParametersBlock
implicit none

! Variables
real*8 x_new(1:N_particles,1:N_dimension), x_current(1:N_particles,1:N_dimension), x_old(1:N_particles,1:N_dimension)   ! Particles' direction
real*8 u_current(1:N_particles,1:N_dimension)                                                                           ! Particles' velocity
real*8 acceleration(1:N_particles,1:N_dimension)                                                                        ! Particles' acceleration
real*8 potential_energy(1:N_particles), kinetic_energy(1:N_particles)                                                   ! particle's Potential energy and Kinetic energy

! Index
integer i,j,ierror
integer:: t = 0

! Initial all varibales 
do i = 1, N_particles
  do j = 1, N_dimension
      x_new(i,j) = 0.d0
      x_current(i,j) = 0.d0
      x_old(i,j) = 0.d0
      u_current(i,j) = 0.d0
      acceleration(i,j) = 0.d0 
  end do
end do

do i = 1, N_particles
    potential_energy(i) = 0.d0
    kinetic_energy(i) = 0.d0 
end do

! Initial particles' state
! particles' initial position
x_current(1,1) =  x1_ini
x_current(1,2) =  y1_ini

x_current(2,1) = x2_ini
x_current(2,2) = y2_ini

x_current(3,1) = x3_ini
x_current(3,2) = y3_ini

call compute_acceleration(x_current,acceleration)          ! Compute acceleration at first time step

! time step evolution
do i = 1,N_particles
   do j = 1,N_dimension
      x_new(i,j) = x_current(i,j) + delta_t*u_current(i,j) + 0.5d0*(delta_t**2)*acceleration(i,j) 
   end do
end do

! compute potential energy
call compute_potnetialEnergy(x_current,potential_energy)    ! Compute potential energy

! open the output file
open(unit = 10, file = Trim(OutputFile), IOSTAT = ierror)
if (ierror.EQ.0) then
  print*, 'File result.txt successfully opened!...'
endif


! output result in first time step 
Write(10,100) (t-1)*delta_t,&
& x_current(1,1), x_current(1,2), u_current(1,1), u_current(1,2), kinetic_energy(1), potential_energy(1),& 
& x_current(2,1), x_current(2,2), u_current(2,1), u_current(2,2), kinetic_energy(2), potential_energy(2),&
& x_current(3,1), x_current(3,2), u_current(3,1), u_current(3,2), kinetic_energy(3), potential_energy(3)
100 format(1x,19ES24.15,/)

! Output simulation setup parameters
call OutputParameters(x_current,u_current,kinetic_energy,potential_energy)
Write(*,*) 'Are these parameters property?... Programing Paused, Continue?...'
pause

! new to old
x_old = x_current
x_current = x_new


! Time evolution
do t = 2,N_time

   ! Compute acceleration
     call compute_acceleration(x_current,acceleration)
   ! particles' new pisition
     call compute_position(x_current,x_old,acceleration,x_new)    
   
   ! Output    
       if(mod(t,Output_interval).EQ.0) then   
          ! compute particle's celocity   
          call compute_velocity(x_new,x_old,u_current)
          ! compute particle's potential energy
          call compute_potnetialEnergy(x_current,potential_energy)
          ! compute particles' kinetic energy
          call compute_kineticEnergy(u_current,kinetic_energy)
          print 200, (t-1)       
          Write(10,100) (t-1)*delta_t,&
&            x_current(1,1), x_current(1,2), u_current(1,1), u_current(1,2), kinetic_energy(1), potential_energy(1),& 
&            x_current(2,1), x_current(2,2), u_current(2,1), u_current(2,2), kinetic_energy(2), potential_energy(2),&
&            x_current(3,1), x_current(3,2), u_current(3,1), u_current(3,2), kinetic_energy(3), potential_energy(3)
       end if
    ! time update  
    x_old = x_current
    x_current = x_new

end do
200 format(1x,'Current Time Step = ',I10)
 close(10)

! Output computing log file
call logfile(t,x_current,u_current,kinetic_energy,potential_energy)

END PROGRAM Main
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------
SUBROUTINE OutputParameters(x_current,u_current,kinetic_energy,potential_energy) 
use ParametersBlock
implicit none

! decline variables
real*8, INTENT(IN)::  x_current(1:N_particles,1:N_dimension), u_current(1:N_particles,1:N_dimension)          ! Particles' position and velocity 
real*8, INTENT(IN)::  kinetic_energy(1:N_particles), potential_energy(1:N_particles)                          ! Particles' kinetic energy and potential energy

! file name
character(len=20), parameter:: FileName1='SimulationParameters.txt'

! Print parameters in screen 
Write(*,100) N_particles,N_dimension
100 format(1x,'Particles = 'I10,', Problem Dimension = ',I10,/)
Write(*,101) N_time, delta_t, Output_interval
101 format(1x,'Total time step = 'I10,', Time step size = ',ES15.5,', Data output interval = ',I10,/)
Write(*,102)
102 format(1x,'Position     Velocity      Kinetic Energy    Potential Energy',/)
Write(*,103) x_current(1,1), x_current(1,2), u_current(1,1), u_current(1,2), kinetic_energy(1), potential_energy(1)
103 format(1x,', x1 = ',ES15.5,', y1 = ',ES15.5,', u1 = ',ES15.5,', v1 = ',ES15.5,', k1 = ',ES15.5,', p1 = ',ES15.5,/)
Write(*,104) x_current(2,1), x_current(2,2), u_current(2,1), u_current(2,2), kinetic_energy(2), potential_energy(2)
104 format(1x,', x2 = ',ES15.5,', y2 = ',ES15.5,', u2 = ',ES15.5,', v2 = ',ES15.5,', k2 = ',ES15.5,', p2 = ',ES15.5,/)
Write(*,105) x_current(3,1), x_current(3,2), u_current(3,1), u_current(3,2), kinetic_energy(3), potential_energy(3)
105 format(1x,', x3 = ',ES15.5,', y3 = ',ES15.5,', u3 = ',ES15.5,', v3 = ',ES15.5,', k3 = ',ES15.5,', p3 = ',ES15.5,/)
    

! Output to file
Open(unit=20,file=Trim(FileName1))
Write(20,100) N_particles,N_dimension
Write(20,101) N_time, delta_t, Output_interval
Write(20,102)
Write(20,103) x_current(1,1), x_current(1,2), u_current(1,1), u_current(1,2), kinetic_energy(1), potential_energy(1)
Write(20,104) x_current(2,1), x_current(2,2), u_current(2,1), u_current(2,2), kinetic_energy(2), potential_energy(2)
Write(20,105) x_current(3,1), x_current(3,2), u_current(3,1), u_current(3,2), kinetic_energy(3), potential_energy(3)
 close(20)


END SUBROUTINE OutputParameters
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------





!------------------------------------------------------------------------------------------------
! Compute force
SUBROUTINE compute_acceleration(x_current,acceleration)
use ParametersBlock
implicit none

! decline variables
real*8, INTENT(IN)::  x_current(1:N_particles,1:N_dimension)                                     ! Particles' current position
real*8, INTENT(OUT):: acceleration(1:N_particles,1:N_dimension)                                  ! Particles' acceleration

! Temporary varibles
real*8 x_r,y_r,r,f
! real*8 x_r,y_r,z_r,r,f     ! use this for 3-D 

! index
integer i,j,k

! recation with external potential field
do i = 1, N_particles
  do k = 1, N_dimension
      acceleration(i,k) = -x_current(i,k)
  end do
end do

! reaction between particles

do i = 1, N_particles - 1
     do j = i+1, N_particles
        x_r = x_current(i,1) - x_current(j,1)                     ! compute vector x_ij
        y_r = x_current(i,2) - x_current(j,2)                     ! compute vector y_ij
!       z_r = x_current(i,3) - x_current(j,3)                     ! compute vector z_ij      
        r = sqrt( x_r**2 + y_r**2 )                               ! compute norm of vector r_ij
!       r = sqrt( x_r**2 + y_r**2 + z_r**2)                       ! Three dimension
        f = 24.d0 * ( 2.d0*(1.d0/r)**13 - (1.d0/r)**7 )           ! compute force scalar due to LJ 6-12 Potential
        acceleration(i,1) = acceleration(i,1) + x_r / r * f       ! compute acceleration in x direction
        acceleration(i,2) = acceleration(i,2) + y_r / r * f       ! compute acceleration in y direction
        acceleration(j,1) = acceleration(j,1) - x_r / r * f       ! compute acceleration in x direction  (inverse force)
        acceleration(j,2) = acceleration(j,2) - y_r / r * f       ! compute acceleration in y direction  (inverse force)
     end do
end do


END SUBROUTINE compute_acceleration
!------------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------------
! Update particle's Direction
SUBROUTINE compute_position(x_current,x_old,acceleration,x_new)
use ParametersBlock
implicit none

! decline variables
real*8, INTENT(IN):: x_current(1:N_particles,1:N_dimension), x_old(1:N_particles,1:N_dimension)  ! particles' current position and old position 
real*8, INTENT(IN):: acceleration(1:N_particles,1:N_dimension)                                   ! Particles' current acceleration
real*8, INTENT(OUT):: x_new(1:N_particles,1:N_dimension)                                         ! Particles' new direction

! Index
integer i,j

! Verlet Scheme
do i = 1,N_particles
   do j = 1,N_dimension
      x_new(i,j) = 2.d0 * x_current(i,j) - x_old(i,j) + (delta_t**2) * acceleration(i,j) 
   end do
end do


END SUBROUTINE compute_position
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------





!------------------------------------------------------------------------------------------------
SUBROUTINE compute_velocity(x_new,x_old,u_current)
use ParametersBlock
implicit none
! particles' current position and old position 
real*8, INTENT(IN):: x_new(1:N_particles,1:N_dimension), x_old(1:N_particles,1:N_dimension)
real*8, INTENT(OUT):: u_current(1:N_particles,1:N_dimension)                                       ! Particles' new direction

! index
integer i,j

! compute velocity using centered difference scheme
do i = 1, N_particles
  do j = 1, N_dimension
      u_current(i,j) = ( x_new(i,j) - x_old(i,j) ) / (2.d0*delta_t)
  end do
end do


END SUBROUTINE compute_velocity
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------------
SUBROUTINE compute_potnetialEnergy(x_current,potential_energy)
use ParametersBlock
implicit none
real*8, INTENT(IN):: x_current(1:N_particles,1:N_dimension)                                        ! Particles' current position
real*8, INTENT(OUT):: potential_energy(1:N_particles)                                              ! Particles' potential energy

! Index
integer i,j

! temporary varibales
real*8 x_r,y_r,r,phi
! real*8 x_r,y_r,z_r,r,phi     ! 3D 

! external field potential
do i = 1, N_particles
      potential_energy(i) = 0.5d0*( x_current(i,1)**2 + x_current(i,2)**2 ) 
!     Potential_energy(i) = 0.5d0*( x_current(i,1)**2 + x_current(i,2)**2 + x_current(i,3)**2 )      ! 3D   
end do

! Compute potential energy
do i = 1, N_particles - 1
     do j = i+1, N_particles
        x_r = x_current(i,1) - x_current(j,1)                     ! compute vector x_ij
        y_r = x_current(i,2) - x_current(j,2)                     ! compute vector y_ij
!       z_r = x_current(i,3) - x_current(j,3)                     ! compute vector z_ij      
        r = sqrt( x_r**2 + y_r**2 )                               ! compute norm of vector r_ij
!       r = sqrt( x_r**2 + y_r**2 + z_r**2)                       ! Three dimension
        phi = 4.d0 * ( (1.d0/r)**12 - (1.d0/r)**6 )               ! compute LJ 6-12 Potential
        potential_energy(i) = potential_energy(i) + phi
        potential_energy(j) = potential_energy(j) + phi
     end do
end do


END SUBROUTINE compute_potnetialEnergy
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------
SUBROUTINE compute_kineticEnergy(u_current,kinetic_energy)
use ParametersBlock
implicit none
real*8, INTENT(IN):: u_current(1:N_particles,1:N_dimension)                                         ! Particles' current velocity
real*8, INTENT(OUT):: kinetic_energy(1:N_particles)                                                 ! Particles' kinetic energy

! index
integer i,j

! initial
do i = 1,N_particles
    kinetic_energy(i) = 0.d0  
end do

! Sum
do i = 1,N_particles
   do j = 1,N_dimension
      kinetic_energy(i) = kinetic_energy(i) + 0.5d0*(u_current(i,j)**2)  
   end do
end do

END SUBROUTINE compute_kineticEnergy
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------------
SUBROUTINE logfile(t,x_current,u_current,kinetic_energy,potential_energy)
use ParametersBlock
implicit none

! decline variables
integer, INTENT(IN):: t                                                                                        ! Final time step
real*8, INTENT(IN)::  x_current(1:N_particles,1:N_dimension), u_current(1:N_particles,1:N_dimension)         ! Particles' position and velocity 
real*8, INTENT(IN)::  kinetic_energy(1:N_particles), potential_energy(1:N_particles)                       ! Particles' kinetic energy and potential energy

! file name
character(len=20), parameter:: FileName1='Report.log'

! Output to ReportFile
Open(unit=20,file=Trim(FileName1))
Write(20,100) N_particles,N_dimension
Write(20,101) N_time, delta_t, Output_interval
Write(20,102) t, real(t-1)*delta_t
Write(20,103)
Write(20,104) x_current(1,1), x_current(1,2), u_current(1,1), u_current(1,2), kinetic_energy(1), potential_energy(1)
Write(20,105) x_current(2,1), x_current(2,2), u_current(2,1), u_current(2,2), kinetic_energy(2), potential_energy(2)
Write(20,106) x_current(3,1), x_current(3,2), u_current(3,1), u_current(3,2), kinetic_energy(3), potential_energy(3)
100 format(1x,'Particles = 'I10,', Problem Dimension = ',I10,/)
101 format(1x,'Total time step = 'I10,', Time step size = ',ES15.5,', Data output interval = ',I10,/)
102 format(1x,'Final time step= ',I10,', Final time = ',ES15.5,/)    
103 format(1x,'Position     Velocity      Kinetic Energy    Potential Energy',/)
104 format(1x,', x1 = ',ES15.5,', y1 = ',ES15.5,', u1 = ',ES15.5,', v1 = ',ES15.5,', k1 = ',ES15.5,', p1 = ',ES15.5,/)
105 format(1x,', x2 = ',ES15.5,', y2 = ',ES15.5,', u2 = ',ES15.5,', v2 = ',ES15.5,', k2 = ',ES15.5,', p2 = ',ES15.5,/)
106 format(1x,', x3 = ',ES15.5,', y3 = ',ES15.5,', u3 = ',ES15.5,', v3 = ',ES15.5,', k3 = ',ES15.5,', p3 = ',ES15.5,/)
 close(20)
 

END SUBROUTINE logfile
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

