!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Program name: FiveBodyMD.f90
!
! Program Purpose:
!                  Solve ten-bodies probelm with Nose-Hoover Thermostat in molecule simulation.
!  Version message:
!        Date          Programmer            Description of change
!      ========       ============          =========================
!  1.  36/4/2016      Yang Yang :-)               Source Code
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
! Nose-Hoover thermostat: $zeta(t_\delta t) = zeta(t) + frac{\Sigma m_{i}v_{i}^2 - N_{f}T}{Q} \delta t$    ~ O(delta_t)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!------------------------------------------------------------------------------------------------
! Setup computing parameters 
MODULE ParametersBlock
implicit none

integer, parameter:: N_particles = 10                                                      ! Number of particle
integer, parameter:: N_dimension = 3                                                       ! Problem dimension

!real*8, parameter:: Mass = 1.d0                                                            ! Particles' mass (this term is not necessary!)

integer(kind=8), parameter:: N_time = 50000000000                                          ! Number of time step
real*8, parameter::  delta_t = 1.d-7                                                       ! Time step size  (really important)

real*8, parameter:: Temperature = 150.d0                                                   ! Temeprature of thermostat
!real*8, parameter:: Qm = N_particles * N_dimension * (10000000000*delta_t)**2             ! Thermostat's equvalent mass 
real*8, parameter:: Qm = 1.d0                                                             ! Thermostat's equvalent mass 

integer(kind=8), parameter:: Output_interval = 100000                                      ! Output interval

character(len=20), parameter:: OutputFile = 'result.txt'                                   ! Output file name

real*8, parameter:: pi = 3.1415926535897932384626433832795028841971693993751d0             ! Constant ratio of circumference



! distribution in a ball in three-dimension space
real*8, parameter:: x1_ini = cos(0.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi),&
& y1_ini = sin(0.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi), z1_ini = cos(60.d0/180.d0 * pi)        ! Particle No.1's position

real*8, parameter:: x2_ini = cos(34.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi),&
& y2_ini = sin(34.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi), z2_ini = cos(120.d0/180.d0 * pi)     ! Particle No.2's position

real*8, parameter:: x3_ini = cos(70.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi),&
& y3_ini = sin(70.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi), z3_ini = cos(60.d0/180.d0 * pi)     ! Particle No.3's position

real*8, parameter:: x4_ini = cos(110.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi),&
& y4_ini = sin(110.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi), z4_ini = cos(120.d0/180.d0 * pi)    ! Particle No.4's position

real*8, parameter:: x5_ini = cos(142.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi),&
& y5_ini = sin(142.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi), z5_ini = cos(60.d0/180.d0 * pi)    ! Particle No.5's position

real*8, parameter:: x6_ini = cos(181.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi),&
& y6_ini = sin(181.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi), z6_ini = cos(120.d0/180.d0 * pi)  ! Particle No.6's position

real*8, parameter:: x7_ini = cos(214.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi),&
& y7_ini = sin(214.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi), z7_ini = cos(60.d0/180.d0 * pi)  ! Particle No.7's position

real*8, parameter:: x8_ini = cos(256.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi),&
& y8_ini = sin(256.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi), z8_ini = cos(120.d0/180.d0 * pi)  ! Particle No.8's position

real*8, parameter:: x9_ini = cos(290.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi),&
& y9_ini = sin(290.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi), z9_ini = cos(60.d0/180.d0 * pi)  ! Particle No.9's position

real*8, parameter:: x10_ini = cos(324.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi),&
& y10_ini = sin(324.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi), z10_ini = cos(120.d0/180.d0 * pi) ! Particle No.10's position



!real*8, parameter:: x1_ini = cos(0.d0/180.d0 * pi)*sin(0.d0/180.d0 * pi),&
!& y1_ini = sin(0.d0/180.d0 * pi)*sin(0.d0/180.d0 * pi), z1_ini = cos(0.d0/180.d0 * pi)        ! Particle No.1's position
!
!real*8, parameter:: x2_ini = cos(36.d0/180.d0 * pi)*sin(20.d0/180.d0 * pi),&
!& y2_ini = sin(36.d0/180.d0 * pi)*sin(20.d0/180.d0 * pi), z2_ini = cos(20.d0/180.d0 * pi)     ! Particle No.2's position
!
!real*8, parameter:: x3_ini = cos(72.d0/180.d0 * pi)*sin(40.d0/180.d0 * pi),&
!& y3_ini = sin(72.d0/180.d0 * pi)*sin(40.d0/180.d0 * pi), z3_ini = cos(40.d0/180.d0 * pi)     ! Particle No.3's position
!
!real*8, parameter:: x4_ini = cos(108.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi),&
!& y4_ini = sin(108.d0/180.d0 * pi)*sin(60.d0/180.d0 * pi), z4_ini = cos(60.d0/180.d0 * pi)    ! Particle No.4's position
!
!real*8, parameter:: x5_ini = cos(144.d0/180.d0 * pi)*sin(80.d0/180.d0 * pi),&
!& y5_ini = sin(144.d0/180.d0 * pi)*sin(80.d0/180.d0 * pi), z5_ini = cos(80.d0/180.d0 * pi)    ! Particle No.5's position
!
!real*8, parameter:: x6_ini = cos(180.d0/180.d0 * pi)*sin(100.d0/180.d0 * pi),&
!& y6_ini = sin(180.d0/180.d0 * pi)*sin(100.d0/180.d0 * pi), z6_ini = cos(100.d0/180.d0 * pi)  ! Particle No.6's position
!
!real*8, parameter:: x7_ini = cos(216.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi),&
!& y7_ini = sin(216.d0/180.d0 * pi)*sin(120.d0/180.d0 * pi), z7_ini = cos(120.d0/180.d0 * pi)  ! Particle No.7's position
!
!real*8, parameter:: x8_ini = cos(252.d0/180.d0 * pi)*sin(140.d0/180.d0 * pi),&
!& y8_ini = sin(252.d0/180.d0 * pi)*sin(140.d0/180.d0 * pi), z8_ini = cos(140.d0/180.d0 * pi)  ! Particle No.8's position
!
!real*8, parameter:: x9_ini = cos(288.d0/180.d0 * pi)*sin(160.d0/180.d0 * pi),&
!& y9_ini = sin(288.d0/180.d0 * pi)*sin(160.d0/180.d0 * pi), z9_ini = cos(160.d0/180.d0 * pi)  ! Particle No.9's position
!
!real*8, parameter:: x10_ini = cos(324.d0/180.d0 * pi)*sin(180.d0/180.d0 * pi),&
!& y10_ini = sin(324.d0/180.d0 * pi)*sin(180.d0/180.d0 * pi), z10_ini = cos(180.d0/180.d0 * pi)  ! Particle No.10's position



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
real*8 acceleration(1:N_particles,1:N_dimension), acceleration_int(1:N_particles,1:N_dimension)                         ! Particles' acceleration
real*8 potential_energy(1:N_particles), kinetic_energy(1:N_particles)                                                   ! particle's Potential energy and Kinetic energy
real*8 potential_energy_total, kinetic_energy_total                                                                     ! system's total potnetial and kinetic energy
real*8 zeta_new, zeta_current                                                                                           ! damping coefficient

! Index
integer i,j,ierror
integer(kind=8):: t = 0

! Initial all varibales

 
do i = 1, N_particles
  do j = 1, N_dimension
      x_new(i,j) = 0.d0
      x_current(i,j) = 0.d0
      x_old(i,j) = 0.d0
      u_current(i,j) = 0.d0
      acceleration(i,j) = 0.d0
      acceleration_int(i,j) = 0.d0
  end do
end do

do i = 1, N_particles
    potential_energy(i) = 0.d0
    kinetic_energy(i) = 0.d0 
end do

potential_energy_total = 0.d0
kinetic_energy_total = 0.d0
zeta_new = 0.d0
zeta_current = 0.d0   


! Initial particles' state
! particles' initial position
x_current(1,1) = x1_ini
x_current(1,2) = y1_ini
x_current(1,3) = z1_ini

x_current(2,1) = x2_ini
x_current(2,2) = y2_ini
x_current(2,3) = z2_ini

x_current(3,1) = x3_ini
x_current(3,2) = y3_ini
x_current(3,3) = z3_ini

x_current(4,1) = x4_ini
x_current(4,2) = y4_ini
x_current(4,3) = z4_ini

x_current(5,1) = x5_ini
x_current(5,2) = y5_ini
x_current(5,3) = z5_ini

x_current(6,1) = x6_ini
x_current(6,2) = y6_ini
x_current(6,3) = z6_ini

x_current(7,1) = x7_ini
x_current(7,2) = y7_ini
x_current(7,3) = z7_ini

x_current(8,1) = x8_ini
x_current(8,2) = y8_ini
x_current(8,3) = z8_ini

x_current(9,1) = x9_ini
x_current(9,2) = y9_ini
x_current(9,3) = z9_ini

x_current(10,1) = x10_ini
x_current(10,2) = y10_ini
x_current(10,3) = z10_ini


call compute_acceleration(x_current,u_current,zeta_current,acceleration,acceleration_int)   ! Compute acceleration at first time step

! time step evolution
do i = 1,N_particles
   do j = 1,N_dimension
      x_new(i,j) = x_current(i,j) + delta_t*u_current(i,j) + 0.5d0*(delta_t**2)*acceleration(i,j) 
   end do
end do

! compute potential energy
call compute_potnetialEnergy(x_current,potential_energy,potential_energy_total)          ! Compute potential energy

! compute kinetic energy 
call compute_kineticEnergy(u_current,kinetic_energy,kinetic_energy_total)                ! Compute kinetic energy
 
! compute damping coefficient
call compute_zeta(u_current,zeta_current,zeta_new)                                       ! Compute damping coefficient

! open the output file
open(unit = 10, file = Trim(OutputFile), IOSTAT = ierror)
if (ierror.EQ.0) then
  print*, 'File result.txt successfully opened!...'
endif


! output result in first time step 
Write(10,100) (t-1)*delta_t,&
& x_current(1,1), x_current(1,2), x_current(1,3), u_current(1,1), u_current(1,2), u_current(1,3),&
&   acceleration_int(1,1), acceleration_int(1,2), acceleration_int(1,3), kinetic_energy(1), potential_energy(1),& 
& x_current(2,1), x_current(2,2), x_current(2,3), u_current(2,1), u_current(2,2), u_current(2,3),&
&   acceleration_int(2,1), acceleration_int(2,2), acceleration_int(2,3), kinetic_energy(2), potential_energy(2),&
& x_current(3,1), x_current(3,2), x_current(3,3), u_current(3,1), u_current(3,2), u_current(3,3),&
&   acceleration_int(3,1), acceleration_int(3,2), acceleration_int(3,3), kinetic_energy(3), potential_energy(3),&
& x_current(4,1), x_current(4,2), x_current(4,3), u_current(4,1), u_current(4,2), u_current(4,3),&
&   acceleration_int(4,1), acceleration_int(4,2), acceleration_int(4,3), kinetic_energy(4), potential_energy(4),&
& x_current(5,1), x_current(5,2), x_current(5,3), u_current(5,1), u_current(5,2), u_current(5,3),&
&   acceleration_int(5,1), acceleration_int(5,2), acceleration_int(5,3), kinetic_energy(5), potential_energy(5),&
& x_current(6,1), x_current(6,2), x_current(6,3), u_current(6,1), u_current(6,2), u_current(6,3),&
&   acceleration_int(6,1), acceleration_int(6,2), acceleration_int(6,3), kinetic_energy(6), potential_energy(6),& 
& x_current(7,1), x_current(7,2), x_current(7,3), u_current(7,1), u_current(7,2), u_current(7,3),&
&   acceleration_int(7,1), acceleration_int(7,2), acceleration_int(7,3), kinetic_energy(7), potential_energy(7),&
& x_current(8,1), x_current(8,2), x_current(8,3), u_current(8,1), u_current(8,2), u_current(8,3),&
&   acceleration_int(8,1), acceleration_int(8,2), acceleration_int(8,3), kinetic_energy(8), potential_energy(8),&
& x_current(9,1), x_current(9,2), x_current(9,3), u_current(9,1), u_current(9,2), u_current(9,3),&
&   acceleration_int(9,1), acceleration_int(9,2), acceleration_int(9,3), kinetic_energy(9), potential_energy(9),&
& x_current(10,1), x_current(10,2), x_current(10,3), u_current(10,1), u_current(10,2), u_current(10,3),&
&   acceleration_int(10,1), acceleration_int(10,2), acceleration_int(10,3), kinetic_energy(10), potential_energy(10)

100 format(1x,111ES24.15,/)

! Output simulation setup parameters
call OutputParameters(x_current,u_current,kinetic_energy,potential_energy,kinetic_energy_total, potential_energy_total,zeta_current)
Write(*,*) 'Are these parameters property?... Programing Paused, Continue?...'
pause

! new to old
x_old = x_current
x_current = x_new
zeta_current = zeta_new

! Time evolution
do t = 2,N_time

   ! Compute acceleration
     call compute_acceleration(x_current,u_current,zeta_current,acceleration,acceleration_int)
   ! particles' new pisition
     call compute_position(x_current,x_old,acceleration,x_new)    
   ! compute particle's velocity    
     call compute_velocity(x_new,x_old,u_current)
   ! compute damping coefficient  
     call compute_zeta(u_current,zeta_current,zeta_new)                                                    
   ! Output    
       if(mod(t,Output_interval).EQ.0) then   
          ! compute particle's potential energy
          call compute_potnetialEnergy(x_current,potential_energy,potential_energy_total)          
          ! compute particles' kinetic energy
          call compute_kineticEnergy(u_current,kinetic_energy,kinetic_energy_total)                
          print 200, (t-1)
          print 201, kinetic_energy_total, potential_energy_total,  kinetic_energy_total/(0.5*N_particles*N_dimension)
          print 202, zeta_current
          print *,  ''     
          Write(10,100) (t-1)*delta_t,&
&         x_current(1,1), x_current(1,2), x_current(1,3), u_current(1,1), u_current(1,2), u_current(1,3),&
&           acceleration_int(1,1), acceleration_int(1,2), acceleration_int(1,3), kinetic_energy(1), potential_energy(1),& 
&         x_current(2,1), x_current(2,2), x_current(2,3), u_current(2,1), u_current(2,2), u_current(2,3),&
&           acceleration_int(2,1), acceleration_int(2,2), acceleration_int(2,3), kinetic_energy(2), potential_energy(2),&
&         x_current(3,1), x_current(3,2), x_current(3,3), u_current(3,1), u_current(3,2), u_current(3,3),&
&           acceleration_int(3,1), acceleration_int(3,2), acceleration_int(3,3), kinetic_energy(3), potential_energy(3),&
&         x_current(4,1), x_current(4,2), x_current(4,3), u_current(4,1), u_current(4,2), u_current(4,3),&
&           acceleration_int(4,1), acceleration_int(4,2), acceleration_int(4,3), kinetic_energy(4), potential_energy(4),&
&         x_current(5,1), x_current(5,2), x_current(5,3), u_current(5,1), u_current(5,2), u_current(5,3),&
&           acceleration_int(5,1), acceleration_int(5,2), acceleration_int(5,3), kinetic_energy(5), potential_energy(5),&
&         x_current(6,1), x_current(6,2), x_current(6,3), u_current(6,1), u_current(6,2), u_current(6,3),&
&           acceleration_int(6,1), acceleration_int(6,2), acceleration_int(6,3), kinetic_energy(6), potential_energy(6),& 
&         x_current(7,1), x_current(7,2), x_current(7,3), u_current(7,1), u_current(7,2), u_current(7,3),&
&           acceleration_int(7,1), acceleration_int(7,2), acceleration_int(7,3), kinetic_energy(7), potential_energy(7),&
&         x_current(8,1), x_current(8,2), x_current(8,3), u_current(8,1), u_current(8,2), u_current(8,3),&
&           acceleration_int(8,1), acceleration_int(8,2), acceleration_int(8,3), kinetic_energy(8), potential_energy(8),&
&         x_current(9,1), x_current(9,2), x_current(9,3), u_current(9,1), u_current(9,2), u_current(9,3),&
&           acceleration_int(9,1), acceleration_int(9,2), acceleration_int(9,3), kinetic_energy(9), potential_energy(9),&
&         x_current(10,1), x_current(10,2), x_current(10,3), u_current(10,1), u_current(10,2), u_current(10,3),&
&           acceleration_int(10,1), acceleration_int(10,2), acceleration_int(10,3), kinetic_energy(10), potential_energy(10)
       !pause
       end if
    ! time update
    x_old = x_current
    x_current = x_new
    zeta_current = zeta_new 
!    zeta_current = 1.d0 

end do


200 format(1x,'Current Time Step = ',I20)
201 format(1x,'SystemKineticEnergy = ',ES15.5,', SystemPotentialEnergy=',ES15.5,', SystemTemperature = 'ES15.5)
202 format(1x,'DampingCoefficient = ',ES15.5,/)

 close(10)

! Output computing log file
call logfile(t,x_current,u_current,kinetic_energy,potential_energy,potential_energy_total, kinetic_energy_total,zeta_current)

END PROGRAM Main
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------
SUBROUTINE OutputParameters(x_current,u_current,kinetic_energy,potential_energy,kinetic_energy_total, potential_energy_total,zeta)
use ParametersBlock
implicit none

! decline variables
real*8, INTENT(IN):: x_current(1:N_particles,1:N_dimension), u_current(1:N_particles,1:N_dimension)          ! Particles' position and velocity 
real*8, INTENT(IN):: kinetic_energy(1:N_particles), potential_energy(1:N_particles)                          ! Particles' kinetic energy and potential energy
real*8, INTENT(IN):: potential_energy_total, kinetic_energy_total                                            ! system's total potnetial and kinetic energy
real*8, INTENT(IN):: zeta                                                                                    ! Damping coefficient

! file name
character(len=20), parameter:: FileName1='SimulationParameters.txt'

! Print parameters in screen 
Write(*,100) N_particles,N_dimension
100 format(1x,'Particles = 'I10,', Problem Dimension = ',I10,/)

Write(*,101) N_time, delta_t, Output_interval
101 format(1x,'Total time step = 'I20,', Time step size = ',ES15.5,', Data output interval = ',I20,/)

Write(*,102)
102 format(1x,'Position     Velocity      Kinetic Energy    Potential Energy',/)

Write(*,103) x_current(1,1), x_current(1,2), x_current(1,3), u_current(1,1), u_current(1,2), u_current(1,3),&
&            kinetic_energy(1), potential_energy(1)
103 format(1x,', x1 = ',ES15.5,', y1 = ',ES15.5,', z1 = ',ES15.5,', u1 = ',ES15.5,', v1 = ',ES15.5,', w1 = ',ES15.5,&
&          ', k1 = ',ES15.5,', p1 = ',ES15.5,/)

Write(*,104) x_current(2,1), x_current(2,2), x_current(2,3), u_current(2,1), u_current(2,2), u_current(2,3),&
&            kinetic_energy(2), potential_energy(2)
104 format(1x,', x2 = ',ES15.5,', y2 = ',ES15.5,', z2 = ',ES15.5,', u2 = ',ES15.5,', v2 = ',ES15.5,', w2 = ',ES15.5,&
&          ', k2 = ',ES15.5,', p2 = ',ES15.5,/)

Write(*,105) x_current(3,1), x_current(3,2), x_current(3,3), u_current(3,1), u_current(3,2), u_current(3,3),&
&            kinetic_energy(3), potential_energy(3)
105 format(1x,', x3 = ',ES15.5,', y3 = ',ES15.5,', z3 = ',ES15.5,', u3 = ',ES15.5,', v3 = ',ES15.5,', w3 = ',ES15.5,&
&          ', k3 = ',ES15.5,', p3 = ',ES15.5,/)

Write(*,106) x_current(4,1), x_current(4,2), x_current(4,3), u_current(4,1), u_current(4,2), u_current(4,3),&
&            kinetic_energy(4), potential_energy(4)
106 format(1x,', x4 = ',ES15.5,', y4 = ',ES15.5,', z4 = ',ES15.5,', u4 = ',ES15.5,', v4 = ',ES15.5,', w4 = ',ES15.5,&
&          ', k4 = ',ES15.5,', p4 = ',ES15.5,/)

Write(*,107) x_current(5,1), x_current(5,2), x_current(5,3), u_current(5,1), u_current(5,2), u_current(5,3),&
&            kinetic_energy(5), potential_energy(5)
107 format(1x,', x5 = ',ES15.5,', y5 = ',ES15.5,', z5 = ',ES15.5,', u5 = ',ES15.5,', v5 = ',ES15.5,', w5 = ',ES15.5,&
&          ', k5 = ',ES15.5,', p5 = ',ES15.5,/)

Write(*,108) x_current(6,1), x_current(6,2), x_current(6,3), u_current(6,1), u_current(6,2), u_current(6,3),&
&            kinetic_energy(6), potential_energy(6)
108 format(1x,', x6 = ',ES15.5,', y6 = ',ES15.5,', z6 = ',ES15.5,', u6 = ',ES15.5,', v6 = ',ES15.5,', w6 = ',ES15.5,&
&          ', k6 = ',ES15.5,', p6 = ',ES15.5,/)

Write(*,109) x_current(7,1), x_current(7,2), x_current(7,3), u_current(7,1), u_current(7,2), u_current(7,3),&
&            kinetic_energy(7), potential_energy(7)
109 format(1x,', x7 = ',ES15.5,', y7 = ',ES15.5,', z7 = ',ES15.5,', u7 = ',ES15.5,', v7 = ',ES15.5,', w7 = ',ES15.5,&
&          ', k7 = ',ES15.5,', p7 = ',ES15.5,/)

Write(*,110) x_current(8,1), x_current(8,2), x_current(8,3), u_current(8,1), u_current(8,2), u_current(8,3),&
&            kinetic_energy(8), potential_energy(8)
110 format(1x,', x8 = ',ES15.5,', y8 = ',ES15.5,', z8 = ',ES15.5,', u8 = ',ES15.5,', v8 = ',ES15.5,', w8 = ',ES15.5,&
           ', k8 = ',ES15.5,', p8 = ',ES15.5,/)

Write(*,111) x_current(9,1), x_current(9,2), x_current(9,3), u_current(9,1), u_current(9,2), u_current(9,3),&
&            kinetic_energy(9), potential_energy(9)
111 format(1x,', x9 = ',ES15.5,', y9 = ',ES15.5,', z9 = ',ES15.5,', u9 = ',ES15.5,', v9 = ',ES15.5,', w9 = ',ES15.5,&
&          ', k9 = ',ES15.5,', p9 = ',ES15.5,/)

Write(*,112) x_current(10,1), x_current(10,2), x_current(10,3), u_current(10,1), u_current(10,2), u_current(10,3),&
&            kinetic_energy(10), potential_energy(10)
112 format(1x,', x10 = ',ES15.5,', y10 = ',ES15.5,', z10 = ',ES15.5,', u10 = ',ES15.5,', v10 = ',ES15.5,', w10 = ',ES15.5,&
&          ', k10 = ',ES15.5,', p10 = ',ES15.5,/)

Write(*,113) Temperature, Qm
113 format(1x,'ThermostatTemperatureSetup = ',ES15.5,', Qm = ',ES15.5,/)

Write(*,114) zeta
114 format(1x,'Damping coefficient = ',ES15.5,/)

Write(*,115) potential_energy_total, kinetic_energy_total
115 format(1x,'SystemPotentialEnergy = ',ES15.5,', SystemKineticEnergy = ',ES15.5/)

    


! Output to file
Open(unit=20,file=Trim(FileName1))
Write(20,100) N_particles,N_dimension

Write(20,101) N_time, delta_t, Output_interval

Write(20,102)

Write(20,103) x_current(1,1), x_current(1,2), x_current(1,3), u_current(1,1), u_current(1,2), u_current(1,3),&
&            kinetic_energy(1), potential_energy(1)

Write(20,104) x_current(2,1), x_current(2,2), x_current(2,3), u_current(2,1), u_current(2,2), u_current(2,3),&
&            kinetic_energy(2), potential_energy(2)

Write(20,105) x_current(3,1), x_current(3,2), x_current(3,3), u_current(3,1), u_current(3,2), u_current(3,3),&
&            kinetic_energy(3), potential_energy(3)

Write(20,106) x_current(4,1), x_current(4,2), x_current(4,3), u_current(4,1), u_current(4,2), u_current(4,3),&
&            kinetic_energy(4), potential_energy(4)

Write(20,107) x_current(5,1), x_current(5,2), x_current(5,3), u_current(5,1), u_current(5,2), u_current(5,3),&
&            kinetic_energy(5), potential_energy(5)

Write(20,108) x_current(6,1), x_current(6,2), x_current(6,3), u_current(6,1), u_current(6,2), u_current(6,3),&
&            kinetic_energy(6), potential_energy(6)

Write(20,109) x_current(7,1), x_current(7,2), x_current(7,3), u_current(7,1), u_current(7,2), u_current(7,3),&
&            kinetic_energy(7), potential_energy(7)

Write(20,110) x_current(8,1), x_current(8,2), x_current(8,3), u_current(8,1), u_current(8,2), u_current(8,3),&
&            kinetic_energy(8), potential_energy(8)

Write(20,111) x_current(9,1), x_current(9,2), x_current(9,3), u_current(9,1), u_current(9,2), u_current(9,3),&
&            kinetic_energy(9), potential_energy(9)

Write(20,112) x_current(10,1), x_current(10,2), x_current(10,3), u_current(10,1), u_current(10,2), u_current(10,3),&
&            kinetic_energy(10), potential_energy(10)

Write(20,113) Temperature, Qm

Write(20,114) zeta

Write(20,115) potential_energy_total, kinetic_energy_total


END SUBROUTINE OutputParameters
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------





!------------------------------------------------------------------------------------------------
! Compute force
SUBROUTINE compute_acceleration(x_current,u_current,zeta_current,acceleration,acceleration_int)
use ParametersBlock
implicit none

! decline variables
real*8, INTENT(IN)::  x_current(1:N_particles,1:N_dimension)                                     ! Particles' current position
real*8, INTENT(IN)::  u_current(1:N_particles,1:N_dimension)                                     ! Particles' current velocity
real*8, INTENT(IN)::  zeta_current                                                               ! Damping coefficient
real*8, INTENT(OUT):: acceleration(1:N_particles,1:N_dimension)                                  ! Particles' total acceleration
real*8, INTENT(OUT):: acceleration_int(1:N_particles,1:N_dimension)                              ! Particles' conservation reaction acceleration

! Temporary varibles
!real*8 x_r,y_r,r,f         ! use this for 2-D
real*8 x_r,y_r,z_r,r,f     ! use this for 3-D 

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
        z_r = x_current(i,3) - x_current(j,3)                     ! compute vector z_ij      
!        r = sqrt( x_r**2 + y_r**2 )                               ! compute norm of vector r_ij (2-D)
        r = sqrt( x_r**2 + y_r**2 + z_r**2)                       ! compute norm of vector r_ij (3-D)
        f = 24.d0 * ( 2.d0*(1.d0/r)**13 - (1.d0/r)**7 )           ! compute force scalar due to LJ 6-12 Potential
        acceleration(i,1) = acceleration(i,1) + x_r / r * f       ! compute acceleration in x direction
        acceleration(i,2) = acceleration(i,2) + y_r / r * f       ! compute acceleration in y direction
        acceleration(i,3) = acceleration(i,3) + z_r / r * f       ! compute acceleration in z direction
        acceleration(j,1) = acceleration(j,1) - x_r / r * f       ! compute acceleration in x direction  (inverse force)
        acceleration(j,2) = acceleration(j,2) - y_r / r * f       ! compute acceleration in y direction  (inverse force)
        acceleration(j,3) = acceleration(j,3) - z_r / r * f       ! compute acceleration in z direction  (inverse force)
     end do
end do

acceleration_int = acceleration                                   ! Particles' conservation reaction acceleration                                   

! Nose-hoover thermostat damping coefficient
do i = 1, N_particles
  do k = 1, N_dimension
      acceleration(i,k) = acceleration(i,k) - zeta_current*u_current(i,k)
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
SUBROUTINE compute_zeta(u_current,zeta_current,zeta_new)
use ParametersBlock
implicit none

! velocity imput
real*8, INTENT(IN):: u_current(1:N_particles,1:N_dimension)

! damping coefficient 
real*8, INTENT(IN):: zeta_current
real*8, INTENT(OUT):: zeta_new

! temporary variable
real*8 delta_zeta, temp 

! index
integer i,j

! initial
delta_zeta = 0.d0
temp = 0.d0


! Compute zeta_new
do i = 1, N_particles
   do j = 1, N_dimension
      temp = temp + u_current(i,j)**2
   end do 
end do


delta_zeta = ((temp - (N_particles*N_dimension)*Temperature ) / Qm) * delta_t

zeta_new = zeta_current + delta_zeta


!print 100, temp, (N_particles*N_dimension)*Temperature, delta_zeta
!100 format(1x,3ES15.5,/)

END SUBROUTINE compute_zeta
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------




!------------------------------------------------------------------------------------------------
SUBROUTINE compute_potnetialEnergy(x_current,potential_energy,potential_energy_total)
use ParametersBlock
implicit none
real*8, INTENT(IN):: x_current(1:N_particles,1:N_dimension)                                        ! Particles' current position
real*8, INTENT(OUT):: potential_energy(1:N_particles)                                              ! Particles' potential energy
real*8, INTENT(OUT):: potential_energy_total                                                       ! System's total potential energy

! Index
integer i,j

! temporary varibales
! real*8 x_r,y_r,r,phi       ! Use this for 2D
real*8 x_r,y_r,z_r,r,phi     ! Use this for 3D 
do i = 1, N_particles
      potential_energy = 0.d0                                    ! Initial
end do

! system's total potential energy
potential_energy_total = 0.d0                                    ! Initial


! external field potential
do i = 1, N_particles
!    potential_energy(i) = 0.5d0*( x_current(i,1)**2 + x_current(i,2)**2 + x_current(i,3)**2 )      ! 3D
     do j = 1, N_dimension
      potential_energy(i) = potential_energy(i) + 0.5d0*x_current(i,j)**2
     end do   
     potential_energy_total = potential_energy_total + potential_energy(i)
end do

! Compute potential energy
do i = 1, N_particles - 1
     do j = i+1, N_particles
        x_r = x_current(i,1) - x_current(j,1)                     ! compute vector x_ij
        y_r = x_current(i,2) - x_current(j,2)                     ! compute vector y_ij
        z_r = x_current(i,3) - x_current(j,3)                     ! compute vector z_ij      
!        r = sqrt( x_r**2 + y_r**2 )                               ! compute norm of vector r_ij (2-D)
        r = sqrt( x_r**2 + y_r**2 + z_r**2)                       ! compute norm of vector r_ij (3-D)
        phi = 4.d0 * ( (1.d0/r)**12 - (1.d0/r)**6 )               ! compute LJ 6-12 Potential
        potential_energy(i) = potential_energy(i) + phi
        potential_energy(j) = potential_energy(j) + phi
     end do
     potential_energy_total = potential_energy_total + potential_energy(i)
end do


END SUBROUTINE compute_potnetialEnergy
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------
SUBROUTINE compute_kineticEnergy(u_current,kinetic_energy,kinetic_energy_total)
use ParametersBlock
implicit none
real*8, INTENT(IN):: u_current(1:N_particles,1:N_dimension)                                         ! Particles' current velocity
real*8, INTENT(OUT):: kinetic_energy(1:N_particles)                                                 ! Particles' kinetic energy
real*8, INTENT(OUT):: kinetic_energy_total                                                          ! System's total kinetic energy

! index
integer i,j

! initial
do i = 1,N_particles
    kinetic_energy(i) = 0.d0  
end do

kinetic_energy_total = 0.d0


! Sum
do i = 1,N_particles
   do j = 1,N_dimension
      kinetic_energy(i) = kinetic_energy(i) + 0.5d0*(u_current(i,j)**2)  
   end do
   kinetic_energy_total = kinetic_energy_total + kinetic_energy(i)
end do

END SUBROUTINE compute_kineticEnergy
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------



!------------------------------------------------------------------------------------------------
SUBROUTINE logfile(t,x_current,u_current,kinetic_energy,potential_energy,potential_energy_total, kinetic_energy_total,zeta)
use ParametersBlock
implicit none

! decline variables
integer(kind=8), INTENT(IN):: t                                                                                ! Final time step
real*8, INTENT(IN)::  x_current(1:N_particles,1:N_dimension), u_current(1:N_particles,1:N_dimension)           ! Particles' position and velocity 
real*8, INTENT(IN)::  kinetic_energy(1:N_particles), potential_energy(1:N_particles)                           ! Particles' kinetic energy and potential energy
real*8, INTENT(IN)::  potential_energy_total, kinetic_energy_total                                             ! system's total potnetial and kinetic energy
real*8, INTENT(IN)::  zeta                                                                                     ! Damping coefficient

! file name
character(len=20), parameter:: FileName1='Report.log'

! Output to ReportFile
Open(unit=20,file=Trim(FileName1))
Write(20,100) N_particles,N_dimension
Write(20,101) N_time, delta_t, Output_interval
Write(20,102) t, real(t-1)*delta_t
Write(20,103)
Write(20,104) x_current(1,1), x_current(1,2), x_current(1,3), u_current(1,1), u_current(1,2), u_current(1,3),&
&            kinetic_energy(1), potential_energy(1)
Write(20,105) x_current(2,1), x_current(2,2), x_current(2,3), u_current(2,1), u_current(2,2), u_current(2,3),&
&            kinetic_energy(2), potential_energy(2)
Write(20,106) x_current(3,1), x_current(3,2), x_current(3,3), u_current(3,1), u_current(3,2), u_current(3,3),&
&            kinetic_energy(3), potential_energy(3)
Write(20,107) x_current(4,1), x_current(4,2), x_current(4,3), u_current(4,1), u_current(4,2), u_current(4,3),&
&            kinetic_energy(4), potential_energy(4)
Write(20,108) x_current(5,1), x_current(5,2), x_current(5,3), u_current(5,1), u_current(5,2), u_current(5,3),&
&            kinetic_energy(5), potential_energy(5)
Write(20,109) x_current(6,1), x_current(6,2), x_current(6,3), u_current(6,1), u_current(6,2), u_current(6,3),&
&            kinetic_energy(6), potential_energy(6)
Write(20,110) x_current(7,1), x_current(7,2), x_current(7,3), u_current(7,1), u_current(7,2), u_current(7,3),&
&            kinetic_energy(7), potential_energy(7)
Write(20,111) x_current(8,1), x_current(8,2), x_current(8,3), u_current(8,1), u_current(8,2), u_current(8,3),&
&            kinetic_energy(8), potential_energy(8)
Write(20,112) x_current(9,1), x_current(9,2), x_current(9,3), u_current(9,1), u_current(9,2), u_current(9,3),&
&            kinetic_energy(9), potential_energy(9)
Write(20,113) x_current(10,1), x_current(10,2), x_current(10,3), u_current(10,1), u_current(10,2), u_current(10,3),&
&            kinetic_energy(10), potential_energy(10)
Write(20,114) Temperature, Qm
Write(20,115) zeta
Write(20,116) potential_energy_total, kinetic_energy_total




100 format(1x,'Particles = 'I10,', Problem Dimension = ',I10,/)
101 format(1x,'Total time step = 'I20,', Time step size = ',ES15.5,', Data output interval = ',I20,/)
102 format(1x,'Final time step= ',I20,', Final time = ',ES15.5,/)    
103 format(1x,'Position     Velocity      Kinetic Energy    Potential Energy',/)
104 format(1x,', x1 = ',ES15.5,', y1 = ',ES15.5,', z1 = ',ES15.5,', u1 = ',ES15.5,', v1 = ',ES15.5,', w1 = ',ES15.5,&
&          ', k1 = ',ES15.5,', p1 = ',ES15.5,/)
105 format(1x,', x2 = ',ES15.5,', y2 = ',ES15.5,', z2 = ',ES15.5,', u2 = ',ES15.5,', v2 = ',ES15.5,', w2 = ',ES15.5,&
&          ', k2 = ',ES15.5,', p2 = ',ES15.5,/)
106 format(1x,', x3 = ',ES15.5,', y3 = ',ES15.5,', z3 = ',ES15.5,', u3 = ',ES15.5,', v3 = ',ES15.5,', w3 = ',ES15.5,&
&          ', k3 = ',ES15.5,', p3 = ',ES15.5,/)
107 format(1x,', x4 = ',ES15.5,', y4 = ',ES15.5,', z4 = ',ES15.5,', u4 = ',ES15.5,', v4 = ',ES15.5,', w4 = ',ES15.5,&
&          ', k4 = ',ES15.5,', p4 = ',ES15.5,/)
108 format(1x,', x5 = ',ES15.5,', y5 = ',ES15.5,', z5 = ',ES15.5,', u5 = ',ES15.5,', v5 = ',ES15.5,', w5 = ',ES15.5,&
&          ', k5 = ',ES15.5,', p5 = ',ES15.5,/)
109 format(1x,', x6 = ',ES15.5,', y6 = ',ES15.5,', z6 = ',ES15.5,', u6 = ',ES15.5,', v6 = ',ES15.5,', w6 = ',ES15.5,&
&          ', k6 = ',ES15.5,', p6 = ',ES15.5,/)
110 format(1x,', x7 = ',ES15.5,', y7 = ',ES15.5,', z7 = ',ES15.5,', u7 = ',ES15.5,', v7 = ',ES15.5,', w7 = ',ES15.5,&
&          ', k7 = ',ES15.5,', p7 = ',ES15.5,/)
111 format(1x,', x8 = ',ES15.5,', y8 = ',ES15.5,', z8 = ',ES15.5,', u8 = ',ES15.5,', v8 = ',ES15.5,', w8 = ',ES15.5,&
           ', k8 = ',ES15.5,', p8 = ',ES15.5,/)
112 format(1x,', x9 = ',ES15.5,', y9 = ',ES15.5,', z9 = ',ES15.5,', u9 = ',ES15.5,', v9 = ',ES15.5,', w9 = ',ES15.5,&
&          ', k9 = ',ES15.5,', p9 = ',ES15.5,/)
113 format(1x,', x10 = ',ES15.5,', y10 = ',ES15.5,', z10 = ',ES15.5,', u10 = ',ES15.5,', v10 = ',ES15.5,', w10 = ',ES15.5,&
&          ', k10 = ',ES15.5,', p10 = ',ES15.5,/)
114 format(1x,'ThermostatTemperatureSetup = ',ES15.5,', Qm = ',ES15.5,/)
115 format(1x,'Damping coefficient = ',ES15.5,/)
116 format(1x,'SystemPotentialEnergy = ',ES15.5,', SystemKineticEnergy = ',ES15.5/)

 close(20)
 

END SUBROUTINE logfile
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

