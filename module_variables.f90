module module_variables

!#################################################

!ALLOCATABLE ARRAYS

implicit none

!Velocity model array
real, allocatable, dimension(:,:) :: mod_vel, mod_vs

!Density model array
real, allocatable, dimension(:,:) :: mod_den

!Pressure field arrays
real, allocatable, dimension(:,:) :: P1, P2, P3, P

!Velocity field arrays
real, allocatable, dimension(:,:) :: U,V

!Displacement field arrays
real, allocatable, dimension(:,:) :: Xx,Zz,Tt

!Inverse density array
real, allocatable, dimension(:,:) :: b

!Elastic Properties
real, allocatable, dimension(:,:) :: K_modulus, L, M

!Source array
real, allocatable, dimension(:) :: sg_source

!Cerjan factors array
real, allocatable, dimension(:) :: vetor_cerjan

!Seismogram array
real, allocatable, dimension(:,:) :: seismogram

!Transit Time Matrix(MTT)
real, allocatable, dimension(:,:) :: AMTT
real, allocatable, dimension(:,:) :: MTT

!RTM Image
real, allocatable, dimension(:,:) :: image

!#################################################

!MODELLING VARIABLES

!Velocity model file
character(len=30) :: file_vp_model, file_vs_model

!Density model file
character(len=30) :: file_density_model

!Seismic source variables
integer :: nfx, nfz, passos_fonte
real :: A, fcorte, fc, t0, tm, fonte, fx, fz 

!Numerical solution stability variables
real :: alfa, beta

!Grid, time & spatial operators variables
integer :: Nx, Nz, h, passos, tmax
real :: x, z, dt, delta, vp_max 

!Cerjan layers variables
integer :: pontos_cerjan
real :: fator_cerjan


!Snapshot variables
character(len=3) :: snap_char
character(len=3) :: image_char
character(len=2) :: aux 
character(len=1) :: string="0"
integer :: nsnap,c_snap, dt_snap,c_image
real :: passos_snap 

!Seismogram variables
real :: receptor_z
integer :: nreceptor_z
real :: passos_seismo

!Number of shots
integer :: n_shots

!Streamer variables
real :: streamer_head
real :: streamer_tail
real :: streamer_length
integer :: nstreamer_head
integer :: nstreamer_tail
integer :: nstreamer_length
real :: spacing
integer :: nspacing

!Counter variables
integer :: i,j,k,n,t 

!Modelling constants
integer :: comp_byte=4
real :: pi=3.14159

!Program constants


!#################################################

end module module_variables



