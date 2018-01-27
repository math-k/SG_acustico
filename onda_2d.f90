program onda_2d

!##################################################

!Module containing all necessary program variables
use module_variables

implicit none

real :: start_time,final_time
!#################################################

!#################################################

call CPU_TIME(start_time)
call read_verify_input()
call allocate_variables()
call read_verify_velocity_model()
call read_verify_density_model()
call calculate_source_modelling_parameters()

write(*,*)"Number of snapshots:",nsnap
close(50)
!LOOP DO TEMPO
do k=1, passos
    if (k<=passos_fonte) then
        P(nfz,nfx)=P(nfz,nfx)+sg_source(k)
        !P(nfz,nfx)=P(nfz,nfx)+staggered_grid_source(A,pi,fc,k,dt,t0)
    end if
    
    !call acoustic_second_order_staggered_grid_operator_loop()
    call acoustic_fourth_order_staggered_grid_operator_loop()
    call anti_reflection_cerjan_conditions()
    call transit_time_matrix()
    if(mod(k,dt_snap)==0) then
       c_snap=c_snap+1
       !write(*,*) "Gerando snapshot:", c_snap
       write(snap_char,"(i3)") c_snap
       call take_snapshot()
    end if
    

    do z=nstreamer_tail,nstreamer_head,nspacing
        seismogram(k,z)=P(nreceptor_z,z)
    end do
        
end do

write(*,*) "Escrevendo sismograma..."
open(25,file="seismogram.bin",status="replace",access="direct",form="unformatted",recl=Nx*k)
write(25,rec=1) seismogram 
close(25)

write(*,*) "Gerando imagem..."
P=0
U=0
V=0
do k=passos,1,-1

    do j=nstreamer_tail,nstreamer_head,nspacing
        P(nreceptor_z,j)=P(nreceptor_z,j)+seismogram(k,j)
    end do
    
    !call acoustic_second_order_staggered_grid_operator_loop()
    call acoustic_fourth_order_staggered_grid_operator_loop()
    call anti_reflection_cerjan_conditions()
    call RTM()

end do

call save_image()

call CPU_TIME(final_time)
write(*,*) "Tempo de processamento decorrido:",final_time-start_time,"segundos."

contains
!#################################################

subroutine read_verify_input()
implicit none

open (30, file="input.txt")
read(30,*) x, z, h, tmax
read(30,*) file_vp_model
read(30,*) file_density_model
read(30,*) fx, fz, A, fcorte
read(30,*) alfa, beta
read(30,*) fator_cerjan, pontos_cerjan
read(30,*) nsnap
read(30,*) receptor_z
read(30,*) n_shots
read(30,*) streamer_tail
read(30,*) streamer_head
read(30,*) spacing
read(30,*) streamer_length
Nx=nint(x/h)
Nz=nint(z/h)
nfx=nint(fx/h)
nfz=nint(fz/h)
nreceptor_z=nint(receptor_z/h)
nstreamer_tail=nint(streamer_tail/h)+1
nstreamer_head=nint(streamer_head/h)
nspacing=nint(spacing/h)
nstreamer_length=nint(streamer_length/h)
write(*,*)"Nx=", Nx, "Nz=", Nz, "h=", h, "tmax=", tmax
write(*,*)"KINDS Nx e Nz:", kind(Nx), kind(Nz)
write(*,*) "Modelo de velocidades:", file_vp_model
write(*,*) " Modelo de densidades:", file_density_model
write(*,*)"fx=", nfx, "fz=", nfz, "A=", A, "fcorte=", fcorte
write(*,*)"alfa=", alfa, "beta=", beta
write(*,*)"Fator cerjan=", fator_cerjan, "Tamanho da camada cerjan=",pontos_cerjan
write(*,*)"Number of shots=", n_shots
write(*,*) "Streamer head=", nstreamer_head
write(*,*) "Streamer tail=", nstreamer_tail
close(30)

end subroutine read_verify_input

!#################################################

subroutine read_verify_velocity_model()
implicit none

open(100,file=file_vp_model,status="unknown",access="direct",form="unformatted",recl=(Nx*Nz))
read(100,rec=1) ((mod_vel(i,j),i=1,Nz),j=1,Nx)
close(100)

open(60,file="model_verification.bin",status="replace",access="direct",form="unformatted",recl=(Nx*Nz))
write(60,rec=1) ((mod_vel(i,j),i=1,Nz),j=1,Nx)
close(60)

end subroutine read_verify_velocity_model

!#################################################

subroutine read_verify_density_model()
implicit none

do j=1, Nx
    do i=1, Nz
        if (nint(mod_vel(i,j)).EQ. 1524) then
            mod_den(i,j)= 1000.0
        else if (nint(mod_vel(i,j)).EQ. 1651) then
            mod_den(i,j)= 2000.0
        else if (nint(mod_vel(i,j)).EQ. 1943) then
            mod_den(i,j)= 2100.0
        else if (nint(mod_vel(i,j)).EQ. 2122) then
            mod_den(i,j)= 2200.0
        else if (nint(mod_vel(i,j)).EQ. 2492) then
            mod_den(i,j)= 2350.0
        else if (nint(mod_vel(i,j)).EQ. 2742) then
            mod_den(i,j)= 2500.0
        else if (nint(mod_vel(i,j)).EQ. 2918) then
            mod_den(i,j)= 2550.0
        else if (nint(mod_vel(i,j)).EQ. 2377) then
            mod_den(i,j)= 2400.0
        else if (nint(mod_vel(i,j)).EQ. 3179) then
            mod_den(i,j)= 2650.0
        else if (nint(mod_vel(i,j)).EQ. 2377) then
            mod_den(i,j)= 2400.0
        else if (nint(mod_vel(i,j)).EQ. 3407) then
            mod_den(i,j)= 2250.0
        else if (nint(mod_vel(i,j)).EQ. 3620) then
            mod_den(i,j)= 2700.0
        else if (nint(mod_vel(i,j)).EQ. 4511) then
            mod_den(i,j)= 2950.0
        end if
    end do
end do

open(99,file=file_density_model,status="replace",access="direct",form="unformatted",recl=(Nx*Nz))
write(99,rec=1) ((mod_den(i,j),i=1,Nz),j=1,Nx)
close(99)

open(61,file="density_model_verification.bin",status="replace",access="direct",form="unformatted",recl=(Nx*Nz))
write(61,rec=1) ((mod_den(i,j),i=1,Nz),j=1,Nx)
close(61)

end subroutine read_verify_density_model

!#################################################

subroutine allocate_variables()
implicit none

allocate (P(Nz,Nx))
allocate (U(Nz,Nx),V(Nz,Nx))
allocate (K_modulus(Nz,Nx),b(Nz,Nx))
allocate (mod_vel(Nz, Nx),mod_den(Nz,Nx))
allocate (vetor_cerjan(pontos_cerjan-1))
allocate (AMTT(Nz,Nx),MTT(Nz,Nx))
allocate (image(Nz,Nx))
end subroutine allocate_variables

!#################################################

subroutine calculate_source_modelling_parameters()
implicit none
print*, " "
print*, " "
write(*,*)"Parametros da fonte e de modelagem:"
fc=fcorte/(3*sqrt(pi))
write(*,*)"fc=", fc
t0=2*sqrt(pi)/fcorte
write(*,*)"t0=", t0
tm=2*t0
write(*,*)"tm=", tm
dt=h/(beta*4511.03) !######################################################ATENCAO VP MAX####################################
write(*,*) "dt=", dt
delta=dt/h 
write(*,*) "dt/h=", delta
U=0
V=0
P=0
K_modulus=0
b=0
passos=nint(tmax/dt)
passos_snap=passos
dt_snap=nint(passos_snap/nsnap)
write(*,*) "Numero de passos=", passos 
allocate (seismogram(passos,Nx))
seismogram=0
passos_fonte=nint(tm/dt)
write(*,*) "Fonte aplicada ate", passos_fonte, "passos"
allocate (sg_source(passos_fonte))
sg_source=0
do i=1,passos_fonte
    sg_source(i)=-A*(i*dt-t0)*exp(-pi*(pi*fc*(i*dt-t0))**2)
end do
AMTT=0
MTT=0
image=0
c_snap=0
c_image=0
vetor_cerjan=0


do j=1, Nx
    do i=1, Nz
        K_modulus(i,j)=(mod_vel(i,j)**2)*mod_den(i,j)
        b(i,j)=1/(mod_den(i,j))
    end do
end do

open(42,file="b_model_verification.bin",status="replace",access="direct",form="unformatted",recl=(Nx*Nz))
write(42,rec=1) ((b(i,j),i=1,Nz),j=1,Nx)
close(42)

end subroutine calculate_source_modelling_parameters

!#######################################################

real function staggered_grid_source(A,pi,fc,k,dt,t0) 
implicit none
integer :: k
real :: A, pi, fc, dt, t0

staggered_grid_source=-A*(k*dt-t0)*exp(-pi*(pi*fc*(k*dt-t0))**2)

end function staggered_grid_source

!#######################################################

subroutine acoustic_second_order_staggered_grid_operator_loop()
implicit none

do j=1,Nx-1
    do i=1,Nz-1
        U(i,j)=U(i,j)-b(i,j)*delta*(P(i,j+1)-P(i,j))
        V(i,j)=V(i,j)-b(i,j)*delta*(P(i+1,j)-P(i,j))   
    end do
end do

do j=2, Nx-1
    do i=2, Nz-1
        P(i,j)=P(i,j)-K_modulus(i,j)*delta*(U(i,j)-U(i,j-1)+V(i,j)-V(i-1,j))
    end do
end do

end subroutine acoustic_second_order_staggered_grid_operator_loop

!#######################################################

subroutine acoustic_fourth_order_staggered_grid_operator_loop()
implicit none

do j=2, Nx-2
    do i=2,Nz-2
        U(i,j)=U(i,j)-b(i,j)*(delta/24)*(-P(i,j+2)+27*P(i,j+1)-27*P(i,j)+P(i,j-1))
        V(i,j)=V(i,j)-b(i,j)*(delta/24)*(-P(i+2,j)+27*P(i+1,j)-27*P(i,j)+P(i-1,j))   
    end do
end do

do j=3, Nx-2
    do i=3, Nz-2
        P(i,j)=P(i,j)-K_modulus(i,j)*(delta/24) &
        & *(-U(i,j+1)+27*U(i,j)-27*U(i,j-1)+U(i,j-2)-V(i+1,j)+27*V(i,j)-27*V(i-1,j)+V(i-2,j))
        if(isnan(P(i,j))) then
            write(*,*) "Deu ruim"
         !   image(i,j)=0
        end if
    end do
end do

end subroutine acoustic_fourth_order_staggered_grid_operator_loop

!#######################################################

subroutine anti_reflection_cerjan_conditions()
implicit none

open (51,file="cerjan_check.txt")
do i=1,pontos_cerjan-1
    vetor_cerjan(i)=exp(-((fator_cerjan*(pontos_cerjan-1-i))**2))
    write(51,*) vetor_cerjan(i), i
end do
close(51)

!Applying the cerjan conditions on the left region
do j=1,pontos_cerjan-1
    do i=1,Nz
        P(i,j)=vetor_cerjan(j)*P(i,j)
    end do
end do


!Applying the cerjan conditions on the right region
n=pontos_cerjan
do j=Nx-pontos_cerjan+2, Nx 
    n=n-1
    do i=1,Nz
        P(i,j)=vetor_cerjan(j)*P(i,j)
    end do
end do


!Applying the cerjan conditions to the bottom region
n=pontos_cerjan
do i=Nz-pontos_cerjan+2,Nz
    n=n-1
    do j=1,Nx
        P(i,j)=vetor_cerjan(n)*P(i,j)
    end do
end do

end subroutine anti_reflection_cerjan_conditions

!#######################################################
subroutine transit_time_matrix()
implicit none

do j=1, Nx
    do i=1, Nz
        if (abs(P(i,j))>abs(AMTT(i,j))) then
            AMTT(i,j)=P(i,j)
            MTT(i,j)=k
        end if
    end do
end do

end subroutine transit_time_matrix

!#######################################################
subroutine RTM()
implicit none

do j=1,Nx
    do i=1,Nz
        if (k==MTT(i,j)) then
            image(i,j)=image(i,j)+(P(i,j))
            !image(i,j)=P(i,j)
        end if
        
        !if(isnan(image(i,j))) then
         !   write(*,*) "Deu ruim"
         !   image(i,j)=0
        !end if
    end do
end do

end subroutine RTM

!#######################################################
subroutine take_snapshot()
implicit none
open (20,file="snap"//trim(adjustl(snap_char))//".bin",status="replace",access="direct",form="unformatted", &
& recl=Nx*Nz)
open (27,file="snap"//trim(adjustl(snap_char))//"_U.bin",status="replace",access="direct",form="unformatted", &
& recl=Nx*Nz)
open (28,file="snap"//trim(adjustl(snap_char))//"_V.bin",status="replace",access="direct",form="unformatted", &
& recl=Nx*Nz)

write(20,rec=1) ((P(i,j),i=1,Nz),j=1,Nx)
write(27,rec=1) ((U(i,j),i=1,Nz),j=1,Nx)
write(28,rec=1) ((V(i,j),i=1,Nz),j=1,Nx)

close(20)
close(27)
close(28)


end subroutine take_snapshot

!#######################################################
subroutine save_image()
implicit none
open (37,file="image.bin",status="replace",access="direct",form="unformatted", &
& recl=Nx*Nz)

write(37,rec=1) ((image(i,j),i=1,Nz),j=1,Nx)

close(37)

!open(37,file="image.bin",status="unknown",access="direct",form="unformatted",recl=Nx*Nz)
!write(37,rec=1) image 
!close(37)

end subroutine save_image

!#######################################################
end program onda_2d
