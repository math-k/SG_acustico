program stack

!##################################################

!Module containing all necessary program variables
use module_variables

implicit none

final_image=0
call read_verify_input()
call stack_images()

contains
!#################################################

subroutine stack_images()
implicit none

do i=1,5
    write(image_char,"(i3)"), i
    open (37,file="image"//trim(adjustl(image_char))//".bin",status="replace",access="direct",form="unformatted", &
    & recl=Nx*Nz)

    read(37,rec=1) image
    final_image=final_image+image
    close(37)
end do 

open (39,file="image"//trim(adjustl(image_char))//".bin",status="replace",access="direct",form="unformatted", &
& recl=Nx*Nz)
write(39,rec=1) final_image
close(39)

end subroutine stack_images
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
Nx=nint(x/h)
Nz=nint(z/h)
nfx=nint(fx/h)+1
nfz=nint(fz/h)+1
write(*,*)"Nx=", Nx, "Nz=", Nz, "h=", h, "tmax=", tmax
write(*,*)"KINDS Nx e Nz:", kind(Nx), kind(Nz)
write(*,*) "Modelo de velocidades:", file_vp_model
write(*,*) " Modelo de densidades:", file_density_model
write(*,*)"fx=", nfx, "fz=", nfz, "A=", A, "fcorte=", fcorte
write(*,*)"alfa=", alfa, "beta=", beta
write(*,*)"Fator cerjan=", fator_cerjan, "Tamanho da camada cerjan=",pontos_cerjan
close(30)

end subroutine read_verify_input
!#######################################################
end program stack