program sanity
    use FFT
    use BSPHT
    implicit none
    integer :: i,j,k
    integer,parameter :: lat_division=16,bandwidth=lat_division/2,rad_division=2
    real*8,parameter :: density=1.0,delta_r=1.0
    complex(8),dimension(0:lat_division-1,0:lat_division-1,0:rad_division-1) :: shells
    
    !Create a uniform density grid
    call fill_the_mass(lat_division,rad_division,density,delta_r,shells)
    !Printing the mass Density
    print *,"Printing the mass Density"
    do i=0,rad_division-1
        do j=0,lat_division-1
            do k=0,lat_division-1
                print *,"r = ",(i+1.0/2.0),"m = ",j,"l = ",k,"coeff = ",shells(k,j,i)
            end do
        end do
    end do
    
    
    !Call spherical transform of all the shell
    do i=0,rad_division-1
        call SPH_transform(lat_division,bandwidth,shells(0:lat_division-1,0:lat_division-1,i))
    end do
    !Printing the Spherical Transform
    print *,"Printing the Spherical Transform, for each radius"
    do i=0,rad_division-1
        do j=0,lat_division-1
            do k=0,lat_division-1
                print *,"r = ",(i+1.0/2.0),"m = ",j,"l = ",k,"coeff = ",shells(k,j,i)
            end do
        end do
    end do
    
    
    
end program sanity

subroutine fill_the_mass(lat_division,rad_division,density,delta_r,shells)
    implicit none
    integer,intent(in) :: lat_division,rad_division
    real*8,intent(in) :: density,delta_r
    complex(8),dimension(0:lat_division-1,0:lat_division-1,0:rad_division-1),intent(out):: shells
    
    integer :: i,j,k
    real*8 :: mass_shell,mass_latitude,weight_norm
    real*8,parameter :: pi=4*atan(1.0)
    
    weight_norm=0.0
    do i=0,lat_division-1
        weight_norm=weight_norm+sin(pi*(2*i+1)/(4*(lat_division/2)))
    end do
    
    do i=0,rad_division-1
        mass_shell=density*4*pi*((i+1.0/2.0)*delta_r)**2
        print *,"Radius,Shell Mass : ",i,mass_shell
        do j=0,lat_division-1
            mass_latitude=mass_shell*sin(pi*(2*j+1)/(4*(lat_division/2)))/weight_norm
            do k=0,lat_division-1
                shells(j,k,i)=cmplx(mass_latitude/lat_division,0.0)
            end do
        end do
    end do
end subroutine fill_the_mass

subroutine copy_data_to_file(lat_division,bandwidth,rad_division,shells)
    implicit none
    integer,intent(in) :: lat_division,bandwidth,rad_division
    complex(8),dimension(0:lat_division-1,0:lat_division-1,0:rad_division-1),intent(in) :: shells
    
    integer :: i,j,k
    character*90,parameter :: write_file="spht_file.txt"
    
    !Now openeing the file to write the Spherical Transform data.
    open(unit=1,file=write_file,status='replace',form='formatted',action='write')
    10 format(f10.7,a1)
    11 format(f10.7)
    
    do i=0,lat_division-1
        if(i<=bandwidth)then
            
        else
            
        end if
    end do
    
end subroutine copy_data_to_file