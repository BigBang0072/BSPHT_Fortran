program sanity
    use FFT
    use BSPHT
    implicit none
    integer :: i,j,k
    integer,parameter :: lat_division=8,bandwidth=lat_division/2,rad_division=4
    real*8,parameter :: density=1.0,delta_r=1.0
    complex(8),dimension(0:lat_division-1,0:lat_division-1,0:rad_division-1) :: shells
    
    !Create a uniform density grid
    call fill_the_density(lat_division,rad_division,density,delta_r,shells)
    !Printing the mass Density
    print *,"Printing the mass Density"
    do i=0,rad_division-1
        do j=0,lat_division-1
            do k=0,lat_division-1
                print *,"r = ",(i+1.0/2.0),"phi = ",j,"theta = ",k,"coeff = ",shells(k,j,i)
            end do
        end do
    end do
    
    
    !CALL spherical transform of all the shell
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
    
    !CALL
    call naive_SpBessel_transform(bandwidth,lat_division,rad_division,delta_r,shells)
    
    !Writing the data to file(transform data)
    call copy_data_to_file(lat_division,bandwidth,delta_r,rad_division,shells)
    
    print *,"Printing the Bessel-Spherical Transform, for each k,l,m"
    do i=0,rad_division-1
        do j=0,lat_division-1
            do k=0,lat_division-1
                print *,"f = ",(i+1.0/2.0)/(rad_division-1.0/2.0)/(2*delta_r),"m = ",j,"l = ",k,"coeff = ",shells(k,j,i)
            end do
        end do
    end do
    
end program sanity

!We dont have to fill teh mass, we have to fill the density.So it has to be divided in actual data. CHANGE
subroutine fill_the_density(lat_division,rad_division,density,delta_r,shells)
    implicit none
    integer,intent(in) :: lat_division,rad_division
    real*8,intent(in) :: density,delta_r
    complex(8),dimension(0:lat_division-1,0:lat_division-1,0:rad_division-1),intent(out):: shells
    
    integer :: i,j,k
    real*8 :: mass_shell,mass_latitude,weight_norm
    real*8,parameter :: pi=4*atan(1.0)
    
    weight_norm=0.0
    do i=0,lat_division-1
        !*********************IMPORTANT**************
        !We have to fill density.
        !In real problem we would have to convert mass to density then we have to divide what we multiplied to density here.
        
        weight_norm=weight_norm+sin(pi*(2*i+1)/(4*(lat_division/2)))
        
    end do
    
    do i=0,rad_division-1
        mass_shell=density*4*pi*((i+1.0/2.0)*delta_r)**2
        print *,"Radius,Shell Mass : ",i,mass_shell
        do j=0,lat_division-1
            mass_latitude=mass_shell*sin(pi*(2*j+1)/(4*(lat_division/2)))/weight_norm
            do k=0,lat_division-1
                !shells(j,k,i)=cmplx(mass_latitude/lat_division,0.0)
                shells(j,k,i)=cmplx(density,0.0)     !We have to represrnt density in terms of Eigen function.
            end do
        end do
    end do
end subroutine fill_the_density

subroutine copy_data_to_file(lat_division,bandwidth,delta_r,rad_division,shells)
    implicit none
    integer,intent(in) :: lat_division,bandwidth,rad_division
    real*8,intent(in) :: delta_r
    complex(8),dimension(0:lat_division-1,0:lat_division-1,0:rad_division-1),intent(in) :: shells
    
    integer :: i,j,k,start_index,m_val
    real*8,dimension(0:rad_division-1) :: k_values
    real*8,parameter :: pi=4*atan(1.0)
    character*90,parameter :: write_file="spht_file.txt"
    
    !Generating the k-values
    do i=0,rad_division-1
        k_values(i)=(2*pi)*(i+1.0/2.0)/(rad_division-1.0/2.0)/(2*delta_r)
    end do
    
    !Now openeing the file to write the Spherical Transform data.
    open(unit=1,file=write_file,status='replace',form='formatted',action='write')
    8 format(a100)
    9 format(i3,f10.7)
    10 format(i3,i3,a1)
    11 format(f10.7,a1,f10.7,a1)
    12 format(f10.7,a1,f10.7)
    
    do i=0,lat_division-1
        if(i<bandwidth)then
            start_index=i
            m_val=i
        else if(i>bandwidth)then
            start_index=bandwidth-(i-bandwidth)
            m_val=-(start_index)
        else
            cycle
        end if
        do j=start_index,bandwidth-1
            do k=0,rad_division-1
                !Writing the klm sequence.(header first)
                if(k<rad_division-1)then
                    if(k==0)then
                        write(unit=1,fmt=10,advance='no'),j,m_val,":"
                    end if
                    write(unit=1,fmt=11,advance='no'),real(shells(j,i,k)),",",imag(shells(j,i,k)),";"
                else
                    write(unit=1,fmt=12,advance='yes'),real(shells(j,i,k)),",",imag(shells(j,i,k))
                end if
            end do
        end do
    end do
    
    write(unit=1,fmt=8,advance='yes'),"Writing the available k-values under Nyquist Limit, not orthogonal."
    do i=0,rad_division-1
        write(unit=1,fmt=9,advance='yes'),i,k_values(i)
    end do
    
    end file 1
    rewind 1
    close(unit=1)
    
end subroutine copy_data_to_file