program test
    use FFT
    use BSPHT
    implicit none
    integer :: i,j,nm
    integer,parameter :: l=0,rad_length=4,spherical_band=4
    real*8,parameter :: delta_r=1.0,pi=4*atan(1.0)
    real*8,dimension(0:rad_length-1) :: r_val,k_val
    real*8,dimension(0:rad_length-1,0:rad_length-1,0:spherical_band-1) :: transfrom_mat
    real*8,dimension(0:l) :: sj,dj
    real*8 :: x
    
    complex(8),dimension(0:rad_length-1) :: signal,transform
    
    do i=0,rad_length-1
        r_val(i)=(i+1.0/2.0)*delta_r
        k_val(i)=(2*pi)*(i+1.0/2.0)/(rad_length-1.0/2.0)/(2*delta_r)
        
    end do
    
    signal(0)=cmplx(0.99162858,0.0)
    signal(1)=cmplx(0.92616052,0.0)
    signal(2)=cmplx(0.80300443,0.0)
    signal(3)=cmplx(0.63661977,0.0)
    !x=3.14159265*0.5
    !call sphj(l,x,nm,sj,dj)
    !print *,sj(1)
    
    call create_transfrom_matrix(rad_length,spherical_band,r_val,k_val,transfrom_mat)
    print *,"############"
    do i=0,rad_length-1
        do j=0,rad_length-1
            print *,transfrom_mat(i,j,0)
        end do
        print *,"***"
    end do
end program test