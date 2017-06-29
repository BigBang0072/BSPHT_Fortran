program test_legendre
    use BSPHT
    implicit none
    integer :: i
    integer*8 :: j=9
    !lat: is latitude(direction of theta), long: longitude(direction of phi)
    integer,parameter :: lat_length=32,bandwidth=lat_length/2 
    
    real*8,parameter :: pi=4*atan(1.0),tau=2*pi
    real*8,dimension(0:lat_length-1) :: theta_j,res
    real*8,dimension(1:lat_length-2) :: check
    
    complex,dimension(0:lat_length-1) :: signal
    
    do i=0,lat_length-1
        theta_j(i)=pi*(2*i+1)/(4*bandwidth)
        res(i)=0.0
    end do
    
    do i=1,lat_length-2
        check(i)=14
    end do
    
    !call aLegendre_lm_vector(-9,12,lat_length,cos(theta_j),res)
    call make_quadrature_weight(lat_length,res)
    print *,"########################"
    do i=0,lat_length-1
        print *,i,res(i)
    end do
    
    res(0:lat_length-1)=0
    print *,"########################"
    do i=0,lat_length-1
        print *,i,res(i)
    end do
    
    res(1:lat_length-2)=check(1:lat_length-2)
    print *,"########################"
    do i=0,lat_length-1
        print *,i,res(i)
    end do
    
end program test_legendre