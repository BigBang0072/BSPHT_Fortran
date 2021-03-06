program test_legendre
    use BSPHT
    use FFT
    implicit none
    integer,parameter :: m=5
    integer :: i,j
    !integer*8 :: j=9
    !lat: is latitude(direction of theta), long: longitude(direction of phi)
    integer,parameter :: lat_length=8,bandwidth=lat_length/2 
    
    real*8,parameter :: pi=4*atan(1.0),tau=2*pi
    real*8,dimension(0:lat_length-1) :: theta_j,sig,weigh
    
    !Remember the dimension starts from m. So change.
    real*8,dimension(m:bandwidth-1,0:lat_length-1) :: trans
    real*8,dimension(m:bandwidth-1) :: ans
    complex,dimension(0:lat_length-1) :: signal
    complex,dimension(0:lat_length-1,0:lat_length-1) :: signal_fft,transform_fft
    
    do i=0,lat_length-1
        theta_j(i)=pi*(2*i+1)/(4*bandwidth)
        weigh(i)=0.0
        sig(i)=0.0
    end do
    
    do i=0,lat_length-1
        !print *,res(i)
    end do
    
    !call make_quadrature_weight(lat_length,weigh)
    !call aLegendre_lm_vector(m,15,lat_length,cos(theta_j),sig)
    !call aLegendre_transform_matrix(m,bandwidth,lat_length,cos(theta_j),trans)
        
    
    
    do i=0,lat_length-1
        do j=0,lat_length-1
            signal_fft(i,j)=cmplx(cos(j*tau/lat_length),0.0)
            !transform_fft(i,j)=cmplx(0.0,0.0)
        end do
    end do
    do i=0,lat_length-1
        do j=0,lat_length-1
            print *,i,j,signal_fft(i,j)
            !transform_fft(i,j)=cmplx(0.0,0.0)
        end do
        print *,"********"
    end do
    
    
    call do_FFT(lat_length,bandwidth,signal_fft)
    
    print *,"####################"
    do i=0,lat_length-1
        do j=0,lat_length-1
            print *,signal_fft(i,j)
        end do
        print *,"***********"
    end do
    print *,"#######################"
    !print *,"########################"
    !do i=0,lat_length-1
        !print *,i,sig(i)
    !end do
    
    !print *,"########################"
    !do i=m,bandwidth-1
        !do j=0,lat_length-1
            !print *,i,j,theta_j(j),trans(i,j)
        !end do
        !print *,"**********************"
    !end do
    
    !ans=matmul(trans,weigh*sig)
    
    !do i=m,bandwidth-1
        !print *,i,ans(i)
    !end do
  
    
end program test_legendre