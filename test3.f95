program test3
    use BSPHT
    use FFT
    implicit none
    integer :: i,j
    integer,parameter :: length=4,l=1
    real*8,parameter :: pi=4*atan(1.0)
    real*8,parameter :: delta_r=1.0
    real*8,dimension(0:length-1) :: k_values,r_values
    complex(8),dimension(0:length-1) :: signal,transform
    
    print *,"Printing k_values"
    do i=0,length-1
        k_values(i)=(2*pi)*(i+1.0/2.0)/(length-1.0/2.0)/(2*delta_r)
        r_values(i)=(i+1.0/2.0)*delta_r
        print *,i,r_values(i),k_values(i)
    end do
    print *,"############"
    
    signal(0:length-1)=(/cmplx(0.21439301,0.0),&
                         cmplx(0.43572954,0.0),&
                         cmplx(0.27000044,0.0),&
                         cmplx(-0.04503164,0.0)/)
    print *,"Printing Signal"
    do i=0,length-1
        print *,signal(i)
    end do
    print *,"############"
    
    call SpBessel_transform(l,length,delta_r,signal,transform)
    
    print *,"Printing Transform"
    do i=0,length-1
        print *,transform(i)
    end do
end program test3