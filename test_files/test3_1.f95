program test3
    use BSPHT
    use FFT
    implicit none
    integer :: i,j
    integer,parameter :: length=4,l=2
    real*8,parameter :: delta_r=1.0
    real*8,dimension(0:length-1) :: k_values,r_values
    complex(8),dimension(0:length-1) :: signal,transform
    
    print *,"Printing k_values"
    do i=0,length-1
        k_values(i)=(i+1.0/2.0)/(length-1.0/2.0)/(2*delta_r)
        r_values(i)=(i+1.0/2.0)*delta_r
        print *,i,r_values(i),k_values(i)
    end do
    print *,"############"
    
    signal(0:length-1)=(/cmplx(2,0.0),&
                         cmplx(5,0.0),&
                         cmplx(-13,0.0),&
                         cmplx(-18,0.0)/)
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