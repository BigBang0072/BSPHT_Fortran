program bessel_test
    use FFT
    use BSPHT
    implicit none
    integer,parameter :: length=8
    integer :: i,j,k
    real*8,dimension(0:length-1,0:length-1) :: transform_matrix
    real*8,parameter :: pi=4*atan(1.0)
    real*8,parameter :: delta_r=1
    real*8,dimension(0:length-1) :: k_values,r_values
    complex(8),dimension(0:length-1) :: signal,transform
    complex(8),dimension(0:1,0:length-1,0:2) :: dcst_transforms
    complex(8),dimension(0:length-1,0:3,0:1) :: interpolation_tables
    
    !1. Testing teh doule factorial and factorial function.
    !print *,double_factorial(4)
    !print *,factorial(4)
    
    print *,""
    print *,""
    print *,""
    print *,"Starting the PHASE 2."
    !2. Test: make transform matrix (sin:0, cos: 1 here)
    call make_transform_matrix(1,length,transform_matrix)
        
    !do j=0,length-1
        !i=4
        !signal(j)=cmplx(cos(pi*(i+1.0/2.0)*(j+1.0/2.0)/(length-1.0/2.0)),&
                        !sin(pi*(i+1.0/2.0)*(j+1.0/2.0)/(length-1.0/2.0)))
    !end do
    !making the signal
    signal(0:length-1)=(/cmplx(0.0004627333629297397,0.0),&
                        cmplx(0.00414809773935782,0.0),&
                        cmplx(0.01143123671886824,0.0),&
                        cmplx(0.0221389940687361,0.0),&
                        cmplx(0.03601664614108058,0.0),&
                        cmplx(0.0527337826212236,0.0),&
                        cmplx(0.07189193076641275,0.0),&
                        cmplx(0.0930337440467297,0.0)/)
    
    print *,"Printing Signal!!!"
    do j=0,length-1
        print *,signal(j)
        
    end do
    print *,"########"
    do i=0,length-1
        do j=0,length-1
            !print *,transform_matrix(i,j)
        end do
        !print *,"*******"
    end do
    !transform=matmul(transform_matrix,signal)
    !call customized_DCST1D(length,delta_r,0,signal,transform_matrix,transform)
    do i=0,length-1
        !print *,transform(i)
    end do
    
    print *,""
    print *,""
    print *,""
    print *,"Starting Phase 3 TESTING!!!"
    !3.generate DCST_transforms test.
    
    call generate_DCST_transforms(length,delta_r,signal,dcst_transforms)
    do i=0,1
        do j=0,2
            do k=0,length-1
                print *,i,j,k,dcst_transforms(i,k,j)
            end do
            print *,"******"
        end do
        print *,"/////////"
    end do
    
    
    print *,""
    print *,""
    print *,""
    print *,"Printing Phase 4 TEsting"
    
    do i =0,length-1
        k_values(i)=(i+1.0/2.0)/(length-1.0/2.0)/(2*delta_r)
    end do
    print *,"Printing the k_values available in this transform. acc. to our delta_r"
    do i=0,length-1
        print *,k_values(i)
    end do
    print *,"########### Printing r_values fro testing the transform."
    do i=0,length-1
        r_values(i)=((i+1.0/2.0)*delta_r)
        print *,r_values(i)
    end do
    do i=2,7
        print *,i,sum(signal*(r_values**i))
    end do
    print *,"###########"
    call segment_coefficient(length,delta_r,k_values,signal,dcst_transforms,interpolation_tables)
    print *,""
    do i=0,1
        do j=0,length-1
            do k=0,3
                print *,i,j,k,interpolation_tables(j,k,i)
            end do
            print *,"****"
        end do     
        print *,"***********"
    end do
    
    
    
    print *,""
    print *,""
    print *,""
    print *,"Testing Phase 5 ::: "
    
    call SpBessel_transform(2,length,delta_r,signal,transform)
    
    do i=0,length-1
        print *,k_values(i),transform(i)
    end do
    
    
end program bessel_test