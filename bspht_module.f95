module BSPHT
contains
!This will be used to complete teh spherical Harmnics of particular layer 
!ie. Spherical Layer.
subroutine SPH_transform(lat_length,bWidth,signal_r)
    implicit none
    integer,intent(in) :: lat_length,bWidth
    complex,dimension(0:lat_length-1,0:lat_length-1),intent(inout) :: signal_r 
    
    integer :: i,j
    real*8,parameter :: pi=4*atan(1.0),tau=2*pi
    real*8,dimension(0:lat_length-1) :: theta_j
    
    !This will do the FFT of all the rows. First Step to Successs.
    call do_FFT(lat_length,bWidth,signal_r)
    print *,"########"
    do i=0,lat_length-1
        do j=0,lat_length-1
            print *,signal_r(i,j)
        end do
        print *,"******"
    end do
    print *,"########"
    !now start the Shree Ganesh of our "Legendre" Transform.
    call make_theta_j(lat_length,theta_j)
    !do i=0,lat_length-1
        !print *,theta_j(i)
    !end do
    print *,"HIK"
    call naive_legendre_transform(lat_length,bWidth,cos(theta_j),signal_r)
    print *,"KIH"
    
end subroutine SPH_transform

!THis will be used in above to do the FFT of signal array.
!Its for more abstraction.
subroutine do_FFT(lat_length,bWidth,signal_r)
    use FFT
    implicit none
    integer,intent(in) :: lat_length,bWidth
    complex,dimension(0:lat_length-1,0:lat_length-1),intent(inout):: signal_r
    
    integer :: i,depth
    integer,allocatable,dimension(:,:) :: dump_credential
    real*8,parameter :: pi=4*atan(1.0),tau=2*pi
    complex :: omega
    complex,dimension(0:lat_length-1) :: temp_row
    
    omega=cmplx(cos(tau/lat_length),-sin(tau/lat_length)) !its 2B,So.
    depth=log2(real(lat_length))
    allocate(dump_credential(depth,3))
    
    do i=0,lat_length-1
        call FFT1D(lat_length,depth,1,0,omega,dump_credential,&
                    signal_r(i,0:lat_length-1),temp_row)
        signal_r(i,0:lat_length-1)=temp_row
    end do
    deallocate(dump_credential)
    
end subroutine do_FFT

!This function calculates the Legendre Transform "NAIVELY"
subroutine naive_legendre_transform(lat_length,bWidth,theta_j,signal)
    implicit none
    integer,intent(in) :: lat_length,bWidth
    real*8,dimension(0:lat_length-1),intent(in) :: theta_j
    !you know actually the m goes from -b to b not 0 to 2b-1.(down dim)
    complex,dimension(0:lat_length-1,0:lat_length-1),intent(inout) :: signal
    
    integer :: i
    real*8,allocatable,dimension(:,:) :: transform_matrix
    real*8,dimension(0:lat_length-1) :: weight_vector
    complex,allocatable,dimension(:) :: transform
    
    
    !m:bWidth-1,0:lat_length-1
    !m:bWidth-1
    call make_quadrature_weight(lat_length,weight_vector)
    !Traversing along the m which is already fourier transformed.
    do i=0,bWidth
        print *,i,(lat_length-1)-i+1,bWidth,lat_length
        allocate(transform_matrix(i:bWidth-1,0:lat_length-1))
        !print *,"Prep1"
        allocate(transform(i:bWidth-1))
        !print *,"Prep2"
        call aLegendre_transform_matrix(i,bWidth,lat_length,theta_j,transform_matrix)
        !print *,"Prepped"
        if(i==0)then
            transform=matmul(transform_matrix,&
                        weight_vector*signal(0:lat_length-1,i))
            !print *,"*******"
            !print *,transform
            !print *,"*******"
            signal(0:lat_length-1,i)=cmplx(0.0,0.0)
            !signal(i:bWidth-1,i)=transform(i:bWidth-1)
            call copy_the_transform(i,i,lat_length,bWidth,signal,transform)
        else if(i==bWidth)then
            signal(0:lat_length-1,i)=cmplx(0.0,0.0)
        else
            
            transform=matmul(transform_matrix,&
                        weight_vector*signal(0:lat_length-1,i))
            !print *,"1K"
            signal(0:lat_length-1,i)=cmplx(0.0,0.0)
            !print *,"2K"
            !print *,"******"
            !print *,transform
            !print *,"******"
            !signal(i:bWidth-1,i)=transform(i:bWidth-1)
            call copy_the_transform(i,i,lat_length,bWidth,signal,transform)
            !print *,"3K"
            transform=matmul(transform_matrix,&
                        weight_vector*signal(0:lat_length-1,(lat_length-1)-i+1))
            !print *,"4K"
            transform=transform*(-1)**(i)
            !print *,"5K"
            signal(0:lat_length-1,(lat_length-1)-i+1)=cmplx(0.0,0.0)
            !print *,"******"
            !print *,transform
            !print *,"******"
            !signal(i:bWidth-1,(lat_length-1)-i+1)=transform(i:bWidth-1)
            call copy_the_transform(i,(lat_length-1)-i+1,lat_length,bWidth,signal,transform)
        end if
        deallocate(transform)
        deallocate(transform_matrix)
        !print *,"end",i
    end do
    !For doing the subsampling and smoothing we have to subsample the
    !transform matrix and our weight*signal vector, then multiply 
    !with lesser number of sample.
    
end subroutine naive_legendre_transform

!thIS MODULE IS JUST TO COPY ELEMENTS FROM THE TRANSFORM TO signal
!array. Some probelm is there in direct copyting.
subroutine copy_the_transform(i,m,lat_length,bWidth,signal,transform)
    implicit none
    integer,intent(in) :: i,m,lat_length,bWidth
    complex,dimension(0:lat_length-1,0:lat_length-1),intent(inout) :: signal
    complex,dimension(i:bWidth-1),intent(in) :: transform
    
    integer :: k
    
    do k=i,bWidth-1
        signal(k,m)=transform(k)
    end do
    
end subroutine copy_the_transform

!This function create teh weight for integral approximation
!using Chebyshev Quadrature.(Formula copied)
subroutine make_quadrature_weight(length,weight_vector)
    implicit none
    integer,intent(in) :: length
    integer :: i,j,bWidth
    real*8,dimension(0:length-1),intent(inout) :: weight_vector
    real*8,parameter :: pi=4*atan(1.0)
    real*8 :: temp_sum
    
    bWidth=length/2
    do i=0,length-1
        temp_sum=0
        do j=0,bWidth-1
            temp_sum=temp_sum+1/(2*j+1)*sin((2*i+1)*(2*j+1)*pi/(4*bWidth))
        end do
        weight_vector(i)=temp_sum*(2.0/bWidth)*sin(pi*(2*i+1)/(4*bWidth))
        !print *,(2.0/bWidth),temp_sum,weight_vector(i)
    end do  
end subroutine make_quadrature_weight

!Beware the vector returned and take in for transform are real*8.
!They shoot to large value very fast, so using real*8.
!Make sure you convert to necessary format before multiplying.
subroutine aLegendre_lm_vector(m,l,v_len,arg,aLegendre_vector)
    implicit none
    !l is always >=0.
    integer,intent(in) :: l,m,v_len
    real*8,dimension(0:v_len-1),intent(in) :: arg
    real*8,dimension(0:v_len-1),intent(out) :: aLegendre_vector
    
    integer :: i,ex_recurs
    real*8,dimension(0:v_len-1) :: Pm_m,Pm_m1,Pm_temp
    
    !print *,l,m,v_len
    !Directly using the closed form of Legendre Polynomial to start
    !recurrence formula for getting higher order
    Pm_m=((-1)**abs(m))*(double_factorial(2*abs(m)-1))*(1-arg**2)**(abs(m*1.0)/2)
    !print *,m,l,double_factorial(2*abs(m)-1)
    !print *,double_factorial(5)
    !do i=1,v_len-1
        !print *,arg(i),Pm_m(i)
    !end do
    
    Pm_m1=arg*Pm_m*(2*abs(m)+1)
    
    if(l==abs(m))then
        if(m>=0)then
            aLegendre_vector=Pm_m
        else
            aLegendre_vector=((-1)**abs(m))*(Pm_m)
        end if
    else if(l==abs(m)+1)then
        if(m>=0)then
            aLegendre_vector=Pm_m1
        else
            aLegendre_vector=((-1)**abs(m))*(Pm_m1)
        end if
    else if(l>m)then
        !print *,m
        ! To caluclate all the elemnts of transform array in one go
        ! for a particluar m, just keep addig the p_temp to new row of 
        ! transform matrix. It will be calculated in one go.
        ex_recurs=l-abs(m)-1
        !print *,ex_recurs
        do i=1,ex_recurs
            Pm_temp=(arg*(2*(abs(m)+i+1)-1)*Pm_m1-((abs(m)+i+1)+abs(m)-1)*Pm_m)/((abs(m)+i+1)-abs(m))
            Pm_m=Pm_m1
            Pm_m1=Pm_temp
        end do
            
        if(m<0)then
            aLegendre_vector=((-1)**abs(m))*Pm_temp
        else
            aLegendre_vector=Pm_temp
        end if
    else if (l<m)then
        aLegendre_vector=arg*0
    end if
end subroutine aLegendre_lm_vector

!*******   *******   ******   *******   ******   *******   *******   *******
!This subroutine just give the transform matrix in one shot rahter than calling 
!the above routine again and again.

!The transform for m and -m will have similar co-related transform matrix.
!So they could be done in one shot.
subroutine aLegendre_transform_matrix(m,bWidth,lat_length,arg,transform_matrix)
    implicit none
    !l is always >=0.
    integer,intent(in) :: m,bWidth,lat_length
    real*8,dimension(0:lat_length-1),intent(in) :: arg
    real*8,dimension(m:bWidth-1,0:lat_length-1),intent(out) :: transform_matrix
    
    integer :: i,ex_recurs
    real*8,dimension(0:lat_length-1) :: Pm_m,Pm_m1,Pm_temp
    
    !Directly using the closed form of Legendre Polynomial to start
    !recurrence formula for getting higher order
    !No need for abs cuz it will be called for just for positive make_quadrature_weight
    !For negetive m it will be change there only.
    Pm_m=((-1)**abs(m))*(double_factorial(2*abs(m)-1))*(1-arg**2)**(abs(m*1.0)/2)
    Pm_m1=arg*Pm_m*(2*abs(m)+1)
    
     
    !print *,m,bWidth,double_factorial(2*abs(m)-1)
    !do i=0,lat_length-1
        !print *,arg(i),Pm_m(i)
    !end do
    if(m==bWidth)then
        return
    end if
    transform_matrix(abs(m),0:lat_length-1)=Pm_m(0:lat_length-1)
    if(m==bWidth-1)then
        return
    end if
    transform_matrix(abs(m)+1,0:lat_length-1)=Pm_m1(0:lat_length-1)
    !highest bandwidth ie l=bWidth-1
    ex_recurs=(bWidth-1)-abs(m)-1
    !print *,ex_recurs
    do i=1,ex_recurs
        Pm_temp=(arg*(2*(abs(m)+i+1)-1)*Pm_m1-((abs(m)+i+1)+abs(m)-1)*Pm_m)/((abs(m)+i+1)-abs(m))
        transform_matrix(abs(m)+1+i,0:lat_length-1)=Pm_temp(0:lat_length-1)
        Pm_m=Pm_m1
        Pm_m1=Pm_temp
    end do
    
end subroutine aLegendre_transform_matrix

!This is to create the argument for trasform matrix and all
!Just to prettify the final subroutine,make it look modular.
subroutine make_theta_j(lat_length,theta_j)
    implicit none
    integer :: i
    integer,intent(in) :: lat_length
    real*8,dimension(0:lat_length-1),intent(out) :: theta_j
    real*8,parameter :: pi=4*atan(1.0)
    
    do i=0,lat_length-1
        theta_j(i)=pi*(2*i+1)/(4*(lat_length/2))
    end do
    
end subroutine make_theta_j

!This is a subsidary function used in inner calculation of 
!calculating Legendre polynomial.
integer*8 function double_factorial(num)
    implicit none
    integer :: i
    integer*8 :: temp
    integer,intent(in) :: num
    
    temp=1
    !  DEFINITION :   (-1)!! =1,  (0)!! =1
    do i=1,num
        if(mod(i,2)/=0)then
            temp=temp*i
        end if
    end do
    double_factorial=temp
end function double_factorial

end module BSPHT