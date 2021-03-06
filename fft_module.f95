module FFT
contains
!We have started from n=0 to N-1 in frequency also. So we are basically translating then
!left side (negetive frequency) to the right of the "positive" (See Numerical Receppe) freq.
!ie.(after 0-fn), due to prevailing periodic nature in frequency space 
!after -fn to fn (this block repeats.)

subroutine FFT3D(dim1,dim2,dim3,signal_array,transform_array)
    implicit none
    integer ::i,j,k,depth1,depth2,depth3
    integer,intent(in) :: dim1,dim2,dim3
    integer,allocatable,dimension(:,:) :: dump_credential
    complex(8),dimension(0:dim1-1,0:dim2-1,0:dim3-1),intent(in) :: signal_array
    complex(8),dimension(0:dim1-1,0:dim2-1,0:dim3-1),intent(inout) :: transform_array
    real*8,parameter :: pi=4*atan(1.0),tau=2*pi
    complex(8) :: omega1,omega2,omega3
    complex(8),dimension(0:dim3-1) :: temp_depth
    
    omega3=cmplx(cos(tau/dim3),-sin(tau/dim3))
    depth3=log2(real(dim3))
    
    !Transforming bread slices made of x-y plane and go to z direction.
    do k=0,dim3-1
        call FFT2D(dim1,dim2,signal_array(:,:,k),transform_array(:,:,k))
    end do
    
    !We have to transform the transform array. So again we 
    !need a single French Fries, to hold the cuboid.
    !Now we need to transform along each tunnel for value of n on a 
    !specific value of u,v now present in slice of bread for different n
    allocate(dump_credential(depth3,3))
    do i=0,dim1-1
        do j=0,dim2-1
            !temp_depth=transform_array(i,j,:)
            call FFT1D(dim3,depth3,1,0,omega3,dump_credential,&
                       transform_array(i,j,:),temp_depth)
            transform_array(i,j,:)=temp_depth(:)
        end do
    end do
    deallocate(dump_credential)
    
end subroutine FFT3D


subroutine FFT2D(dim1,dim2,signal_array,transform_array)
    implicit none
    integer :: i,j,depth2,depth1
    integer,intent(in) :: dim1,dim2
    complex(8),dimension(0:dim1-1,0:dim2-1),intent(in) :: signal_array
    complex(8),dimension(0:dim1-1,0:dim2-1),intent(inout) :: transform_array
    integer,allocatable,dimension(:,:) :: dump_credential
    real*8,parameter :: pi=4*atan(1.0),tau=2*pi
    complex(8),dimension(0:dim1-1) :: temp_col !could use it to make code look cleaner
    !but it could slow down and take memory according to dimension
    complex(8) :: omega1,omega2
    
    omega1=cmplx(cos(tau/dim1),-sin(tau/dim1))
    omega2=cmplx(cos(tau/dim2),-sin(tau/dim2))
    depth2=log2(real(dim2)) !when going row by row
    depth1=log2(real(dim1)) !when going col by col
    
    !Transforming : doing all the rows by row first.
    allocate(dump_credential(depth2,3))
    do i=0,dim1-1
        call FFT1D(dim2,depth2,1,0,omega2,dump_credential,&
                   signal_array(i,:),transform_array(i,:))
    end do
    deallocate(dump_credential)
    
    !Transforming : doing all the cols now. Here we hve to transform previous transformed 
    !               matrix.So we need new temp array.
    allocate(dump_credential(depth1,3))
    do j=0,dim2-1
        !temp_col(:)=transform_array(:,j) !BUG:We dont need to copy to temp.
        !Anyway.Its assigned at last level. So gets replaced. But still, no need here.
        call FFT1D(dim1,depth1,1,0,omega1,dump_credential,&
                   transform_array(:,j),temp_col)
        transform_array(:,j)=temp_col(:)
    end do
    deallocate(dump_credential)
    
end subroutine FFT2D

recursive subroutine FFT1D(total_length,total_depth,depth,start_index,omega,dump_credential,signal_array,transform_array)
    implicit none
    integer,intent(in) :: total_length,total_depth,depth,start_index
    complex(8),dimension(0:total_length-1),intent(in) :: signal_array
    integer :: stride_length
    integer,dimension(1:total_depth,3),intent(inout) :: dump_credential
    integer :: start1,end1,start2,end2,i,j,counter !for iterating and combining step
    complex(8),dimension(0:total_length-1),intent(inout) :: transform_array
    complex(8) :: temp_up,temp_down
    complex(8),intent(in) :: omega
    
    stride_length=2**depth    
    
    if(depth==total_depth+1)then
        !print *,"In last depth : ",dump_credential(depth-1,2) 
        !print *,"Start Index : ",start_index
        !We will need sepereate transform array
        !cuz there is only one to one memory mapping. But they will overlap.
        !ie. destroy values of susequent assignment. So we need extra array to store
        !transformed value. finally. We could de allocate it.
        ! Could be helpful for 2D or 3D implementation.
        transform_array(dump_credential(depth-1,2))=signal_array(start_index)*1 !Here 1 is FFT for 1 dimension.
    else
        !Calculating the first half of half-Fourier transform ,again will go recursive.
        !print *," ################# Start depth :",depth," #######################"
        !print *,"depth : ",depth
        dump_credential(depth,1)=1
        if(depth==1) then
            dump_credential(depth,2)=0
            start1=dump_credential(depth,2)
            dump_credential(depth,3)=(total_length/stride_length)-1
            end1=dump_credential(depth,3)
        else
            dump_credential(depth,2)=dump_credential(depth-1,2)
            start1=dump_credential(depth,2)
            dump_credential(depth,3)=dump_credential(depth-1,2)+(total_length/stride_length)-1
            end1=dump_credential(depth,3)
        end if
        !print *,"flag : ",dump_credential(depth,1),"starting : ",start1,&
                !"ending : ",end1
        call FFT1D(total_length,total_depth,depth+1,start_index,omega**2,&
                    dump_credential,signal_array,transform_array)
        
        !Calculating the secong half of the fourier transform.
        dump_credential(depth,1)=-1
        if(depth==1) then
            dump_credential(depth,2)=total_length/stride_length
            start2=dump_credential(depth,2)
            dump_credential(depth,3)=total_length-1
            end2=dump_credential(depth,3)
        else
            dump_credential(depth,2)=dump_credential(depth-1,2)+total_length/stride_length
            start2=dump_credential(depth,2)
            dump_credential(depth,3)=dump_credential(depth-1,3)
            end2=dump_credential(depth,3)
        end if
        !print *,"flag : ",dump_credential(depth,1),"starting : ",start2,&
                !"ending : ",end2
        call FFT1D(total_length,total_depth,depth+1,start_index+(stride_length/2),&
                    omega**2,dump_credential,signal_array,transform_array)
        
        !Combining the two half to backtrack our recursion.
        
        
        !do i=0,total_length-1
            !print *,transform_array(i)
        !end do
        
        counter=1
        !print *,"start1,end1,start2,end2 ",start1,end1,start2,end2
        !print *,"Omega**",omega
        
        do i=start1,end1
            j=end1+counter
            temp_up=transform_array(i)
            temp_down=transform_array(j)
            !Omega to be used is of previous level. eg.
            ! in last recombination ie FTn/2 and FTn/2 we use omega of FTn
            transform_array(i)=temp_up+(omega**(counter-1))*temp_down
            transform_array(j)=temp_up-(omega**(counter-1))*temp_down
            counter=counter+1
        end do
        
        !do i=0,total_length-1
            !print *,transform_array(i)
        !end do
        
        !print *,"###############  End Depth :",depth,"###########################"
    end if
    
end subroutine FFT1D

integer function log2(x)
    !WARNING : Check the output values of requires dimension before using whole code.
    !may cause ewrror if the actual value is less than theoretrical.
    implicit none
    real,intent(in)::x
    
    log2=int(log(x)/log(2.0))
    
end function log2

!This module calcuates the DCT transform but "naively".
!May implement using FFT.
!But naive one will be good if we want to integrate the 
!boundary specified wavenumber restriction on Bessel basis.
!For doing this transform it is specifically multip. by power of r
!based on the degree of the requirement.

!Correctly working individually, both real and complex part.
!Yeah it takes time to think.
subroutine customized_DCST1D(length,delta_r,r_power,signal,transform_matrix,transform)
    implicit none
    integer,intent(in) :: length,r_power
    real*8,intent(in) :: delta_r
    real*8,dimension(0:length-1,0:length-1),intent(inout) :: transform_matrix
    complex(8),dimension(0:length-1),intent(in) :: signal
    complex(8),dimension(0:length-1),intent(out) :: transform
    complex(8),dimension(0:length-1) :: temp0
    
    integer :: i,j
    real*8,dimension(0:length-1) :: power_array
    !Constricting the exter multiplication of r array.
    do i=0,length-1
        power_array(i)=((i+(1.0/2.0))*delta_r)**r_power
    end do
    temp0=power_array*signal
    !Now doing the transform
    transform=matmul(transform_matrix,power_array*signal)
    
    print *,""
    print *,"Inside DCST cant"
    do i=0,length-1
        do j=0,length-1
            !print *,i,j,transform_matrix(i,j)
        end do
    end do
    
    do i=0,length-1
        !print *,i,signal(i),temp0(i),transform(i)
    end do
    print *,"Leaving DCST cant"
end subroutine customized_DCST1D

!***K-USED
!We dont need to make new transform matrix each time. Basically its gonna be same each time
!between different k where we want to find them.
subroutine make_transform_matrix(type_flag,length,transform_matrix)
    !type_flag : which transform we want to do. (sin=0,cos=1)(Remember! :o)
    implicit none
    integer ::  i,j
    integer,intent(in) ::length,type_flag!Its gonnal be fixes for all l.(basically num_sphere)
    !This is equal sized transform. N->N element mapping.
    real*8,dimension(0:length-1,0:length-1),intent(out) :: transform_matrix
    real*8,parameter :: pi=4*atan(1.0),tau=2*pi
    
    !Constructing the transfrom Matrix
    do i=0,length-1
        do j=0,length-1
            if(type_flag==0)then
                transform_matrix(i,j)=sin(pi*(i+(1.0/2.0))*(j+(1.0/2.0))/&
                                                    (length-(1.0/2.0)))
            else if(type_flag==1)then
                transform_matrix(i,j)=cos(pi*(i+(1.0/2.0))*(j+(1.0/2.0))/&
                                                    (length-(1.0/2.0)))
            end if
        end do
    end do
    
end subroutine make_transform_matrix

end module FFT