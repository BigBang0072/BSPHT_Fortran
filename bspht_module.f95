module BSPHT
contains

!X For all the ls in the signal we ahve to do this.
!X For doing complete transfrom call from main program. for each l.
! We can do for all l at teh same place cuz we can create 
!transfrom matrix for all l in one go.
subroutine naive_SpBessel_transform(bandwidth,lat_length,length,delta_r,signal_array)
    implicit none
    !The length here is just the radial length. I am lazy :)
    integer,intent(in) :: lat_length,length,bandwidth!bandwidth is from SPH
    real*8,intent(in) :: delta_r
    complex(8),dimension(0:lat_length-1,0:lat_length-1,0:length-1),intent(inout) :: signal_array
    
    !Workspace Variable
    integer :: i,j,start_index
    real*8,parameter :: pi=4*atan(1.0)
    real*8,dimension(0:length-1) :: f_values,k_values,r_values
    real*8,dimension(0:length-1,0:length-1,0:bandwidth-1) :: super_transform_matrix
    complex(8),dimension(0:length-1) :: temp_transform
    
    !We need these in making of transform matrix.
    do i=0,length-1
        r_values(i)=(i+1.0/2.0)*(delta_r)
        f_values(i)=(i+1.0/2.0)/(length-1.0/2.0)/(2*delta_r)
        k_values(i)=(2*pi)*(i+1.0/2.0)/(length-1.0/2.0)/(2*delta_r)
    end do
    
    call create_transfrom_matrix(length,bandwidth,r_values,k_values,super_transform_matrix)
    !do i=0,length-1
        !do j=0,length-1
            !print *,super_transform_matrix(i,j,0)
        !end do
        !print *,"*******"
    !end do
    
    !Doing transform l by l for each m then puttin tehm in the same place.
    do i=0,lat_length-1 !m ka loop
        if(i<bandwidth)then
            start_index=i
        else if(i>bandwidth)then
            start_index=bandwidth-(i-bandwidth)
        else
            cycle
        end if
        do j=start_index,bandwidth-1
            temp_transform=matmul(super_transform_matrix(0:length-1,0:length-1,j),&
                            signal_array(j,i,0:length-1))
            signal_array(j,i,0:length-1)=temp_transform(0:length-1)
        end do
    end do
    
end subroutine naive_SpBessel_transform

subroutine create_transfrom_matrix(length,bandwidth,r_values,k_values,transform_matrix)
    implicit none
    !BEWARAE the radial abndwidth is same as length.
    integer,intent(in) :: bandwidth,length
    real*8,dimension(0:length-1),intent(in) :: r_values,k_values
    real*8,dimension(0:length-1,0:length-1,0:bandwidth-1),intent(out) :: transform_matrix
    
    integer :: i,j,dummy
    real*8,dimension(0:bandwidth-1) :: temp_arr,dummy_arr
    
    do i=0,length-1
        do j=0,length-1
            call sphj(bandwidth-1,k_values(i)*r_values(j),dummy,temp_arr,dummy_arr)
            !print *,temp_arr
            transform_matrix(i,j,0:bandwidth-1)=temp_arr
        end do
    end do
    
end subroutine create_transfrom_matrix

! @ PARTICULAR "l"
!This subroutine contains the Bessel transform to be done after 
!Spherical Transform for a particular m,l.
subroutine SpBessel_transform(l,length,delta_r,signal_array,transform_array)
    use FFT
    implicit none
    integer,intent(in) :: l,length
    real*8,intent(in) :: delta_r
    complex(8),dimension(0:length-1),intent(in) ::  signal_array
    complex(8),dimension(0:length-1),intent(out) :: transform_array
    
    integer :: i,j,k,small_n,small_n_max
    real*8,parameter :: pi=4*atan(1.0)
    real*8,dimension(0:length-1) :: k_values
    complex(8),allocatable,dimension(:,:,:) :: interpolation_table
    complex(8),dimension(0:1,0:length-1,0:2) :: dcst_transforms
    complex(8),allocatable,dimension(:,:,:) :: integration_table
    
    small_n=l/2
    print *,"small_n : ",small_n
    !small_n_max=l_max
    !allocate(current_segment(0:small_n+1))
    !allocate(interpolation_table(0:length-1,0:3,0:1)) !There is cubic interpolation.
    allocate(interpolation_table(0:length-1,0:5,0:1)) !Lets try quintic interpolation.
    allocate(integration_table(0:length-1,0:small_n,0:1))
    
    !*********ROOT OF THE PROBLEM************
    ! All the (N) k values, well within the Nyquist frequency
    do i=0,length-1
        k_values(i)=(2*pi)*((i+(1.0/2.0))/(length-(1.0/2.0)))/(2*delta_r)
    end do
    
    !Call for the dcst_transforms ,and calculate for intepolating also,the differentiated
    !transform and save it in 3 d matrix. sin,cos ;sin',cos' ;sin'',cos''
    !We will be going with P=3 degree interpolation.
    call generate_DCST_transforms(length,delta_r,signal_array,dcst_transforms)
    
    !Getting the coefficient of local interpolation 
    call segment_coefficient(length,delta_r,k_values,signal_array,dcst_transforms,interpolation_table)
    !Have to integrate segments.
    call segment_integration(small_n,length,k_values,interpolation_table,integration_table)
    
    !Now we have to do the weighted average over all n to get the transform for a 
    !particular "K", for this particular "l" !!!!!!!! (Remember).
    print *,"Small _n",small_n
    do i=0,length-1
        transform_array(i)=0
        do j=0,small_n
            if(mod(l,2)==0)then
                transform_array(i)=transform_array(i)+((-1)**j)*double_factorial(2*small_n+2*j-1)/&
                    (factorial(2*j))/(double_factorial(2*small_n-2*j))/(k_values(i)**(2*j+1))*&
                    integration_table(i,j,0)
            else if(mod(l,2)/=0)then
                !print *,j,"ODD"
                !print *,transform_array(i)
                transform_array(i)=transform_array(i)+((-1)**j)*double_factorial(2*small_n+2*j+1)/&
                    (factorial(2*j+1))/double_factorial(2*small_n-2*j)/(k_values(i)**(2*j+2))*&
                    integration_table(i,j,1)
                !print *,transform_array(i),k_values(i)
            end if
        end do
    end do
    
end subroutine SpBessel_transform

!Integrating the segments, for different n and different k in one go.
subroutine segment_integration (small_n,length,k_values,interpolation_table,integration_table)
    implicit none
    integer,intent(in) :: small_n,length
    real*8,dimension(0:length-1) :: k_values
    complex(8),dimension(0:length-1,0:small_n,0:1),intent(out) :: integration_table
    !complex(8),dimension(0:length-1,0:3,0:1),intent(in) :: interpolation_table
    complex(8),dimension(0:length-1,0:5,0:1),intent(in) :: interpolation_table
    
    integer :: i,j,k
    do i=0,length-1
        do j=0,small_n
            if(i==0)then
                !EVEN l (00:gamma2,10:g4,20:g6)
                integration_table(i,j,0)=interpolation_table(i,0,0)*(k_values(i)**(2*j+1))/&
                    (2*j+1)-interpolation_table(i,1,0)*(k_values(i)**(2*j+3))/(2*(2*j+3))+&
                    interpolation_table(i,2,0)*(k_values(i)**(2*j+5))/(24*(2*j+5))
                !ODD l (01:g3,11:g5,21:g7)
                integration_table(i,j,1)=interpolation_table(i,0,1)*(k_values(i)**(2*j+1+2))/&
                    (2*j+1+2)-interpolation_table(i,1,1)*(k_values(i)**(2*j+1+4))/(6*(2*j+1+4))+&
                    interpolation_table(i,2,1)*(k_values(i)**(2*j+1+6))/(120*(2*j+1+6))
            else
                !EVEN
                integration_table(i,j,0)=integration_table(i-1,j,0)+&
                         (integrate_func(i,j,0,length,k_values,interpolation_table))
                !ODD
                integration_table(i,j,1)=integration_table(i-1,j,1)+&
                         (integrate_func(i,j,1,length,k_values,interpolation_table))
            end if
        end do
    end do
    
    print *,"Moving inside integrating the segment"
    do k=0,1
        do i=0,length-1
            do j=0,small_n
                !print *,k,i,j,integration_table(i,j,k)
            end do
        end do
    end do
    print *,"Coming out of integration segment"
end subroutine segment_integration

!Just to ease off above calcuation.
complex function integrate_func(i,j,OE_flag,length,k_values,interpolation_table) result(final)
    !E=0,O=1 : the OE flag.
    implicit none
    integer,intent(in) :: i,j,OE_flag,length
    real*8,dimension(0:length-1),intent(in) :: k_values
    !complex(8),dimension(0:length-1,0:3,0:1),intent(in) :: interpolation_table
    complex(8),dimension(0:length-1,0:5,0:1),intent(in) :: interpolation_table
    integer :: power
    real*8 :: k0,k1
    complex(8) :: coff0,coff1,coff2,coff3,coff4,coff5
    
    power=2*j+OE_flag
    k0=k_values(i-1)! i=1: k values zeroth entry is k0. there i is refering to interval k0-k1
    k1=k_values(i)
    
    !Finding the coefficient of the interpolation function.
    coff0=interpolation_table(i,0,OE_flag)-interpolation_table(i,1,OE_flag)*k0&
            +interpolation_table(i,2,OE_flag)*(k0**2)-interpolation_table(i,3,OE_flag)*(k0**3)&
            +interpolation_table(i,4,OE_flag)*(k0**4)-interpolation_table(i,5,OE_flag)*(k0**5)
    coff1=interpolation_table(i,1,OE_flag)-2*interpolation_table(i,2,OE_flag)*(k0)&
            +3*interpolation_table(i,3,OE_flag)*(k0**2)-4*interpolation_table(i,4,OE_flag)*(k0**3)&
            +5*interpolation_table(i,5,OE_flag)*(k0**4)
    coff2=interpolation_table(i,2,OE_flag)-3*interpolation_table(i,3,OE_flag)*(k0)&
            +6*interpolation_table(i,4,OE_flag)*(k0**2)-10*interpolation_table(i,5,OE_flag)*(k0**3)
    coff3=interpolation_table(i,3,OE_flag)-4*interpolation_table(i,4,OE_flag)*(k0)&
            +10*interpolation_table(i,5,OE_flag)*(k0**2)
    coff4=interpolation_table(i,4,OE_flag)-5*interpolation_table(i,5,OE_flag)*k0
    coff5=interpolation_table(i,5,OE_flag)
    
    print *,"Printing the coefficient of interpolation"
    !Integrating and putting the Limits k0 to k1
    final=coff0/(power+1)*(k1**(power+1)-k0**(power+1))+&
            coff1/(power+2)*(k1**(power+2)-k0**(power+2))+&
            coff2/(power+3)*(k1**(power+3)-k0**(power+3))+&
            coff3/(power+4)*(k1**(power+4)-k0**(power+4))+&
            coff4/(power+5)*(k1**(power+5)-k0**(power+5))+&
            coff5/(power+6)*(k1**(power+6)-k0**(power+6))
    
end function  integrate_func

!get integral done for all j upto small_n(not here)
!Just interpolating the functions between the segment using cubic polynomial
!So findgin the coefficient of teh interpolation.
subroutine segment_coefficient(length,delta_r,k_values,signal,dcst_transforms,interpolation_table)
    implicit none
    integer :: i,j,k
    integer,intent(in) :: length
    real*8 :: o_gamma,temp !For first segment coefficient.
    real*8,dimension(0:length-1,0:5) :: r_power
    
    real*8,intent(in) :: delta_r
    real*8,dimension(0:length-1),intent(in) :: k_values
    complex(8),dimension(0:length-1),intent(in) :: signal
    complex(8),dimension(0:1,0:length-1,0:2),intent(in) :: dcst_transforms
    !complex(8),dimension(0:length-1,0:3,0:1),intent(out) :: interpolation_table
    complex(8),dimension(0:length-1,0:5,0:1),intent(out) :: interpolation_table
    
    !Construct r-array for the zeros-k1 interpolation
    do i=0,length-1
        temp=(i+(1.0/2.0))*delta_r
        !for both even and odd case.
        r_power(i,0:5)=(/temp**2,temp**4,temp**6,temp**3,temp**5,temp**7/)
    end do
    
    print *,""
    print *,"Entering the domain of coefficient"
    
    print *,"checking r_values"
    do i =0,length-1
        !print *,r_power(i,0),r_power(i,1),r_power(i,2),r_power(i,3),r_power(i,4),r_power(i,5)
    end do
    
    do i=0,length-1
        if(i==0)then
            !Only three coefficient are there in case of this interval of integ.
            !There is difference in position of cos and sine here and in dcst matrix
            !There 0: sine(accordingg to 1st val) here 0:cos (according to n)
            !Even case
            !delta is missiing . Multiply later.
            interpolation_table(i,0,0)=sum(signal*r_power(0:length-1,0))
            interpolation_table(i,1,0)=sum(signal*r_power(0:length-1,1))
            interpolation_table(i,2,0)=sum(signal*r_power(0:length-1,2))
            interpolation_table(i,3,0)=cmplx(0.0,0.0)
            interpolation_table(i,4,0)=cmplx(0.0,0.0)
            interpolation_table(i,5,0)=cmplx(0.0,0.0)
            !Odd case
            interpolation_table(i,0,1)=sum(signal*r_power(0:length-1,3))
            interpolation_table(i,1,1)=sum(signal*r_power(0:length-1,4))
            interpolation_table(i,2,1)=sum(signal*r_power(0:length-1,5))
            interpolation_table(i,3,1)=cmplx(0.0,0.0)
            interpolation_table(i,4,1)=cmplx(0.0,0.0)
            interpolation_table(i,5,1)=cmplx(0.0,0.0)
        else
            !Even Case
            !delta from previous transform sin /cosine is missing.Multiply ltr.
            interpolation_table(i,0,0)=dcst_transforms(1,i-1,0)
            interpolation_table(i,1,0)=dcst_transforms(1,i-1,1)
            interpolation_table(i,2,0)=dcst_transforms(1,i-1,2)
            !interpolation_table(i,2,0)=3*(dcst_transforms(1,i,0)-dcst_transforms(1,i-1,0))/&
                    !(k_values(i)-k_values(i-1))**2-(dcst_transforms(1,i,1)-dcst_transforms(1,i-1,1))/&
                    !(k_values(i)-k_values(i-1))
            !interpolation_table(i,3,0)=(-2)*(dcst_transforms(1,i,0)-dcst_transforms(1,i-1,0))/&
                    !(k_values(i)-k_values(i-1))**3+(dcst_transforms(1,i,1)+dcst_transforms(1,i-1,1))/&
                    !(k_values(i)-k_values(i-1))**2
            interpolation_table(i,3,0)=10*(dcst_transforms(1,i,0)-dcst_transforms(1,i-1,0))/(k_values(i)-k_values(i-1))**3-&
                                        (4*dcst_transforms(1,i,1)+6*dcst_transforms(1,i-1,1))/(k_values(i)-k_values(i-1))**2+&
                                        (dcst_transforms(1,i,2)-3*dcst_transforms(1,i-1,2))/(2*(k_values(i)-k_values(i-1)))
            interpolation_table(i,4,0)=(-15)*(dcst_transforms(1,i,0)-dcst_transforms(1,i-1,0))/(k_values(i)-k_values(i-1))**4+&
                                        (7*dcst_transforms(1,i,1)+8*dcst_transforms(1,i-1,1))/(k_values(i)-k_values(i-1))**3-&
                                        (2*dcst_transforms(1,i,1)-3*dcst_transforms(1,i-1,1))/(2*(k_values(i)-k_values(i-1))**2)
            interpolation_table(i,5,0)=6*(dcst_transforms(1,i,0)-dcst_transforms(1,i-1,0))/(k_values(i)-k_values(i-1))**5-&
                                        3*(dcst_transforms(1,i,1)+dcst_transforms(1,i-1,1))/(k_values(i)-k_values(i-1))**4+&
                                        (dcst_transforms(1,i,1)-dcst_transforms(1,i-1,1))/(2*(k_values(i)-k_values(i-1))**3)
            
            !Odd case
            interpolation_table(i,0,1)=dcst_transforms(0,i-1,0)
            interpolation_table(i,1,1)=dcst_transforms(0,i-1,1)
            interpolation_table(i,2,1)=dcst_transforms(0,i-1,2)
            !interpolation_table(i,2,1)=3*(dcst_transforms(0,i,0)-dcst_transforms(0,i-1,0))/&
                    !(k_values(i)-k_values(i-1))**2-(dcst_transforms(0,i,1)-dcst_transforms(0,i-1,1))/&
                    !(k_values(i)-k_values(i-1))
            !interpolation_table(i,3,1)=(-2)*(dcst_transforms(0,i,0)-dcst_transforms(0,i-1,0))/&
                    !(k_values(i)-k_values(i-1))**3+(dcst_transforms(0,i,1)+dcst_transforms(0,i-1,1))/&
                    !(k_values(i)-k_values(i-1))**2
            interpolation_table(i,3,1)=10*(dcst_transforms(0,i,0)-dcst_transforms(0,i-1,0))/(k_values(i)-k_values(i-1))**3-&
                                        (4*dcst_transforms(0,i,1)+6*dcst_transforms(0,i-1,1))/(k_values(i)-k_values(i-1))**2+&
                                        (dcst_transforms(0,i,2)-3*dcst_transforms(0,i-1,2))/(2*(k_values(i)-k_values(i-1)))
            interpolation_table(i,4,1)=-15*(dcst_transforms(0,i,0)-dcst_transforms(0,i-1,0))/(k_values(i)-k_values(i-1))**4+&
                                        (7*dcst_transforms(0,i,1)+8*dcst_transforms(0,i-1,1))/(k_values(i)-k_values(i-1))**3-&
                                        (2*dcst_transforms(0,i,1)-3*dcst_transforms(1,i-1,1))/(2*(k_values(i)-k_values(i-1))**2)
            interpolation_table(i,5,1)=6*(dcst_transforms(0,i,0)-dcst_transforms(0,i-1,0))/(k_values(i)-k_values(i-1))**5-&
                                        3*(dcst_transforms(0,i,1)+dcst_transforms(0,i-1,1))/(k_values(i)-k_values(i-1))**4+&
                                        (dcst_transforms(0,i,1)-dcst_transforms(0,i-1,1))/(2*(k_values(i)-k_values(i-1))**3)
        end if
    end do
    
    do i=0,1
        do j=0,length-1
            do k=0,5
                !print *,i,j,k,interpolation_table(j,k,i)
            end do
            !print *,"XXXXXXXXX"
        end do
        !print *,"XXXXXXXXXXXXXXXX"
    end do
    print *,"Ending and leaving interpolation table subroutine."
    
end subroutine segment_coefficient

subroutine generate_DCST_transforms(length,delta_r,signal,dcst_transforms)
    use FFT
    implicit none
    integer :: i,j,k
    integer,intent(in) ::length
    real*8,intent(in) :: delta_r
    complex(8),dimension(0:1,0:length-1,0:2),intent(out) :: dcst_transforms
    complex(8),dimension(0:length-1),intent(in) :: signal
    
    real*8,dimension(0:length-1,0:length-1) :: transform_matrix
    
    !Now generating the transform matrix to do the transform of the data.
    !The sequence is maintained. THe derivative of first layer is just behind it.
    
    !First filling the sine part.
    call make_transform_matrix(0,length,transform_matrix)
    call customized_DCST1D(length,delta_r,2,signal,transform_matrix,&
                            dcst_transforms(0,0:length-1,0))
    call customized_DCST1D(length,delta_r,3,signal,transform_matrix,&
                            dcst_transforms(1,0:length-1,1))
    dcst_transforms(1,0:length-1,1)=dcst_transforms(1,0:length-1,1)*(-1)
    call customized_DCST1D(length,delta_r,4,signal,transform_matrix,&
                            dcst_transforms(0,0:length-1,2))  
    dcst_transforms(0,0:length-1,2)=dcst_transforms(0,0:length-1,2)*(-1)
    
    ! Now Filling the cos part
    call make_transform_matrix(1,length,transform_matrix)
    call customized_DCST1D(length,delta_r,2,signal,transform_matrix,&
                            dcst_transforms(1,0:length-1,0))
    call customized_DCST1D(length,delta_r,3,signal,transform_matrix,&
                            dcst_transforms(0,0:length-1,1))
    call customized_DCST1D(length,delta_r,4,signal,transform_matrix,&
                            dcst_transforms(1,0:length-1,2))
    dcst_transforms(1,0:length-1,2)=dcst_transforms(1,0:length-1,2)*(-1)
    
    print *,"Printing Dcst Matrix"
    do i=0,1
        do j=0,2
            do k=0,length-1
                !print *,i,j,k,dcst_transforms(i,k,j)
            end do
            !print *,"******"
        end do
        !print *,"*************"
    end do
    print *,"Going out of DCST SUBROUTINE."
end subroutine generate_DCST_transforms





!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!|||||||||||||||||||||||| SPHERICAL TRANSFORM ||||||||||||||||||||||||||||||||
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


!This will be used to complete teh spherical Harmnics of particular layer 
!ie. Spherical Layer.
subroutine SPH_transform(lat_length,bWidth,signal_r)
    implicit none
    integer,intent(in) :: lat_length,bWidth
    complex(8),dimension(0:lat_length-1,0:lat_length-1),intent(inout) :: signal_r 
    
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
    !multiplying the final normalization constant after the 
    !constants directly into the transform matrix according to l and m.
    signal_r=signal_r*((2*pi)**(1.0/2.0)/(2*bWidth))
    print *,"KIH"
    
end subroutine SPH_transform

!THis will be used in above to do the FFT of signal array.
!Its for more abstraction.
subroutine do_FFT(lat_length,bWidth,signal_r)
    use FFT
    implicit none
    integer,intent(in) :: lat_length,bWidth
    complex(8),dimension(0:lat_length-1,0:lat_length-1),intent(inout):: signal_r
    
    integer :: i,depth
    integer,allocatable,dimension(:,:) :: dump_credential
    real*8,parameter :: pi=4*atan(1.0),tau=2*pi
    complex(8) :: omega
    complex(8),dimension(0:lat_length-1) :: temp_row
    
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
    complex(8),dimension(0:lat_length-1,0:lat_length-1),intent(inout) :: signal
    
    integer :: i
    real*8,allocatable,dimension(:,:) :: transform_matrix
    real*8,dimension(0:lat_length-1) :: weight_vector
    complex(8),allocatable,dimension(:) :: transform
    
    
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
    complex(8),dimension(0:lat_length-1,0:lat_length-1),intent(inout) :: signal
    complex(8),dimension(i:bWidth-1),intent(in) :: transform
    
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
    real*8,parameter :: pi=4*atan(1.0)
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
    !******IMP******
    !Adding the normalization constant here only rather than multiplying it later in teh signal.
    transform_matrix(abs(m),0:lat_length-1)=Pm_m(0:lat_length-1)*(((2*abs(m)+1)/(4*pi)/factorial(2*abs(m)))**(1.0/2.0))
    if(m==bWidth-1)then
        return
    end if
    transform_matrix(abs(m)+1,0:lat_length-1)=Pm_m1(0:lat_length-1)*(((2*(abs(m)+1)+1)/(4*pi)/factorial(2*abs(m)+1))**(1.0/2.0))
    !highest bandwidth ie l=bWidth-1
    ex_recurs=(bWidth-1)-abs(m)-1
    !print *,ex_recurs
    do i=1,ex_recurs
        Pm_temp=(arg*(2*(abs(m)+i+1)-1)*Pm_m1-((abs(m)+i+1)+abs(m)-1)*Pm_m)/((abs(m)+i+1)-abs(m))
        transform_matrix(abs(m)+1+i,0:lat_length-1)=Pm_temp(0:lat_length-1)*(((2*(abs(m)+i+1)+1)/(4*pi)*&
                                                       factorial((abs(m)+i+1)-abs(m))/factorial((abs(m)+i+1)-abs(m)))**(1.0/2.0))
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
        if(mod(num,2)/=0 .and. mod(i,2)/=0)then
            temp=temp*i
        else if(mod(num,2)==0 .and. mod(i,2)==0)then
            temp=temp*i
        end if
    end do
    double_factorial=temp
    
end function double_factorial

integer*8 function factorial(num)
    implicit none
    integer,intent(in) :: num
    integer :: i
    
    factorial=1
    do i=1,num
        factorial=factorial*i
    end do
end function factorial

end module BSPHT