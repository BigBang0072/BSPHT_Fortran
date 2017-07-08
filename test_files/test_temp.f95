program test_temp
    implicit none
    integer :: i
    integer,allocatable,dimension(:) :: ak
    
    allocate(ak(2))
    do i=1,2
        ak(i)=i
    end do
    do i=1,2
        print *,ak(i)
    end do
    
    deallocate(ak)
    allocate(ak(3))
    
    do i=1,3
        ak(i)=i
    end do
    do i=1,3
        print *,ak(i)
    end do
    
    deallocate(ak)
    
    
end program test_temp