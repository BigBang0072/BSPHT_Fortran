program naive
    implicit none
    real,parameter :: pi=4*atan(1.0)
    !real,external :: bessjy
    real :: Jp,dumm1,dumm2,dumm3
    
    
    call bessjy(2*pi*0.5*0.5,1,Jp,dumm1,dumm2,dumm3)
    print *,Jp
    
end program naive