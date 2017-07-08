for i in range(8):
	for k in range(8):
		a[i,k]=math.sqrt(math.sin(math.pi*(2*i+1)/(4*4))**2)*(math.cos(2*math.pi*k/(2*4))+math.sin(2*math.pi*k/(2*4))j)
        
        
        
    for i in range(8):
	for k in range(8):
		a[i,k]=-math.sqrt(math.sin(math.pi*(2*i+1)/(4*4))**2)
    
    for i in range(8):
	for k in range(8):
		a[i,k]=a[i,k]*(math.cos(2*math.pi*k/(2*4))+math.sin(2*math.pi*k/(2*4))*1j)
        
b=np.zeros((8,8),dtype=np.complex_)
for i in range(8):
	b[i,:]=np.fft.fft(a[i,:])
    
    
for i in range(16):
	for k in range(16):
		a[i,k]=(105.0/8.0)*(math.sin(math.pi*(2*i+1)/(4*8))**2)*(5*math.cos(math.pi*(2*i+1)/(4*8))+3*math.cos(math.pi*3*(2*i+1)/(4*8)))
        
        
for i in range(16):
	for k in range(16):
		a[i,k]=a[i,k]*(math.cos(2*math.pi*3*k/(2*8))+math.sin(2*math.pi*3*k/(2*8))*1j)
        
        
