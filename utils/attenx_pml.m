function [atten_x,atten_z ] =attenx_pml( nx,nz,PML,par)
   nx_PML=nx+2*PML;
   nz_PML=nz+2*PML;
   atten_x=zeros(nz_PML,nx_PML);
   atten_z=zeros(nz_PML,nx_PML);
   
 for(ix=1:1:nx_PML)
	
    for (iz=1:1:nz_PML)
		
        if (ix<=PML)
			
				atten_x(iz,ix)=par*(PML-ix)*(PML-ix);
			
        elseif (ix>=nx+PML)
			
				atten_x(iz,ix)=par*(ix-nx-PML)*(ix-nx-PML);
			
         end
     end
 end
	
	
	for (ix=1:1:nx_PML)
	
		for (iz=1:1:nz_PML)
		
			if (iz<=PML)
			
				atten_z(iz,ix)=par*(PML-iz)*(PML-iz);
			
            elseif (iz>=nz+PML)
			
				atten_z(iz,ix)=par*(iz-nz-PML)*(iz-nz-PML);
            end
            
        end
    end
  %figure;
   %imagesc(atten_x);
    %figure;
   %imagesc(atten_z);
end

