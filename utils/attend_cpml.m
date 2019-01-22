function [atten_ax,atten_az atten_bx,atten_bz atten_dx,atten_dz ] =attend_cpml( nx,nz,PML,attend_a0, attend_b0, attend_d0)
   nx_PML=nx+2*PML;
   nz_PML=nz+2*PML;
   atten_ax=zeros(nz_PML,nx_PML);
   atten_az=zeros(nz_PML,nx_PML);
   atten_bx=zeros(nz_PML,nx_PML);
   atten_bz=zeros(nz_PML,nx_PML);
   atten_dx=zeros(nz_PML,nx_PML);
   atten_dz=zeros(nz_PML,nx_PML);
 for(ix=1:1:nx_PML)
	
    for (iz=1:1:nz_PML)
		
        if (ix<=PML)
			
				atten_dx(iz,ix)=attend_d0*(PML-ix)*(PML-ix)/(PML^2);
                atten_bx(iz,ix)=1+(attend_b0-1)*(PML-ix)*(PML-ix)/(PML^2);
			    atten_ax(iz,ix)=attend_a0*(1-(PML-ix)/PML);
        elseif (ix>=nx+PML)
			
				atten_dx(iz,ix)=attend_d0*(ix-nx-PML)*(ix-nx-PML)/(PML^2);
			   atten_bx(iz,ix)=1+(attend_b0-1)*(ix-nx-PML)*(ix-nx-PML)/(PML^2);
               atten_ax(iz,ix)=attend_a0*(1-(ix-nx-PML)/PML);
         end
     end
 end
	
	
	for (ix=1:1:nx_PML)
	
		for (iz=1:1:nz_PML)
		
			if (iz<=PML)
			
				atten_dz(iz,ix)=attend_d0*(PML-iz)*(PML-iz)/(PML^2);
                atten_bz(iz,ix)=1+(attend_b0-1)*(PML-iz)*(PML-iz)/(PML^2);
                atten_az(iz,ix)=attend_a0*(1-(PML-iz)/PML);
			
            elseif (iz>=nz+PML)
			
				atten_dz(iz,ix)=attend_d0*(iz-nz-PML)*(iz-nz-PML)/(PML^2);
                atten_bz(iz,ix)=1+(attend_b0-1)*(iz-nz-PML)*(iz-nz-PML)/(PML^2);
               atten_az(iz,ix)=attend_a0*(1-(iz-nz-PML)/PML);
            end
            
        end
    end
  %figure;
   %imagesc(atten_x);
    %figure;
   %imagesc(atten_z);
end
