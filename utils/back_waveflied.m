function [ U_record] = back_waveflied( nx,nz,dx,dz,pml,nxpml,nzpml,nt,dt,Vp_PML,atten_x,atten_z,nz_location,amp, riker_fsource_back)
nx_pml=nx+2*nxpml;
nz_pml=nz+2*nzpml;
U1_pre=zeros(nz_pml,nx_pml);
U2_pre=zeros(nz_pml,nx_pml);
U1_now=zeros(nz_pml,nx_pml);
U2_now=zeros(nz_pml,nx_pml);
%U_1=zeros(nx_pml*nz_pml,1);
%record=zeros(nt,nx);
U_record=zeros(nt,nx*nz);
fsource_location_back=zeros(nz_pml,nx_pml);
fsource_location_back(nz_location+pml+1,pml+1:nx+pml)=amp;
%Vp_PML=  velocity_pml( nxpml,nzpml,nx,nz,vel );
%Vp_PML_real=  myfun2( nxpml,nzpml,nx,nz,vel_1 );
%[atten_x,atten_z] = attenx_pml( nx,nz,pml,par);
   coef_1=1.0/(2.0*dx*dx); 
    coef_2=1.0/(2.0*dz*dz); 
	omiga0=-5.8544444444;
	omiga1=3.33333333333;
	omiga2=-0.4761904762;
	omiga3=0.0793650794;
	omiga4= -0.0099206349;
	omiga5=0.0006349206;
	nStartx =6;
	nEndx =nx_pml-6;
	nStartz =6;
	nEndz =nz_pml-6;
    U=zeros(nz_pml,nx_pml);
    Uxx=zeros(nz_pml,nx_pml);
    Uzz=zeros(nz_pml,nx_pml); 
for(it=nt:-1:1)
    %it
    fsource_location_back(nz_location+pml+1,pml+1:nx+pml)=riker_fsource_back(it,:);
      %[Uxx,Uzz] = difference_ten(nx,nz,dx,dz,pml,U_1);
     Uxx(:,nStartx:nEndx)=coef_1*(omiga0*U(:,nStartz:nEndz)+ omiga1*(U(:,nStartz-1:nEndz-1)+U(:,nStartz+1:nEndz+1))+ omiga2*(U(:,nStartz-2:nEndz-2)+U(:,nStartz+2:nEndz+2)) + omiga3*(U(:,nStartz-3:nEndz-3)+U(:,nStartz+3:nEndz+3))+ omiga4*(U(:,nStartz-4:nEndz-4)+U(:,nStartz+4:nEndz+4))+ omiga5*(U(:,nStartz-5:nEndz-5)+U(:,nStartz+5:nEndz+5)));
     Uzz(nStartz:nEndz,:)=coef_2*(omiga0*U(nStartz:nEndz,:)+ omiga1*(U(nStartz-1:nEndz-1,:)+U(nStartz+1:nEndz+1,:))+ omiga2*(U(nStartz-2:nEndz-2,:)+U(nStartz+2:nEndz+2,:)) + omiga3*(U(nStartz-3:nEndz-3,:)+U(nStartz+3:nEndz+3,:))+ omiga4*(U(nStartz-4:nEndz-4,:)+U(nStartz+4:nEndz+4,:))+ omiga5*(U(nStartz-5:nEndz-5,:)+U(nStartz+5:nEndz+5,:)));

			U1_next=1.0./(1+atten_x*dt).*((2-atten_x.*atten_x*dt*dt).*U1_now+(atten_x*dt-1).*U1_pre+ Vp_PML.*Vp_PML.*Uxx*dt*dt);

			U2_next=1.0./(1+atten_z*dt).*((2-atten_z.*atten_z*dt*dt).*U2_now+(atten_z*dt-1).*U2_pre+ Vp_PML.*Vp_PML.*Uzz*dt*dt);
			U = U1_next+U2_next+fsource_location_back;
            
            U_save=U(nzpml+1:nzpml+nz,nxpml+1:nxpml+nx);
        
  
               %for (ix=1:1:nx)
		
			  % record(it,ix)=U_save(nz_location,ix);
               
              % end
                U1_pre=U1_now;
                U1_now=U1_next;
				
                U2_pre=U2_now;
                U2_now=U2_next;
                    
               U_record(it,:)=U_save(:);
               
          
          
end 