function [record_obs U_obs] = forward_waveflied( nx,nz,dx,dz,pml,nxpml,nzpml,nt,dt,Vp_PML,atten_x,atten_z,ishot,nz_location,amp,riker_fsource)
nx_pml=nx+2*nxpml;
nz_pml=nz+2*nzpml;
U1_pre=zeros(nz_pml,nx_pml);
U2_pre=zeros(nz_pml,nx_pml);
U1_now=zeros(nz_pml,nx_pml);
U2_now=zeros(nz_pml,nx_pml);
U_record=zeros(nt,nx*nz);
 record=zeros(nt,nx);
nx_location=4*ishot-2;
fsource_location=zeros(nz_pml,nx_pml);
fsource_location(nz_location+pml+1,nx_location+pml+1)=amp;
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
for(it=1:1:nt)
    %it
     %[Uxx,Uzz] = difference_ten(nx,nz,dx,dz,pml,U_1);
    
%for(i=nStartz:1:nEndz) 	
    
	%for(k=nStartx:1:nEndx) 
			%Uxx(i,k) = coef_1*(omiga0*U(i,k)+ omiga1*(U(i,k+1)+U(i,k-1))+ omiga2*(U(i,k+2)+U(i,k-2)) + omiga3*(U(i,k+3)+U(i,k-3))+ omiga4*(U(i,k+4)+U(i,k-4))+ omiga5*(U(i,k+5)+U(i,k-5)));
            %Uzz(i,k) = coef_2*(omiga0*U(i,k)+ omiga1*(U(i+1,k)+U(i-1,k))+ omiga2*(U(i+2,k)+U(i-2,k)) + omiga3*(U(i+3,k)+U(i-3,k))+ omiga4*(U(i+4,k)+U(i-4,k))+ omiga5*(U(i+5,k)+U(i-5,k))); 
   % end

%end
     Uxx(:,nStartx:nEndx)=coef_1*(omiga0*U(:,nStartz:nEndz)+ omiga1*(U(:,nStartz-1:nEndz-1)+U(:,nStartz+1:nEndz+1))+ omiga2*(U(:,nStartz-2:nEndz-2)+U(:,nStartz+2:nEndz+2)) + omiga3*(U(:,nStartz-3:nEndz-3)+U(:,nStartz+3:nEndz+3))+ omiga4*(U(:,nStartz-4:nEndz-4)+U(:,nStartz+4:nEndz+4))+ omiga5*(U(:,nStartz-5:nEndz-5)+U(:,nStartz+5:nEndz+5)));
     Uzz(nStartz:nEndz,:)=coef_2*(omiga0*U(nStartz:nEndz,:)+ omiga1*(U(nStartz-1:nEndz-1,:)+U(nStartz+1:nEndz+1,:))+ omiga2*(U(nStartz-2:nEndz-2,:)+U(nStartz+2:nEndz+2,:)) + omiga3*(U(nStartz-3:nEndz-3,:)+U(nStartz+3:nEndz+3,:))+ omiga4*(U(nStartz-4:nEndz-4,:)+U(nStartz+4:nEndz+4,:))+ omiga5*(U(nStartz-5:nEndz-5,:)+U(nStartz+5:nEndz+5,:)));

			U1_next=1.0./(1+atten_x*dt).*((2-atten_x.*atten_x*dt*dt).*U1_now+(atten_x*dt-1).*U1_pre+ Vp_PML.*Vp_PML.*Uxx*dt*dt);

			U2_next=1.0./(1+atten_z*dt).*((2-atten_z.*atten_z*dt*dt).*U2_now+(atten_z*dt-1).*U2_pre+ Vp_PML.*Vp_PML.*Uzz*dt*dt);
			U= U1_next+U2_next+fsource_location*riker_fsource(it);
            
            U_save=U(nzpml+1:nzpml+nz,nxpml+1:nxpml+nx);
        
  
        
		
	     record(it,:)=U_save(nz_location,:);
               
         
                U1_pre=U1_now;
                U1_now=U1_next;
				
                U2_pre=U2_now;
                U2_now=U2_next;
                    
               U_record(it,:)=U_save(:);
               
              %U_1=U(:);
   
end
     
             for(ix=1:1:nx)
              
                  record(1:floor((abs(ix-nx_location)*dx)/(Vp_PML(pml+2,pml+1)*dt))+250,ix)=0;
              
             end
              record_obs=reshape(record,nt,nx,1);
              U_obs=reshape(U_record,nt,nx*nz,1);
              
     
%figure;
%imagesc(reshape(U_record(200,:),nz,nx))
clear U1_pre U1_now  U1_next U2_pre U2_now U2_next Uxx Uzz U_save  record  atten_x  atten_z U_record VP_PML
               
end 
