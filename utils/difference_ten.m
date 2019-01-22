function[D1,D2] = difference_ten(nx,nz,dx,dz,pml,U_1) 
    coef_1=1.0/(2.0*dx*dx); 
    coef_2=1.0/(2.0*dz*dz); 
	omiga0=-5.8544444444;
	omiga1=3.33333333333;
	omiga2=-0.4761904762;
	omiga3=0.0793650794;
	omiga4= -0.0099206349;
	omiga5=0.0006349206;
    nx_pml=nx+2*pml;
	nz_pml=nz+2*pml;
	nStartx =6;
	nEndx =nx_pml-6;
	nStartz =6;
	nEndz =nz_pml-6;
    U=reshape(U_1,nz_pml,nx_pml);
    D1=zeros(nz_pml,nx_pml);
    D2=zeros(nz_pml,nx_pml);
for(i=nStartz:1:nEndz) 	
    
	for(k=nStartx:1:nEndx) 
			D1(i,k) = coef_1*(omiga0*U(i,k)+ omiga1*(U(i,k+1)+U(i,k-1))+ omiga2*(U(i,k+2)+U(i,k-2)) + omiga3*(U(i,k+3)+U(i,k-3))+ omiga4*(U(i,k+4)+U(i,k-4))+ omiga5*(U(i,k+5)+U(i,k-5)));
            D2(i,k) = coef_2*(omiga0*U(i,k)+ omiga1*(U(i+1,k)+U(i-1,k))+ omiga2*(U(i+2,k)+U(i-2,k)) + omiga3*(U(i+3,k)+U(i-3,k))+ omiga4*(U(i+4,k)+U(i-4,k))+ omiga5*(U(i+5,k)+U(i-5,k))); 
    end

end

end
