function [ ob_data_Vz Vz_all_fft ] = time_frequency_forward(nx,nz, pml, dx, dz, Vp_PML )

m=4;
kappa_max=1;
alpha_max=0;
bj=1;
MP_Tyz_ya=ones(nz_pml,nx_pml+1); 

sig_max=bj.*(m+1)./(150.*pi.*dz).*sqrt(max(max(Vp_PML)));
sigy=sig_max.*repmat(([pml-1:-1:1 1:pml-1]./pml).^m,[nz_pml,1]); %pml������1
kappay=1+(kappa_max-1).*([pml-1:-1:1 1:pml-1]./pml).^m; %pml������1
kappay=repmat(kappay,[nz_pml,1]);
alphay=alpha_max.*([pml-1:-1:1 1:pml-1]./pml); %=0 %pml������1
alphay=repmat(alphay,[nz_pml,1]);
MP_Tyz_ya(:,[2:pml nz_pml-pml+2:nz_pml])=1./kappay;
MP_Tyz_yb=exp(-((sigy./kappay)+alphay));
MP_Tyz_yc=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Tyz_yb-1);
clear  sig_max sigy kappay alphay


MP_Txz_xa=ones(nz_pml+1,nx_pml);

sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(Vp_PML)));
sigx=sig_max.*repmat(([pml-1:-1:1 1:pml-1]'./pml).^m,[1,nx_pml]); %pml������1
kappax=1+(kappa_max-1).*([pml-1:-1:1 1:pml-1]'./pml).^m; %pml������1
kappax=repmat(kappax,[1,nx_pml]);
alphax=alpha_max.*([pml-1:-1:1 1:pml-1]'./pml); %=0 %pml������1
alphax=repmat(alphax,[1,nx_pml]);
MP_Txz_xa([2:pml nz_pml-pml+2:nz_pml],:)=1./kappax; 
MP_Txz_xb=exp(-((sigx./kappax)+alphax));
MP_Txz_xc=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Txz_xb-1);


MP_Vz=dt.*Vp_PML.^2;
MP_Vz_xa=ones(size(Vp_PML)); 
MP_Vz_ya=ones(size(Vp_PML));

sig_max=bj.*(m+1)./(150.*pi.*dz).*sqrt(max(max(Vp_PML)));
sigx=sig_max.*repmat((([pml:-1:1 1:pml]'-0.5)./pml).^m,[1,nx_pml]);
kappax=1+(kappa_max-1).*(([pml:-1:1 1:pml]'-0.5)./pml).^m;
kappax=repmat(kappax,[1,nx_pml]);
alphax=alpha_max.*(([pml:-1:1 1:pml]'-0.5)./pml); %=0
alphax=repmat(alphax,[1,nx_pml]);
MP_Vz_xa([1:pml nz_pml-pml+1:nz_pml],:)=1./kappax;
MP_Vz_xb=exp(-((sigx./kappax)+alphax));
MP_Vz_xc=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vz_xb-1);

sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(Vp_PML)));
sigy=sig_max.*repmat((([pml:-1:1 1:pml]-0.5)./pml).^m,[nz_pml,1]);
kappay=1+(kappa_max-1).*(([pml:-1:1 1:pml]-0.5)./pml).^m;
kappay=repmat(kappay,[nz_pml,1]);
alphay=alpha_max.*(([pml:-1:1 1:pml]-0.5)./pml); %=0
alphay=repmat(alphay,[nz_pml,1]);
MP_Vz_ya(:,[1:pml nx_pml-pml+1:nx_pml])=1./kappay; 
MP_Vz_yb=exp(-((sigy./kappay)+alphay));
MP_Vz_yc=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vz_yb-1);

clear sig_max sigx sigy alphay kappay

Tyz=zeros(nz_pml,nx_pml+1);
UV_Tyz_y=zeros(nz_pml,2*(pml-1));
Txz=zeros(nz_pml+1,nx_pml);
UV_Txz_x=zeros(2*(pml-1),nx_pml);
Vz=zeros(nz_pml,nx_pml);
UV_Vz_x=zeros(2*pml,nx_pml);
UV_Vz_y=zeros(nz_pml,2*pml);
Vz_all=zeros(nt,nz,nx);
Vz_all_fft=zeros(nt,nx*nz);
ob_data_Vz=zeros(length(recv_x),length(recv_y),nt);
%myavi=VideoWriter('2D_forward_Acoustic_waves.avi','Motion JPEG AVI');
%myavi.FrameRate=5;
%open(myavi);
%tic

for j=1:nt/2
	j
    % Vz(fsx,fsy)=Vz(fsx,fsy) + source(j); %obmode_1
%     Vz(fsx,fsy)=Vz(fsx,fsy) + source(:,j); % obmode_2
  %  Vz(fsx,fsy)=Vz(fsx,fsy) + repmat(source(:,j)',length(fsx),1); % obmode_3
     Vz(fsx,fsy)= source(j); 
	Txz(2:nz_pml,:)=Txz(2:nz_pml,:) + dt.*MP_Txz_xa(2:nz_pml,:).*(Vz(2:nz_pml,:)-Vz(1:nz_pml-1,:))./dz;
	UV_Txz_x=MP_Txz_xb.*UV_Txz_x + MP_Txz_xc.*(Vz([2:pml nz_pml-pml+2:nz_pml],:)-Vz([1:pml-1 nz_pml-pml+1:nz_pml-1],:))./dz;
	Txz([2:pml nz_pml-pml+2:nz_pml],:)=Txz([2:pml nz_pml-pml+2:nz_pml],:) + dt.*UV_Txz_x;
	Txz([1 nz_pml+1],:)=Txz([2 nz_pml],:);
		
	Tyz(:,2:nx_pml)=Tyz(:,2:nx_pml) + dt.*MP_Tyz_ya(:,2:nx_pml).*(Vz(:,2:nx_pml)-Vz(:,1:nx_pml-1))./dx;
	UV_Tyz_y=MP_Tyz_yb.*UV_Tyz_y + MP_Tyz_yc.*(Vz(:,[2:pml nx_pml-pml+2:nx_pml])-Vz(:,[1:pml-1 nx_pml-pml+1:nx_pml-1]))./dx;
	Tyz(:,[2:pml nx_pml-pml+2:nx_pml])=Tyz(:,[2:pml nx_pml-pml+2:nx_pml]) + dt.*UV_Tyz_y;
    Tyz(:,[1 nx_pml+1])=Tyz(:,[2 nx_pml]);
    
	Vz=Vz + MP_Vz.*MP_Vz_xa.*(Txz(2:nz_pml+1,:)-Txz(1:nz_pml,:))./dz + MP_Vz.*MP_Vz_ya.*(Tyz(:,2:nx_pml+1)-Tyz(:,1:nx_pml))./dx;	
	UV_Vz_x=MP_Vz_xb.*UV_Vz_x + MP_Vz_xc.*(Txz([2:pml+1 nz_pml-pml+2:nz_pml+1],:)-Txz([1:pml nz_pml-pml+1:nz_pml],:))./dz;
	Vz([1:pml nz_pml-pml+1:nz_pml],:)=Vz([1:pml nz_pml-pml+1:nz_pml],:) + MP_Vz([1:pml nz_pml-pml+1:nz_pml],:).*UV_Vz_x;
	UV_Vz_y=MP_Vz_yb.*UV_Vz_y + MP_Vz_yc.*(Tyz(:,[2:pml+1 nx_pml-pml+2:nx_pml+1])-Tyz(:,[1:pml nx_pml-pml+1:nx_pml]))./dx;
	Vz(:,[1:pml nx_pml-pml+1:nx_pml])=Vz(:,[1:pml nx_pml-pml+1:nx_pml]) + MP_Vz(:,[1:pml nx_pml-pml+1:nx_pml]).*UV_Vz_y;
    
    ob_data_Vz(:,:,j)=Vz(recv_x,recv_y);
    Vz_all(j,:,:)=Vz(pml+1:nz_pml-pml,1+pml:nx_pml-pml);
    %if mod(j,50)==0 %��ͼ�Ͷ���
        %figure(3);
        %set(gcf,'outerposition',get(0,'screensize')); %ȫ��
        %imagesc(squeeze(Vz(pml+1:nz_pml-pml,1+pml:nx_pml-pml)));
        %line([1,tunnel_length],[round((nz_pml-tunnel_width)/2-pml),round((nz_pml-tunnel_width)/2-pml)],'LineWidth',2,'Color',[1 1 1]); % ����Ͻ���
        %line([1,tunnel_length],[round((nz_pml+tunnel_width)/2-pml),round((nz_pml+tunnel_width)/2-pml)],'LineWidth',2,'Color',[1 1 1]); % ����½���
        %line([tunnel_length,tunnel_length],[round((nz_pml-tunnel_width)/2-pml),round((nz_pml+tunnel_width)/2-pml)],'LineWidth',2,'Color',[1 1 1]); % ������
       % line([301-npml,301-npml],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % ��ֱ����        
       %line([326,235],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 60����б����
        %line([80,160],[1,nz_pml-2*pml],'LineWidth',2,'Color',[1 1 1]); % 45����б����
%         line([360,141],[35,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 30����б����
%         pos = [237 12 40 40]; rectangle('Position',pos,'Curvature',[1 1]);
        %colorbar;
       % shading interp;
        %grid on;
        %axis equal tight;
       % set(gca,'fontsize',12,'fontweight', 'Bold');        
       % h=suptitle(['Vz at step = ',num2str(j), ',  time = ',num2str(j*dt),' sec']);
       % set(h,'fontsize',18,'fontweight', 'Bold');
       % writeVideo(myavi,getframe(gcf));
  %  end
end
Vz_all_1=Vz_all(:,:);
for(i=1:nx*nz)
  Vz_all_fft(:,i)=fft(Vz_all_1(:,i));
end

end

