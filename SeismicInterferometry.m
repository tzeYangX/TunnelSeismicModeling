% Complex frequency shifted convolution PML for FDTD modelling of elastic waves

clear;
close all;
clc;

time_window=50;%�ܼ���ʱ��
dt=5e-5;
% dt=single(1/(sqrt(1/dx^2+1/dy^2)*4000))/3;%����ʱ�䲽��
nt=round(time_window/dt); %����ʱ���������

dx=1;%������
dy=dx;
%---------------------
%����ͼ�����ݣ�ȷ��ģ�ͽṹ
a=imread('Model_Tunnel&fault.bmp');
a=double(a);
a=a(:,:,1)+a(:,:,2).*1e3+a(:,:,3).*1e6;
b=unique(a);
index=zeros(size(a));
for i=1:length(b)
	index(a==b(i))=i;
end
[nx,ny]=size(index);


file=dir('PSV_*.m');  % ����ǰĿ¼�µ��ļ�������Ϣ�浽file����
name_ind=strfind(file.name,'m.m');
name_tunnelface=str2num([file.name([name_ind-5,name_ind-4,name_ind-3]),'.',file.name([name_ind-1])]);
excavation_speed=1; % ģ��TBMÿ�ο��ڣ�������ÿ��ǰ��1����
excavation_num=(name_tunnelface-100)/(excavation_speed*dx)+1; % ���μ�����ģ��ڼ��ο��ڣ���������100mʱΪ���
tunnel_length=100/dx+(excavation_num-1)*excavation_speed; %tunnel length ������ 100.0+(excavation_num-1)*0.5 m
tunnel_width=12/dx; % tunnel thickness 12 m

%********************************************************************
%���ʲ���
%********************************************************************
param=zeros(length(b),4);
rho_0=1000; % reference density=1000 kg/m^3
Vp_0=1000;  % reference factor, compressional velocity Vp=1000 m/s
Vs_0=1000;  % reference factor, shear velocity Vs=1000 m/s
% disp('�������������(��λ10^3)����ʽΪ[rho_r ����ϵ�� Vp_r Vs_r]��');
% for i=single(1):single(length(b))
	% p{i}=input([num2str(i) '=']);
	% param(i,:)=p{i};
% end
% clear a b p;

%���������쳣�����Բ���
param(:)=[5 0 4 2.388;
          4 0 3 1.594];

param1=rho_0.*reshape(param(index,1),size(index));
% param2=reshape(param(index,2),size(index));
param3=Vp_0.*reshape(param(index,3),size(index));
param4=Vs_0.*reshape(param(index,4),size(index));
% ǿ�и�������򸳿���
param1((nx-tunnel_width)/2:(nx+tunnel_width)/2,1:tunnel_length)=rho_0*0.1; 
% param2((nx-tunnel_width)/2:(nx+tunnel_width)/2,1:tunnel_length)=0; 
param3((nx-tunnel_width)/2:(nx+tunnel_width)/2,1:tunnel_length)=Vp_0*0.34; 
param4((nx-tunnel_width)/2:(nx+tunnel_width)/2,1:tunnel_length)=Vs_0*0.002; 


%********************************************************************
% Sources and observation mode ����Դ�͹۲ⷽʽ
%********************************************************************
load('source.mat');
source=double(source);

dgeophone=1;         % �첨�����,m
dsource=1;   % ��Դ���,m
% %---------------------
fsx=round(nx/2+tunnel_width/2):-dsource/dx:round(nx/2-tunnel_width/2); % ÿ���������һ����Դ
fsy=tunnel_length+1; % ������ǰ��1����
% %---------------------
ntrac=120; % ���մ���������
irX0=90/dx;     % ��һ���첨���̶���90m��λ�ô�
recv_x=[round(nx/2+tunnel_width/2)+1 round(nx/2-tunnel_width/2)-1];
recv_y=irX0:-dgeophone/dx:irX0-(ntrac/2-1)*dgeophone/dx;
pilot_x=round(nx/2); % �ȵ�������
pilot_y=tunnel_length+1;


%********************************************************************
%PML parameters
%********************************************************************
npml=5;
m=4;
kappa_max=1;
alpha_max=0;
bj=0.8;

%********************************************************************
%�����м����
%********************************************************************
tic;

%��÷ϵ��
param5=param1.*(param3.^2-2.*param4.^2); %lambda=rou*(Vp^2-2*Vs^2)
param6=param1.*param4.^2; %mu=rou*Vs^2

%���˥������
%���ף�Convolutional PML (CPML): an efficient FDTD implementation of the CFS-PML for arbitrary media
sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(param3)));

disp('Computing middle paramters for ...');
disp('Vx'); %Vx nx+1,ny
rou=zeros(nx+1,ny);
rou(2:nx,:)=1/2.*(param1(1:nx-1,:)+param1(2:nx,:));
rou([1 nx+1],:)=rou([2 nx],:);%�������ܶ���ȣ�����PML��
MP_Vx=dt./rou;


%����Bj��Cj
sigx=sig_max.*repmat(([npml-1:-1:1,1:npml-1]./npml)'.^m,[1,ny]); %pml������1
alphax=repmat(alpha_max.*(([npml-1:-1:1,1:npml-1]./npml)'.^m),[1,ny]);
kappax=repmat(1+(kappa_max-1).*(([npml-1:-1:1,1:npml-1]./npml)'.^m),[1,ny]);
MP_Vx_xb=exp(-((sigx./kappax)+alphax));%Bj
MP_Vx_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vx_xb-1);%Cj

sigy=sig_max.*repmat(([npml-0.5:-1:0.5,0.5:npml-0.5]./npml).^m,[nx-1,1]);
kappay=1+(kappa_max-1).*repmat(([npml-0.5:-1:0.5,0.5:npml-0.5]./npml).^m,[nx-1,1]);
alphay=alpha_max.*repmat(([npml-0.5:-1:0.5,0.5:npml-0.5]./npml),[nx-1,1]);
MP_Vx_yb=exp(-((sigy./kappay)+alphay));
MP_Vx_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vx_yb-1);

clear rou sigx kappax alphax sigy kappay alphay

disp('Vy'); %Vy nx,ny+1
rou=single(zeros(nx,ny+1));
rou(:,2:ny)=1/2.*(param1(:,1:ny-1)+param1(:,2:ny));
rou(:,[1 ny+1])=rou(:,[2 ny]);
MP_Vy=dt./rou;

sigx=sig_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml).^m,[1,ny-1]);
kappax=1+(kappa_max-1).*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml).^m,[1,ny-1]);
alphax=alpha_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml),[1,ny-1]);
MP_Vy_xb=exp(-((sigx./kappax)+alphax));
MP_Vy_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vy_xb-1);

sigy=sig_max.*repmat(([npml-1:-1:1 1:npml-1]./npml).^m,[nx,1]); 
kappay=1+(kappa_max-1).*repmat(([npml-1:-1:1 1:npml-1]./npml).^m,[nx,1]);
alphay=alpha_max.*repmat(([npml-1:-1:1 1:npml-1]./npml),[nx,1]);
MP_Vy_yb=exp(-((sigy./kappay)+alphay));
MP_Vy_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vy_yb-1);
clear rou sigx kappax alphax sigy kappay alphay

disp('Txx'); %Txx nx,ny
lambda=param5; 
mu=param6;
MP_Txx1=dt.*(lambda+2.*mu);
MP_Txx2=dt.*lambda;

sigx=sig_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml).^m,[1,ny]);
kappax=1+(kappa_max-1).*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml).^m,[1,ny]);
alphax=alpha_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml),[1,ny]);
MP_Txx_xb=exp(-((sigx./kappax)+alphax));
MP_Txx_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Txx_xb-1);

sigy=sig_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m,[nx,1]);
kappay=1+(kappa_max-1).*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m,[nx,1]);
alphay=alpha_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml),[nx,1]);
MP_Txx_yb=exp(-((sigy./kappay)+alphay));
MP_Txx_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Txx_yb-1);
clear sigx kappax alphax sigy kappay alphay

disp('Tyy'); %Tyy nx,ny
MP_Tyy1=dt.*(lambda+2.*mu);
MP_Tyy2=dt.*lambda;

sigx=sig_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml).^m,[1,ny]);
kappax=1+(kappa_max-1).*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml).^m,[1,ny]);
alphax=alpha_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]'./npml),[1,ny]);
MP_Tyy_xb=exp(-((sigx./kappax)+alphax));
MP_Tyy_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Tyy_xb-1);

sigy=sig_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m,[nx,1]);
kappay=1+(kappa_max-1).*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m,[nx,1]);
alphay=alpha_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml),[nx,1]);
MP_Tyy_yb=exp(-((sigy./kappay)+alphay));
MP_Tyy_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Tyy_yb-1);
clear lambda mu sigx kappax alphax sigy kappay alphay

disp('Txy'); %Txy nx+1,ny+1
mu=single(zeros(nx+1,ny+1));
mu(2:nx,2:ny)=1/4.*(param6(1:nx-1,1:ny-1)+param6(1:nx-1,2:ny)+param6(2:nx,1:ny-1)+param6(2:nx,2:ny));
mu([1 nx+1],2:ny)=mu([2 nx],2:ny); mu(:,[1 ny+1])=mu(:,[2 ny]);
MP_Txy=dt.*mu;

sigx=sig_max.*repmat(([npml-1:-1:1 1:npml-1]./npml)'.^m,[1,ny-1]); %pml������1
alphax=repmat(alpha_max.*(([npml-1:-1:1 1:npml-1]./npml)'.^m),[1,ny-1]);
kappax=repmat(1+(kappa_max-1).*(([npml-1:-1:1 1:npml-1]./npml)'.^m),[1,ny-1]);
MP_Txy_xb=exp(-((sigx./kappax)+alphax));
MP_Txy_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Txy_xb-1);

sigy=sig_max.*repmat(([npml-1:-1:1 1:npml-1]./npml).^m,[nx-1,1]); 
kappay=1+(kappa_max-1).*repmat(([npml-1:-1:1 1:npml-1]./npml).^m,[nx-1,1]);
alphay=alpha_max.*repmat(([npml-1:-1:1 1:npml-1]./npml),[nx-1,1]);
MP_Txy_yb=exp(-((sigy./kappay)+alphay));
MP_Txy_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Txy_yb-1);
clear mu sigx kappax alphax sigy kappay alphay
% clear sig_max param1 param2 param3 param4 param5 param6



%********************************************************************
%��ʼ������
Vx =zeros(nx+1,ny);
Vy =zeros(nx,ny+1);
Txy =zeros(nx+1,ny+1); 
Txx =zeros(nx,ny);
Tyy =zeros(nx,ny); 


% ��ʼ��PML�߽�
UV_Txy_x =zeros(2*npml-2,ny-1);
UV_Txy_y =zeros(nx-1,2*npml-2);
UV_Txx_x =zeros(2*npml,ny);
UV_Txx_y =zeros(nx,2*npml);
UV_Tyy_x =zeros(2*npml,ny);
UV_Tyy_y =zeros(nx,2*npml);
UV_Vx_x =zeros(2*npml-2,ny);
UV_Vx_y =zeros(nx-1,2*npml);
UV_Vy_x =zeros(2*npml,ny-1);
UV_Vy_y =zeros(nx,2*npml-2);

ob_data_Vx=zeros(nt,length(recv_y),length(recv_x));
ob_data_Vx_pilot=zeros(nt,length(pilot_y),length(pilot_x));
ob_data_Vy=zeros(nt,length(recv_y),length(recv_x));
ob_data_Vy_pilot=zeros(nt,length(pilot_y),length(pilot_x));


tunnelface=sprintf('%.1f',tunnel_length*dx);
tunnel_face=['_',tunnelface(1:3),'_',tunnelface(5)];
disp(['TunnelFace at ',tunnelface,'m  Program is now working.......']);



for j=1:nt
	
	if mod(j,10000)==0
		j
    end
    
    
    if j<=size(source,1)
        Vx(fsx,fsy)=Vx(fsx,fsy) + source(j,:);
        Vy(fsx,fsy)=Vy(fsx,fsy) + source(j,:);
    end	
	
	Vx(2:nx,:)=Vx(2:nx,:)+ MP_Vx(2:nx,:).* ( (Txx(2:nx,:)-Txx(1:nx-1,:))./dx + (Txy(2:nx,2:ny+1)-Txy(2:nx,1:ny))./dy );	
	UV_Vx_x=MP_Vx_xb.*UV_Vx_x+MP_Vx_xa.*(Txx([2:npml nx-npml+2:nx],:)-Txx([1:npml-1 nx-npml+1:nx-1],:))./dx;
	Vx([2:npml nx-npml+2:nx],:)=Vx([2:npml nx-npml+2:nx],:) + MP_Vx([2:npml nx-npml+2:nx],:).*UV_Vx_x; %npml-1��pml,pml�߽��Vxδ˥��
	UV_Vx_y=MP_Vx_yb.*UV_Vx_y+MP_Vx_ya.*(Txy(2:nx,[2:npml+1 ny-npml+2:ny+1])-Txy(2:nx,[1:npml ny-npml+1:ny]))./dy;
	Vx(2:nx,[1:npml ny-npml+1:ny])=Vx(2:nx,[1:npml ny-npml+1:ny]) + MP_Vx(2:nx,[1:npml ny-npml+1:ny]).*UV_Vx_y;
	Vx([1 nx+1],:)=Vx([2 nx],:);
	
	Vy(:,2:ny)=Vy(:,2:ny)+MP_Vy(:,2:ny).* ( (Txy(2:nx+1,2:ny)-Txy(1:nx,2:ny))./dx + (Tyy(:,2:ny)-Tyy(:,1:ny-1))./dy );
	UV_Vy_x=MP_Vy_xb.*UV_Vy_x+MP_Vy_xa.*(Txy([2:npml+1 nx-npml+2:nx+1],2:ny)-Txy([1:npml nx-npml+1:nx],2:ny))./dx;
	Vy([1:npml,nx-npml+1:nx],2:ny)=Vy([1:npml nx-npml+1:nx],2:ny) + MP_Vy([1:npml nx-npml+1:nx],2:ny).*UV_Vy_x;    
	UV_Vy_y=MP_Vy_yb.*UV_Vy_y+MP_Vy_ya.*(Tyy(:,[2:npml ny-npml+2:ny])-Tyy(:,[1:npml-1 ny-npml+1:ny-1]))./dy;
	Vy(:,[2:npml ny-npml+2:ny])=Vy(:,[2:npml,ny-npml+2:ny]) + MP_Vy(:,[2:npml ny-npml+2:ny]).*UV_Vy_y;
	Vy(:,[1 ny+1])=Vy(:,[2 ny]);
	
	Txx=Txx + MP_Txx1.*( Vx(2:nx+1,:)-Vx(1:nx,:) )./dx + MP_Txx2.* ( (Vy(:,2:ny+1)-Vy(:,1:ny))./dy );
	UV_Txx_x=MP_Txx_xb.*UV_Txx_x+MP_Txx_xa.*(Vx([2:npml+1 nx-npml+2:nx+1],:)-Vx([1:npml nx-npml+1:nx],:))./dx;
	Txx([1:npml nx-npml+1:nx],:)=Txx([1:npml nx-npml+1:nx],:)+MP_Txx1([1:npml nx-npml+1:nx],:).*UV_Txx_x;    
	UV_Txx_y=MP_Txx_yb.*UV_Txx_y+MP_Txx_ya.*(Vy(:,[2:npml+1 ny-npml+2:ny+1])-Vy(:,[1:npml ny-npml+1:ny]))./dy;
	Txx(:,[1:npml ny-npml+1:ny])=Txx(:,[1:npml ny-npml+1:ny])+MP_Txx2(:,[1:npml ny-npml+1:ny]).*UV_Txx_y;
	
	Tyy=Tyy + MP_Tyy1.*( Vy(:,2:ny+1)-Vy(:,1:ny) )./dy + MP_Tyy2.* ( (Vx(2:nx+1,:)-Vx(1:nx,:))./dx );
	UV_Tyy_x=MP_Tyy_xb.*UV_Tyy_x+MP_Tyy_xa.*(Vx([2:npml+1 nx-npml+2:nx+1],:)-Vx([1:npml nx-npml+1:nx],:))./dx;
	Tyy([1:npml nx-npml+1:nx],:)=Tyy([1:npml nx-npml+1:nx],:)+MP_Tyy2([1:npml nx-npml+1:nx],:).*UV_Tyy_x;
	UV_Tyy_y=MP_Tyy_yb.*UV_Tyy_y+MP_Tyy_ya.*(Vy(:,[2:npml+1 ny-npml+2:ny+1])-Vy(:,[1:npml ny-npml+1:ny]))./dy;
	Tyy(:,[1:npml ny-npml+1:ny])=Tyy(:,[1:npml ny-npml+1:ny])+MP_Tyy1(:,[1:npml ny-npml+1:ny]).*UV_Tyy_y;
	
	Txy(2:nx,2:ny)=Txy(2:nx,2:ny)+ MP_Txy(2:nx,2:ny).*( (Vx(2:nx,2:ny)-Vx(2:nx,1:ny-1))./dy + (Vy(2:nx,2:ny)-Vy(1:nx-1,2:ny))./dx );
	UV_Txy_x=MP_Txy_xb.*UV_Txy_x+MP_Txy_xa.*(Vy([2:npml nx-npml+2:nx],2:ny)-Vy([1:npml-1 nx-npml+1:nx-1],2:ny))./dx;
	Txy([2:npml nx-npml+2:nx],2:ny)=Txy([2:npml nx-npml+2:nx],2:ny)+MP_Txy([2:npml nx-npml+2:nx],2:ny).*UV_Txy_x;
	UV_Txy_y=MP_Txy_yb.*UV_Txy_y+MP_Txy_ya.*(Vx(2:nx,[2:npml ny-npml+2:ny])-Vx(2:nx,[1:npml-1 ny-npml+1:ny-1]))./dy;
	Txy(2:nx,[2:npml ny-npml+2:ny])=Txy(2:nx,[2:npml ny-npml+2:ny])+MP_Txy(2:nx,[2:npml ny-npml+2:ny]).*UV_Txy_y;
	Txy([1 nx+1],2:ny)=Txy([2 nx],2:ny); Txy(:,[1 ny+1])=Txy(:,[2 ny]);
    
    
    
    	%--------------------
	% ��������¼
	ob_data_Vx(j,:,:)=Vx(recv_x,recv_y)';
	ob_data_Vx_pilot(j,:,:)=(Vx(pilot_x,pilot_y))';
	ob_data_Vy(j,:,:)=Vy(recv_x,recv_y)';
	ob_data_Vy_pilot(j,:,:)=(Vy(pilot_x,pilot_y))';

	


	

end


