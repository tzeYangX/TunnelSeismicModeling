
%PML based on the paper 弹性波正演模拟中pml吸收边界条件的改进（秦臻 et al., 2009）
%CFS-CPML based on the paper: a.基于CFS-CPML边界处理的LOVE面波有限差分模拟
%   b.Complex frequency shifted convolution PML for FDTD modelling of elastic waves
clear;
close all;
clc;

%*********************************************************************** 
% fundamental parameters 基本参数
%***********************************************************************

rho_0=1000; % reference density=1000 kg/m^3
Vc_0=1000;  % reference factor, compressional velocity Vp=1000 m/s
Vs_0=1000;  % reference factor, shear velocity Vs=1000 m/s

freq=200;   % central frequency of source in Hz
dgeophone=1;         % 检波器间距

time_window=0.08;    %总计算时窗
dt=1.0*1.0e-5;
nt=round(time_window/dt);           %计算时间迭代次数

%********************************************************************
% mesh parameters 网格参数
%********************************************************************

dx=single(0.5);				%网格宽度
dy=dx;

%********************************************************************
%PML parameters
%********************************************************************

npml=single(20);
m=single(4);
kappa_max=single(1);
alpha_max=single(0);
bj=single(1);

%********************************************************************
% Sources 激励源 
%********************************************************************

source=zeros(1,nt);   %Ricker wavelet雷克子波
for i=1:nt
     A=1;t=i*dt-1/freq;
     source(i)=A*(1-2*(pi*freq*t)^2)*exp(-(pi*freq*t)^2);
end
figure(1);
plot(source);

%********************************************************************
%介质参数
%********************************************************************

a=imread('Model_Tunnel&90interface.bmp');
a=single(a);
a=a(:,:,1)+a(:,:,2).*1e3+a(:,:,3).*1e6;
b=unique(a);
index=single(zeros(size(a)));
for i=single(1):single(length(b))
	index(a==b(i))=single(i);
end
figure(2);
imagesc(index);
colorbar;
colormap('jet');
title('原模型');
getframe(gcf);
[nx_tunnel,ny_tunnel]=find(a<100); nx_tunnel=unique(nx_tunnel); ny_tunnel=unique(ny_tunnel); %隧道空间所占网格
tunnel_width=length(nx_tunnel);tunnel_length=length(ny_tunnel);

param=single(zeros(i,4));
disp('请输入介质性质，格式为[rho_r 阻尼系数 Vc_r Vs_r]：');
for i=single(1):single(length(b))
	p{i}=input([num2str(i) '=']);
	param(i,:)=p{i};
end
clear a b p;
param=single(param);
%tunnel space:[0.1 0 0.34 0.002]
%surrounding rock:[5 0 3.4 2]
%fault:[3 0 2.4 1.4]
%water-bearing cave:[1 0 1.5 0]
%sediment:[2 0 2 1.15]

[nx,ny]=size(index);
nx=single(nx);
ny=single(ny);

param1=single(reshape(param(index,1),size(index)));
% param2=single(reshape(param(index,2),size(index)));
param3=single(reshape(param(index,3),size(index)));
param4=single(reshape(param(index,4),size(index)));

%observation mode 1
fsx=round(nx/2); fsy=tunnel_length+1; %发射点位置 
recv_x=[round(nx/2+tunnel_width/2+1) round(nx/2-tunnel_width/2)]; %接收点位置
recv_y=tunnel_length-10:-dgeophone/dx:tunnel_length-10-60;

%observation mode 2 
% fsx=[round(nx/2-tunnel_width/2)  round(nx/2+tunnel_width/2+1)]; %发射点位置 
% fsy=(tunnel_length-8:2:tunnel_length-4);
% recv_x=[round(nx/2+tunnel_width/2+1) round(nx/2-tunnel_width/2)]; %接收点位置
% recv_y=tunnel_length-20:-dgeophone/dx:tunnel_length-20-60;

tic;

%********************************************************************
% computing middle parameters 计算中间矩阵
%********************************************************************
disp('Computing middle paramters for ...');
disp('Tyz');
rho_r = zeros(nx,ny+1);
Vc_r = zeros(nx,ny+1);
Vs_r = zeros(nx,ny+1);
rho_r(:,1:ny)=param1; rho_r(:,ny+1)=rho_r(:,ny);
Vc_r(:,1:ny)=param3;  Vc_r(:,ny+1)=Vc_r(:,ny);
Vs_r(:,1:ny)=param4;  Vs_r(:,ny+1)=Vs_r(:,ny);
rho=rho_r.*rho_0; Vc=Vc_r.*Vc_0; Vs=Vs_r.*Vs_0;

MP_Tyz=Vs.^2.*rho.*dt;
MP_Tyz_ya=ones(size(rho_r)); 

sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(Vc)));
sigy=sig_max.*repmat(([npml-1:-1:1 1:npml-1]./npml).^m,[nx,1]); %pml层数少1
kappay=1+(kappa_max-1).*([npml-1:-1:1 1:npml-1]./npml).^m; %pml层数少1
kappay=repmat(kappay,[nx,1]);
alphay=alpha_max.*([npml-1:-1:1 1:npml-1]./npml); %=0 %pml层数少1
alphay=repmat(alphay,[nx,1]);
MP_Tyz_ya(:,[2:npml ny-npml+2:ny])=1./kappay;
MP_Tyz_yb=exp(-((sigy./kappay)+alphay));
MP_Tyz_yc=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Tyz_yb-1);
clear rho_r Vs_r Vc_r rho Vs Vc sig_max sigy kappay alphay

disp('Txz');
rho_r=zeros(nx+1,ny);
Vc_r=zeros(nx+1,ny);
Vs_r=zeros(nx+1,ny);
rho_r(1:nx,:)=param1; rho_r(nx+1,:)=rho_r(nx,:);
Vc_r(1:nx,:)=param3;  Vc_r(nx+1,:)=Vc_r(nx,:);
Vs_r(1:nx,:)=param4;  Vs_r(nx+1,:)=Vs_r(nx,:);
rho=rho_r.*rho_0; Vc=Vc_r.*Vc_0; Vs=Vs_r.*Vs_0;

MP_Txz=Vs.^2.*rho.*dt;
MP_Txz_xa=ones(size(rho_r));

sig_max=bj.*(m+1)./(150.*pi.*dy).*sqrt(max(max(Vc)));
sigx=sig_max.*repmat(([npml-1:-1:1 1:npml-1]'./npml).^m,[1,ny]); %pml层数少1
kappax=1+(kappa_max-1).*([npml-1:-1:1 1:npml-1]'./npml).^m; %pml层数少1
kappax=repmat(kappax,[1,ny]);
alphax=alpha_max.*([npml-1:-1:1 1:npml-1]'./npml); %=0 %pml层数少1
alphax=repmat(alphax,[1,ny]);
MP_Txz_xa([2:npml nx-npml+2:nx],:)=1./kappax; 
MP_Txz_xb=exp(-((sigx./kappax)+alphax));
MP_Txz_xc=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Txz_xb-1);
clear rho_r Vs_r Vc_r rho Vs Vc sig_max sigx kappay alphay

disp('Vz');
rho_r=param1; Vc_r=param3;  
rho=rho_r.*rho_0; Vc=Vc_r.*Vc_0;

MP_Vz=dt./rho;
MP_Vz_xa=ones(size(rho_r)); 
MP_Vz_ya=ones(size(rho_r));

sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(Vc)));
sigx=sig_max.*repmat((([npml:-1:1 1:npml]'-0.5)./npml).^m,[1,ny]);
kappax=1+(kappa_max-1).*(([npml:-1:1 1:npml]'-0.5)./npml).^m;
kappax=repmat(kappax,[1,ny]);
alphax=alpha_max.*(([npml:-1:1 1:npml]'-0.5)./npml); %=0
alphax=repmat(alphax,[1,ny]);
MP_Vz_xa([1:npml nx-npml+1:nx],:)=1./kappax;
MP_Vz_xb=exp(-((sigx./kappax)+alphax));
MP_Vz_xc=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vz_xb-1);

sig_max=bj.*(m+1)./(150.*pi.*dy).*sqrt(max(max(Vc)));
sigy=sig_max.*repmat((([npml:-1:1 1:npml]-0.5)./npml).^m,[nx,1]);
kappay=1+(kappa_max-1).*(([npml:-1:1 1:npml]-0.5)./npml).^m;
kappay=repmat(kappay,[nx,1]);
alphay=alpha_max.*(([npml:-1:1 1:npml]-0.5)./npml); %=0
alphay=repmat(alphay,[nx,1]);
MP_Vz_ya(:,[1:npml ny-npml+1:ny])=1./kappay; 
MP_Vz_yb=exp(-((sigy./kappay)+alphay));
MP_Vz_yc=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vz_yb-1);
clear rho_r Vs_r Vc_r rho sig Vs Vc  
clear sig_max sigx sigy alphay kappay
clear param1 param2 param3 param4

%********************************************************************
%初始化
%********************************************************************

Tyz=single(zeros(nx,ny+1));
UV_Tyz_y=single(zeros(nx,2*(npml-1)));
Txz=single(zeros(nx+1,ny));
UV_Txz_x=single(zeros(2*(npml-1),ny));
Vz=single(zeros(nx,ny));
UV_Vz_x=single(zeros(2*npml,ny));
UV_Vz_y=single(zeros(nx,2*npml));

ob_data_Vz=zeros(length(recv_x),length(recv_y),nt);
% frwd_field=zeros(nx-2*npml,ny-2*npml,nt);

%********************************************************************
% computing forward wavefiled
%********************************************************************
myavi=VideoWriter('2D_forward_SH_waves.avi','Motion JPEG AVI');
myavi.FrameRate=5;
open(myavi);

for j=single(1):single(nt)
	
    Vz(fsx,fsy)=Vz(fsx,fsy) + source(j);	
    
	Txz(2:nx,:)=Txz(2:nx,:) + MP_Txz(2:nx,:).*MP_Txz_xa(2:nx,:).*(Vz(2:nx,:)-Vz(1:nx-1,:))./dx;
	UV_Txz_x=MP_Txz_xb.*UV_Txz_x + MP_Txz_xc.*(Vz([2:npml nx-npml+2:nx],:)-Vz([1:npml-1 nx-npml+1:nx-1],:))./dx;
	Txz([2:npml nx-npml+2:nx],:)=Txz([2:npml nx-npml+2:nx],:) + MP_Txz([2:npml nx-npml+2:nx],:).*UV_Txz_x;
	Txz([1 nx+1],:)=Txz([2 nx],:);
		
	Tyz(:,2:ny)=Tyz(:,2:ny) + MP_Tyz(:,2:ny).*MP_Tyz_ya(:,2:ny).*(Vz(:,2:ny)-Vz(:,1:ny-1))./dy;
	UV_Tyz_y=MP_Tyz_yb.*UV_Tyz_y + MP_Tyz_yc.*(Vz(:,[2:npml ny-npml+2:ny])-Vz(:,[1:npml-1 ny-npml+1:ny-1]))./dy;
	Tyz(:,[2:npml ny-npml+2:ny])=Tyz(:,[2:npml ny-npml+2:ny]) + MP_Tyz(:,[2:npml ny-npml+2:ny]).*UV_Tyz_y;
    Tyz(:,[1 ny+1])=Tyz(:,[2 ny]);
    
	Vz=Vz + MP_Vz.*MP_Vz_xa.*(Txz(2:nx+1,:)-Txz(1:nx,:))./dx + MP_Vz.*MP_Vz_ya.*(Tyz(:,2:ny+1)-Tyz(:,1:ny))./dy;	
	UV_Vz_x=MP_Vz_xb.*UV_Vz_x + MP_Vz_xc.*(Txz([2:npml+1 nx-npml+2:nx+1],:)-Txz([1:npml nx-npml+1:nx],:))./dx;
	Vz([1:npml nx-npml+1:nx],:)=Vz([1:npml nx-npml+1:nx],:) + MP_Vz([1:npml nx-npml+1:nx],:).*UV_Vz_x;
	UV_Vz_y=MP_Vz_yb.*UV_Vz_y + MP_Vz_yc.*(Tyz(:,[2:npml+1 ny-npml+2:ny+1])-Tyz(:,[1:npml ny-npml+1:ny]))./dy;
	Vz(:,[1:npml ny-npml+1:ny])=Vz(:,[1:npml ny-npml+1:ny]) + MP_Vz(:,[1:npml ny-npml+1:ny]).*UV_Vz_y;
    
%     ob_data_Vz(:,:,j)=Vz(recv_x,recv_y);
%     frwd_field(:,:,j)=Vz(npml+1:nx-npml,1+npml:ny-npml);
    
	if mod(j,50)==0 %出图和动画
        figure(3);
        set(gcf,'outerposition',get(0,'screensize')); %全屏
        imagesc(squeeze(Vz(npml+1:nx-npml,1+npml:ny-npml)));
        line([1,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx-tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 隧道上界面
        line([1,tunnel_length-npml],[round((nx+tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 隧道下界面
        line([tunnel_length-npml,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 掌子面
        line([301-npml,301-npml],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 垂直界面        
%         line([326,235],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 60度倾斜界面
%         line([360,141],[35,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 30度倾斜界面
%         pos = [237 12 40 40]; rectangle('Position',pos,'Curvature',[1 1]);
        colorbar;
        shading interp;
        grid on;
        axis equal tight;
        set(gca,'fontsize',12,'fontweight', 'Bold');        
        h=suptitle(['Vz at step = ',num2str(j), ',  time = ',num2str(j*dt),' sec']);
        set(h,'fontsize',18,'fontweight', 'Bold');
        writeVideo(myavi,getframe(gcf));
	end
end
close(myavi);
toc;
save ob_data.mat ob_data_Vz
close all