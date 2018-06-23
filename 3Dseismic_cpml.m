% Complex frequency shifted convolution PML for FDTD modelling of elastic waves
%P-SV wave propagation in heterogeneous media_velocity-stress finite-difference method
%三维各向同性介质弹性波方程交错网格高阶有限差分法模拟_裴正林
%基于CFS_CPML边界处理的LOVE面波有限差分模拟_孙成禹

clear;
close all;
clc;
%*********************************************************************** 
% fundamental parameters 基本参数
%***********************************************************************
rho_0=1000;           % density=1000 kg/m^3
Vc_0=1000; Vs_0=1000; % Scaling factor, compressional velocity Vp=1000 m/s, shear velocity Vs=1000 m/s
freq=200;		      % (central frequency in Hz) 中心频率Hz

dt=1.0*1.0e-5;
time_window=0.04;
nt=round(time_window/dt);           %计算时间迭代次数
%uniform surrounding rock
param=[ 4     4    4    4    4    %density rho_r
        0     0    0    0    0    %电导率sigma or pi/(Q*eta) 都设为0.
	   3.4   3.4   3.4  3.4  3.4  %Vc_r 
	    2     2    2    2    2];  %Vs_r

%pml parameters
npml=20;
m=4;
kappa_max=1;
alpha_max=0;
bj=0.8;

%***********************************************************************
% mesh parameters 网格参数
%***********************************************************************
   domain_x=80;                %网格的总尺寸in m
   domain_y=80;
   domain_z=80;
   dx=1;
   dy=1;
   dz=1;
   nx=round(domain_x/dx);      %横向的网格数目
   ny=round(domain_y/dx);      %竖向的网格数目
   nz=round(domain_z/dz);      %纵向的网格数目

   fsx=nx/2; fsy=ny/2; fsz=nz/2;	%发射点位置  按nx*ny*nz矩阵定义
%*********************************************************************** 
%激励源 
%*********************************************************************** 
source=zeros(1,nt);   %雷克子波
for i=1:nt
     A=1;t=i*dt-1/freq;
     source(i)=A*(1-2*(pi*freq*t)^2)*exp(-(pi*freq*t)^2);
end

% figure(1);
% plot(source);
% getframe(gca);
% ***********************************************************************
% 介质参数
% ***********************************************************************   
index=ones(nx,ny,nz);
param1=rho_0.*reshape(param(1,index),size(index));
param2=0.*reshape(param(2,index),size(index));
param3=Vc_0.*reshape(param(3,index),size(index));
param4=Vs_0.*reshape(param(4,index),size(index));
param5=param1.*(param3.^2-2.*param4.^2); %lambda=rou*(Vp^2-2*Vs^2)
param6=param1.*param4.^2; %mu=rou*Vs^2
sig_max=max(bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(param3))));
%*********************************************************************** 
%初始化
%*********************************************************************** 
Vx  =zeros(nx+1,ny,nz);   Vy  =zeros(nx,ny+1,nz);   Vz  =zeros(nx,ny,nz+1); 
Txy =zeros(nx+1,ny+1,nz); Tyz =zeros(nx,ny+1,nz+1); Txz =zeros(nx+1,ny,nz+1);
Txx =zeros(nx,ny,nz);     Tyy =zeros(nx,ny,nz);     Tzz =zeros(nx,ny,nz);

% Update Variables in PML
UV_Txy_x =zeros(2*npml-2,ny-1,nz);   UV_Txy_y =zeros(nx-1,2*npml-2,nz);   UV_Txy_z =zeros(nx-1,ny-1,2*npml);
UV_Tyz_x =zeros(2*npml,ny-1,nz-1); UV_Tyz_y =zeros(nx,2*npml-2,nz-1);   UV_Tyz_z =zeros(nx,ny-1,2*npml-2);
UV_Txz_x =zeros(2*npml-2,ny,nz-1);   UV_Txz_y =zeros(nx-1,2*npml,nz-1); UV_Txz_z =zeros(nx-1,ny,2*npml-2);
UV_Txx_x =zeros(2*npml,ny,nz);     UV_Txx_y =zeros(nx,2*npml,nz);     UV_Txx_z =zeros(nx,ny,2*npml);
UV_Tyy_x =zeros(2*npml,ny,nz);     UV_Tyy_y =zeros(nx,2*npml,nz);     UV_Tyy_z =zeros(nx,ny,2*npml);
UV_Tzz_x =zeros(2*npml,ny,nz);     UV_Tzz_y =zeros(nx,2*npml,nz);     UV_Tzz_z =zeros(nx,ny,2*npml);
UV_Vx_x =zeros(2*npml-2,ny,nz);      UV_Vx_y =zeros(nx-1,2*npml,nz);    UV_Vx_z =zeros(nx-1,ny,2*npml);
UV_Vy_x =zeros(2*npml,ny-1,nz);    UV_Vy_y =zeros(nx,2*npml-2,nz);      UV_Vy_z =zeros(nx,ny-1,2*npml);
UV_Vz_x =zeros(2*npml,ny,nz-1);    UV_Vz_y =zeros(nx,2*npml,nz-1);    UV_Vz_z =zeros(nx,ny,2*npml-2);

tic

%*********************************************************************** 
%计算中间矩阵
%***********************************************************************
disp('Computing middle paramters for ...');
disp('Vx');
rou=zeros(nx+1,ny,nz);
rou(2:nx,:,:)=1/2.*(param1(1:nx-1,:,:)+param1(2:nx,:,:));
rou([1 nx+1],:,:)=rou([2 nx],:,:);
MP_Vx=dt./rou;

%PML
sigx=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny,nz]); 
alphax=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny,nz]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny,nz]);
MP_Vx_xb=exp(-((sigx./kappax)+alphax));
MP_Vx_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vx_xb-1);

sigy=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx-1,1,nz]); 
alphay=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx-1,1,nz]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx-1,1,nz]);
MP_Vx_yb=exp(-((sigy./kappay)+alphay));
MP_Vx_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vx_yb-1);

sigz=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx-1,ny,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx-1,ny,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx-1,ny,1]);
MP_Vx_zb=exp(-((sigz./kappaz)+alphaz));
MP_Vx_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Vx_zb-1);
clear rou sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

disp('Vy');
rou=zeros(nx,ny+1,nz);
rou(:,2:ny,:)=1/2.*(param1(:,1:ny-1,:)+param1(:,2:ny,:));
rou(:,[1 ny+1],:)=rou(:,[2 ny],:);
MP_Vy=dt./rou;
%PML
sigx=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny-1,nz]); 
alphax=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny-1,nz]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny-1,nz]);
MP_Vy_xb=exp(-((sigx./kappax)+alphax));
MP_Vy_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vy_xb-1);

sigy=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx,1,nz]); 
alphay=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx,1,nz]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx,1,nz]);
MP_Vy_yb=exp(-((sigy./kappay)+alphay));
MP_Vy_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vy_yb-1);

sigz=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny-1,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny-1,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny-1,1]);
MP_Vy_zb=exp(-((sigz./kappaz)+alphaz));
MP_Vy_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Vy_zb-1);
clear rou sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

disp('Vz');
rou=zeros(nx,ny,nz+1);
rou(:,:,2:nz)=1/2.*(param1(:,:,1:nz-1)+param1(:,:,2:nz));
rou(:,:,[1 nz+1])=rou(:,:,[2 nz]);
MP_Vz=dt./rou;
%PML
sigx=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz-1]); 
alphax=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz-1]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz-1]);
MP_Vz_xb=exp(-((sigx./kappax)+alphax));
MP_Vz_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vz_xb-1);

sigy=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz-1]); 
alphay=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz-1]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz-1]);
MP_Vz_yb=exp(-((sigy./kappay)+alphay));
MP_Vz_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Vz_yb-1);

sigz=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx,ny,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx,ny,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx,ny,1]);
MP_Vz_zb=exp(-((sigz./kappaz)+alphaz));
MP_Vz_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Vz_zb-1);
clear rou sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

disp('Txy');
mu=zeros(nx+1,ny+1,nz);
mu(2:nx,2:ny,:)=1/4.*(param6(1:nx-1,1:ny-1,:)+param6(1:nx-1,2:ny,:)+param6(2:nx,1:ny-1,:)+param6(2:nx,2:ny,:));
mu([1 nx+1],2:ny,:)=mu([2 nx],2:ny,:); mu(2:nx,[1 ny+1],:)=mu(2:nx,[2 ny],:);
MP_Txy=dt.*mu;

sigx=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny-1,nz]); 
alphax=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny-1,nz]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny-1,nz]);
MP_Txy_xb=exp(-((sigx./kappax)+alphax));
MP_Txy_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Txy_xb-1);

sigy=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx-1,1,nz]); 
alphay=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx-1,1,nz]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx-1,1,nz]);
MP_Txy_yb=exp(-((sigy./kappay)+alphay));
MP_Txy_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Txy_yb-1);

sigz=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx-1,ny-1,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx-1,ny-1,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx-1,ny-1,1]);
MP_Txy_zb=exp(-((sigz./kappaz)+alphaz));
MP_Txy_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Txy_zb-1);
clear mu sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

disp('Tyz'); %nx,ny+1,nz+1
mu=zeros(nx,ny+1,nz+1);
mu(:,2:ny,2:nz)=1/4.*(param6(:,1:ny-1,1:nz-1)+param6(:,2:ny,1:nz-1)+param6(:,1:ny-1,2:nz)+param6(:,2:ny,2:nz));
mu(:,[1 ny+1],2:nz)=mu(:,[2 ny],2:nz); mu(:,2:ny,[1 nz+1])=mu(:,2:ny,[2 nz]);
MP_Tyz=dt.*mu;

sigx=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny-1,nz-1]); 
alphax=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny-1,nz-1]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny-1,nz-1]);
MP_Tyz_xb=exp(-((sigx./kappax)+alphax));
MP_Tyz_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Tyz_xb-1);

sigy=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx,1,nz-1]); 
alphay=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx,1,nz-1]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,2*npml-2,1]),[nx,1,nz-1]);
MP_Tyz_yb=exp(-((sigy./kappay)+alphay));
MP_Tyz_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Tyz_yb-1);

sigz=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx,ny-1,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx,ny-1,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx,ny-1,1]);
MP_Tyz_zb=exp(-((sigz./kappaz)+alphaz));
MP_Tyz_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Tyz_zb-1);
clear mu sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

disp('Txz');%nx+1,ny,nz+1
mu=zeros(nx+1,ny,nz+1);
mu(2:nx,:,2:nz)=1/4.*(param6(1:nx-1,:,1:nz-1)+param6(1:nx-1,:,2:nz)+param6(2:nx,:,1:nz-1)+param6(2:nx,:,2:nz));
mu(2:nx,:,[1 nz+1])=mu(2:nx,:,[2 nz]); mu([1 nx+1],:,2:nz)=mu([2 nx],:,2:nz);
MP_Txz=dt.*mu;

%CPML
sigx=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny,nz-1]); 
alphax=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny,nz-1]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[2*npml-2,1,1]),[1,ny,nz-1]);
MP_Txz_xb=exp(-((sigx./kappax)+alphax));
MP_Txz_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Txz_xb-1);

sigy=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx-1,1,nz-1]); 
alphay=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx-1,1,nz-1]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx-1,1,nz-1]);
MP_Txz_yb=exp(-((sigy./kappay)+alphay));
MP_Txz_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Txz_yb-1);

sigz=repmat(sig_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx-1,ny,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx-1,ny,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-1:-1:1 1:npml-1]./npml).^m),[1,1,2*npml-2]),[nx-1,ny,1]);
MP_Txz_zb=exp(-((sigz./kappaz)+alphaz));
MP_Txz_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Txz_zb-1);
clear mu sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

disp('Txx');
lambda=param5; 
mu=param6;
MP_Txx1=dt.*(lambda+2.*mu);
MP_Txx2=dt.*lambda;

sigx=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]); 
alphax=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]);
MP_Txx_xb=exp(-((sigx./kappax)+alphax));
MP_Txx_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Txx_xb-1);

sigy=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]); 
alphay=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]);
MP_Txx_yb=exp(-((sigy./kappay)+alphay));
MP_Txx_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Txx_yb-1);

sigz=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]);
MP_Txx_zb=exp(-((sigz./kappaz)+alphaz));
MP_Txx_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Txx_zb-1);
clear sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

disp('Tyy');
MP_Tyy1=dt.*(lambda+2.*mu);
MP_Tyy2=dt.*lambda;

sigx=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]); 
alphax=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]);
MP_Tyy_xb=exp(-((sigx./kappax)+alphax));
MP_Tyy_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Tyy_xb-1);

sigy=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]); 
alphay=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]);
MP_Tyy_yb=exp(-((sigy./kappay)+alphay));
MP_Tyy_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Tyy_yb-1);

sigz=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]);
MP_Tyy_zb=exp(-((sigz./kappaz)+alphaz));
MP_Tyy_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Tyy_zb-1);
clear sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

disp('Tzz');
MP_Tzz1=dt.*(lambda+2.*mu);
MP_Tzz2=dt.*lambda;

sigx=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]); 
alphax=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]);
kappax=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[2*npml,1,1]),[1,ny,nz]);
MP_Tzz_xb=exp(-((sigx./kappax)+alphax));
MP_Tzz_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Tzz_xb-1);

sigy=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]); 
alphay=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]);
kappay=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,2*npml,1]),[nx,1,nz]);
MP_Tzz_yb=exp(-((sigy./kappay)+alphay));
MP_Tzz_ya=sigy./(sigy.*kappay+kappay.^2.*alphay).*(MP_Tzz_yb-1);

sigz=repmat(sig_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]); 
alphaz=repmat(alpha_max.*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]);
kappaz=repmat(1+(kappa_max-1).*reshape((([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m),[1,1,2*npml]),[nx,ny,1]);
MP_Tzz_zb=exp(-((sigz./kappaz)+alphaz));
MP_Tzz_za=sigz./(sigz.*kappaz+kappaz.^2.*alphaz).*(MP_Tzz_zb-1);
clear sigx kappax alphax sigy kappay alphay sigz kappaz alphaz

%*********************************************************************** 
%更新系数
%*********************************************************************** 
myavi=VideoWriter('3D_Elastic_waves.avi','Motion JPEG AVI');
myavi.FrameRate=5;
% myavi.Height=1024;
% myavi.Width=1920;
open(myavi);
for j=1:nt %Time iteration
    disp([num2str(j) '/' num2str(nt)]);
    
    Vx(fsx,fsy,fsz)=Vx(fsx,fsy,fsz) + source(j);
    Vy(fsx,fsy,fsz)=Vy(fsx,fsy,fsz) + source(j);
    Vz(fsx,fsy,fsz)=Vz(fsx,fsy,fsz) + source(j);
    
    Txx=Txx + MP_Txx1.*( Vx(2:nx+1,:,:)-Vx(1:nx,:,:) )./dx ...
        + MP_Txx2.* ( (Vy(:,2:ny+1,:)-Vy(:,1:ny,:))./dy + (Vz(:,:,2:nz+1)-Vz(:,:,1:nz))./dz );
    UV_Txx_x=MP_Txx_xb.*UV_Txx_x+MP_Txx_xa.*(Vx([2:npml+1 nx-npml+2:nx+1],:,:)-Vx([1:npml nx-npml+1:nx],:,:))./dx;
    Txx([1:npml nx-npml+1:nx],:,:)=Txx([1:npml nx-npml+1:nx],:,:) - MP_Txx1([1:npml nx-npml+1:nx],:,:).*UV_Txx_x;
    UV_Txx_y=MP_Txx_yb.*UV_Txx_y+MP_Txx_ya.*(Vy(:,[2:npml+1 ny-npml+2:ny+1],:)-Vy(:,[1:npml ny-npml+1:ny],:))./dy;
    Txx(:,[1:npml ny-npml+1:ny],:)=Txx(:,[1:npml ny-npml+1:ny],:) - MP_Txx2(:,[1:npml ny-npml+1:ny],:).*UV_Txx_y;
    UV_Txx_z=MP_Txx_zb.*UV_Txx_z+MP_Txx_za.*(Vz(:,:,[2:npml+1 nz-npml+2:nz+1])-Vz(:,:,[1:npml nz-npml+1:nz]))./dz;
    Txx(:,:,[1:npml nz-npml+1:nz])=Txx(:,:,[1:npml nz-npml+1:nz]) - MP_Txx2(:,:,[1:npml nz-npml+1:nz]).*UV_Txx_z;
    
	Tyy=Tyy + MP_Tyy1.*( Vy(:,2:ny+1,:)-Vy(:,1:ny,:) )./dy...
        + MP_Tyy2.* ( (Vx(2:nx+1,:,:)-Vx(1:nx,:,:))./dx + (Vz(:,:,2:nz+1)-Vz(:,:,1:nz))./dz );
    UV_Tyy_x=MP_Tyy_xb.*UV_Tyy_x+MP_Tyy_xa.*(Vx([2:npml+1 nx-npml+2:nx+1],:,:)-Vx([1:npml nx-npml+1:nx],:,:))./dx;
    Tyy([1:npml nx-npml+1:nx],:,:)=Tyy([1:npml nx-npml+1:nx],:,:) - MP_Tyy2([1:npml nx-npml+1:nx],:,:).*UV_Tyy_x;
    UV_Tyy_y=MP_Tyy_yb.*UV_Tyy_y+MP_Tyy_ya.*(Vy(:,[2:npml+1 ny-npml+2:ny+1],:)-Vy(:,[1:npml ny-npml+1:ny],:))./dy;
    Tyy(:,[1:npml ny-npml+1:ny],:)=Tyy(:,[1:npml ny-npml+1:ny],:) - MP_Tyy1(:,[1:npml ny-npml+1:ny],:).*UV_Tyy_y;
    UV_Tyy_z=MP_Tyy_zb.*UV_Tyy_z+MP_Tyy_za.*(Vz(:,:,[2:npml+1 nz-npml+2:nz+1])-Vz(:,:,[1:npml nz-npml+1:nz]))./dz;
    Tyy(:,:,[1:npml nz-npml+1:nz])=Tyy(:,:,[1:npml nz-npml+1:nz]) - MP_Tyy2(:,:,[1:npml nz-npml+1:nz]).*UV_Tyy_z;
    
	Tzz=Tzz + MP_Tzz1.*( Vz(:,:,2:nz+1)-Vz(:,:,1:nz) )./dz...
        + MP_Tzz2.* ( (Vx(2:nx+1,:,:)-Vx(1:nx,:,:))./dx + (Vy(:,2:ny+1,:)-Vy(:,1:ny,:))./dy );
    UV_Tzz_x=MP_Tzz_xb.*UV_Tzz_x+MP_Tzz_xa.*(Vx([2:npml+1 nx-npml+2:nx+1],:,:)-Vx([1:npml nx-npml+1:nx],:,:))./dx;
    Tzz([1:npml nx-npml+1:nx],:,:)=Tzz([1:npml nx-npml+1:nx],:,:) - MP_Tzz2([1:npml nx-npml+1:nx],:,:).*UV_Tzz_x;
    UV_Tzz_y=MP_Tzz_yb.*UV_Tzz_y+MP_Tzz_ya.*(Vy(:,[2:npml+1 ny-npml+2:ny+1],:)-Vy(:,[1:npml ny-npml+1:ny],:))./dy;
    Tzz(:,[1:npml ny-npml+1:ny],:)=Tzz(:,[1:npml ny-npml+1:ny],:) - MP_Tzz2(:,[1:npml ny-npml+1:ny],:).*UV_Tzz_y;
    UV_Tzz_z=MP_Tzz_zb.*UV_Tzz_z+MP_Tzz_za.*(Vz(:,:,[2:npml+1 nz-npml+2:nz+1])-Vz(:,:,[1:npml nz-npml+1:nz]))./dz;
    Tzz(:,:,[1:npml nz-npml+1:nz])=Tzz(:,:,[1:npml nz-npml+1:nz]) - MP_Tzz1(:,:,[1:npml nz-npml+1:nz]).*UV_Tzz_z;
	
    Txy(2:nx,2:ny,:)=Txy(2:nx,2:ny,:)+ MP_Txy(2:nx,2:ny,:).*( (Vx(2:nx,2:ny,:)-Vx(2:nx,1:ny-1,:))./dy + (Vy(2:nx,2:ny,:)-Vy(1:nx-1,2:ny,:))./dx );
    UV_Txy_x=MP_Txy_xb.*UV_Txy_x+MP_Txy_xa.*(Vy([2:npml nx-npml+2:nx],2:ny,:)-Vy([1:npml-1 nx-npml+1:nx-1],2:ny,:))./dx;
    Txy([2:npml nx-npml+2:nx],2:ny,:)=Txy([2:npml nx-npml+2:nx],2:ny,:) - MP_Txy([2:npml nx-npml+2:nx],2:ny,:).*UV_Txy_x;
    UV_Txy_y=MP_Txy_yb.*UV_Txy_y+MP_Txy_ya.*(Vx(2:nx,[2:npml ny-npml+2:ny],:)-Vx(2:nx,[1:npml-1 ny-npml+1:ny-1],:))./dy;
    Txy(2:nx,[2:npml ny-npml+2:ny],:)=Txy(2:nx,[2:npml ny-npml+2:ny],:) - MP_Txy(2:nx,[2:npml ny-npml+2:ny],:).*UV_Txy_y;
    Txy([1 nx+1],:,:)=Txy([2 nx],:,:); Txy(:,[1 ny+1],:)=Txy(:,[2 ny],:);

    Tyz(:,2:ny,2:nz)=Tyz(:,2:ny,2:nz)+ MP_Tyz(:,2:ny,2:nz).*( (Vy(:,2:ny,2:nz)-Vy(:,2:ny,1:nz-1))./dz + (Vz(:,2:ny,2:nz)-Vz(:,1:ny-1,2:nz))./dy );
    UV_Tyz_y=MP_Tyz_yb.*UV_Tyz_y+MP_Tyz_ya.*(Vz(:,[2:npml ny-npml+2:ny],2:nz)-Vz(:,[1:npml-1 ny-npml+1:ny-1],2:nz))./dy;
    Tyz(:,[2:npml ny-npml+2:ny],2:nz)=Tyz(:,[2:npml ny-npml+2:ny],2:nz) - MP_Tyz(:,[2:npml ny-npml+2:ny],2:nz).*UV_Tyz_y;    
    UV_Tyz_z=MP_Tyz_zb.*UV_Tyz_z+MP_Tyz_za.*(Vy(:,2:ny,[2:npml nz-npml+2:nz])-Vy(:,2:ny,[1:npml-1 nz-npml+1:nz-1]))./dz;
    Tyz(:,2:ny,[2:npml nz-npml+2:nz])=Tyz(:,2:ny,[2:npml nz-npml+2:nz]) - MP_Tyz(:,2:ny,[2:npml nz-npml+2:nz]).*UV_Tyz_z;
    Tyz(:,:,[1 nz+1])=Tyz(:,:,[2 nz]); Tyz(:,[1 ny+1],:)=Tyz(:,[2 ny],:);
    
	Txz(2:nx,:,2:nz)=Txz(2:nx,:,2:nz)+ MP_Txz(2:nx,:,2:nz).*( (Vx(2:nx,:,2:nz)-Vx(2:nx,:,1:nz-1))./dz + (Vz(2:nx,:,2:nz)-Vz(1:nx-1,:,2:nz))./dx );
    UV_Txz_x=MP_Txz_xb.*UV_Txz_x+MP_Txz_xa.*(Vz([2:npml nx-npml+2:nx],:,2:nz)-Vz([1:npml-1 nx-npml+1:nx-1],:,2:nz))./dx;
    Txz([2:npml nx-npml+2:nx],:,2:nz)=Txz([2:npml nx-npml+2:nx],:,2:nz) - MP_Txz([2:npml nx-npml+2:nx],:,2:nz).*UV_Txz_x;
    UV_Txz_z=MP_Txz_zb.*UV_Txz_z+MP_Txz_za.*(Vx(2:nx,:,[2:npml nz-npml+2:nz])-Vx(2:nx,:,[1:npml-1 nz-npml+1:nz-1]))./dz;
    Txz(2:nx,:,[2:npml nz-npml+2:nz])=Txz(2:nx,:,[2:npml nz-npml+2:nz]) - MP_Txz(2:nx,:,[2:npml nz-npml+2:nz]).*UV_Txz_z;
    Txz([1 nx+1],:,:)=Txz([2 nx],:,:); Txz(:,:,[1 nz+1])=Txz(:,:,[2 nz]);    
        
    Vx(2:nx,:,:)=Vx(2:nx,:,:)+ MP_Vx(2:nx,:,:).* ( (Txx(2:nx,:,:)-Txx(1:nx-1,:,:))./dx ...
        + (Txy(2:nx,2:ny+1,:)-Txy(2:nx,1:ny,:))./dy + (Txz(2:nx,:,2:nz+1)-Txz(2:nx,:,1:nz))./dz );	
    UV_Vx_x=MP_Vx_xb.*UV_Vx_x+MP_Vx_xa.*(Txx([2:npml ny-npml+2:ny],:,:)-Txx([1:npml-1 nz-npml+1:nz-1],:,:))./dx;
    Vx([2:npml nx-npml+2:nx],:,:)=Vx([2:npml nx-npml+2:nx],:,:) - MP_Vx([2:npml nx-npml+2:nx],:,:).*UV_Vx_x;    
    UV_Vx_y=MP_Vx_yb.*UV_Vx_y+MP_Vx_ya.*(Txy(2:nx,[2:npml+1 ny-npml+2:ny+1],:)-Txy(2:nx,[1:npml ny-npml+1:ny],:))./dy;
    Vx(2:nx,[1:npml ny-npml+1:ny],:)=Vx(2:nx,[1:npml ny-npml+1:ny],:) - MP_Vx(2:nx,[1:npml ny-npml+1:ny],:).*UV_Vx_y;
    UV_Vx_z=MP_Vx_zb.*UV_Vx_z+MP_Vx_za.*(Txz(2:nx,:,[2:npml+1 nz-npml+2:nz+1])-Txz(2:nx,:,[1:npml nz-npml+1:nz]))./dz;
    Vx(2:nx,:,[1:npml nz-npml+1:nz])=Vx(2:nx,:,[1:npml nz-npml+1:nz]) - MP_Vx(2:nx,:,[1:npml nz-npml+1:nz]).*UV_Vx_z;
    Vx([1 nx+1],:,:)=Vx([2 nx],:,:);
    
    Vy(:,2:ny,:)=Vy(:,2:ny,:)+MP_Vy(:,2:ny,:).* ( (Txy(2:nx+1,2:ny,:)-Txy(1:nx,2:ny,:))./dx...
        + (Tyy(:,2:ny,:)-Tyy(:,1:ny-1,:))./dy + (Tyz(:,2:ny,2:nz+1)-Tyz(:,2:ny,1:nz))./dz );
    UV_Vy_x=MP_Vy_xb.*UV_Vy_x+MP_Vy_xa.*(Txy([2:npml+1 nx-npml+2:nx+1],2:ny,:)-Txy([1:npml nx-npml+1:nx],2:ny,:))./dx;
    Vy([1:npml nx-npml+1:nx],2:ny,:)=Vy([1:npml nx-npml+1:nx],2:ny,:) - MP_Vy([1:npml nx-npml+1:nx],2:ny,:).*UV_Vy_x;    
    UV_Vy_y=MP_Vy_yb.*UV_Vy_y+MP_Vy_ya.*(Tyy(:,[2:npml ny-npml+2:ny],:)-Tyy(:,[1:npml-1 nz-npml+1:nz-1],:))./dy;
    Vy(:,[2:npml ny-npml+2:ny],:)=Vy(:,[2:npml ny-npml+2:ny],:) - MP_Vy(:,[2:npml ny-npml+2:ny],:).*UV_Vy_y;    
    UV_Vy_z=MP_Vy_zb.*UV_Vy_z+MP_Vy_za.*(Tyz(:,2:ny,[2:npml+1 nz-npml+2:nz+1])-Tyz(:,2:ny,[1:npml nz-npml+1:nz]))./dz;
    Vy(:,2:ny,[1:npml nz-npml+1:nz])=Vy(:,2:ny,[1:npml nz-npml+1:nz]) - MP_Vy(:,2:ny,[1:npml nz-npml+1:nz]).*UV_Vy_z;
    Vy(:,[1 ny+1],:)=Vy(:,[2 ny],:);
    
    Vz(:,:,2:nz)=Vz(:,:,2:nz)+MP_Vz(:,:,2:nz).* ( (Txz(2:nx+1,:,2:nz)-Txz(1:nx,:,2:nz))./dx...
        + (Tyz(:,2:ny+1,2:nz)-Tyz(:,1:ny,2:nz))./dy + (Tzz(:,:,2:nz)-Tzz(:,:,1:nz-1))./dz );  
    UV_Vz_x=MP_Vz_xb.*UV_Vz_x+MP_Vz_xa.*(Txz([2:npml+1 nx-npml+2:nx+1],:,2:nz)-Txz([1:npml nx-npml+1:nx],:,2:nz))./dx;
    Vz([1:npml nx-npml+1:nx],:,2:nz)=Vz([1:npml nx-npml+1:nx],:,2:nz) - MP_Vz([1:npml nx-npml+1:nx],:,2:nz).*UV_Vz_x;    
    UV_Vz_y=MP_Vz_yb.*UV_Vz_y+MP_Vz_ya.*(Tyz(:,[2:npml+1 ny-npml+2:ny+1],2:nz)-Tyz(:,[1:npml ny-npml+1:ny],2:nz))./dy;
    Vz(:,[1:npml ny-npml+1:ny],2:nz)=Vz(:,[1:npml ny-npml+1:ny],2:nz) - MP_Vz(:,[1:npml ny-npml+1:ny],2:nz).*UV_Vz_y;
    UV_Vz_z=MP_Vz_zb.*UV_Vz_z+MP_Vz_za.*(Tzz(:,:,[2:npml nz-npml+2:nz])-Tzz(:,:,[2:npml nz-npml+2:nz]))./dz;
    Vz(:,:,[2:npml nz-npml+2:nz])=Vz(:,:,[2:npml nz-npml+2:nz]) - MP_Vz(:,:,[2:npml nz-npml+2:nz]).*UV_Vz_z;
    Vz(:,:,[1 nz+1])=Vz(:,:,[2 nz]);
        
    if mod(j,100)==0 %出图和动画
        figure(2);
        set(gcf,'outerposition',get(0,'screensize')); %全屏
        subplot(1,2,1);
        imagesc(squeeze(Vx(npml+1:nx+1-npml,1+npml:ny-npml,fsz)));
        colorbar;
        shading interp;
        axis equal tight;
        set(gca,'fontsize',12,'fontweight', 'Bold');
        
        subplot(1,2,2);
        imagesc(squeeze(Vx(fsx,1+npml:ny-npml,1+npml:nz-npml)));
        colorbar('horiz');
        set(gca,'fontsize',12,'fontweight', 'Bold');
        axis equal tight;
        h=suptitle(['Vx at step = ',num2str(j), ',  time = ',num2str(j*dt),' sec']);
        set(h,'fontsize',18,'fontweight', 'Bold');
        writeVideo(myavi,getframe(gcf));
    end
end
close(myavi);
toc

   
   