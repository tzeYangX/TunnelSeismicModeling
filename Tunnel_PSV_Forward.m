
%CFS-CPML based on the paper: a.基于CFS-CPML边界处理的LOVE面波有限差分模拟
%   b.Complex frequency shifted convolution PML for FDTD modelling of elastic waves
clear;
close all;
clc;

%*********************************************************************** 
% fundamental parameters 基本参数
%***********************************************************************

rho_0=1000; % reference density=1000 kg/m^3
Vp_0=1000;  % reference factor, compressional velocity Vp=1000 m/s
Vs_0=1000;  % reference factor, shear velocity Vs=1000 m/s

freq=single(200);   % central frequency of source in Hz
dgeophone=1;         % 检波器间距

time_window=single(0.08);    %总计算时窗
dt=1.0*1.0e-5;
% dt=single(1/(sqrt(1/dx^2+1/dy^2)*4000))/3;         %计算时间步长
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
bj=single(0.8);

%********************************************************************
% Sources 激励源 
%********************************************************************

source=single(zeros(1,nt));   %Ricker wavelet雷克子波
for i=single(1):single(nt)
     A=1;t=i*dt-1/freq;
     source(i)=A*(1-2*(pi*freq*t)^2)*exp(-(pi*freq*t)^2);
end
figure(1);
plot(source);

%********************************************************************
%介质参数
%********************************************************************

a=imread('Model_Tunnel&fault.bmp');
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
disp('请输入介质性质(单位10^3)，格式为[rho_r 阻尼系数 Vp_r Vs_r]：');
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

param1=rho_0.*single(reshape(param(index,1),size(index)));
% param2=single(reshape(param(index,2),size(index)));
param3=Vp_0.*single(reshape(param(index,3),size(index)));
param4=Vs_0.*single(reshape(param(index,4),size(index)));

fsx=round(nx/2); %发射点位置 
% fsx=[round(nx/2-tunnel_width/2) round(nx/2) round(nx/2+tunnel_width/2)];  
fsy=tunnel_length+1;
recv_x=[round(nx/2+tunnel_width/2+1) round(nx/2-tunnel_width/2)]; %接收点位置
% recv_x=round(nx/2+tunnel_width/2+1);
recv_y=tunnel_length-10:-dgeophone/dx:tunnel_length-10-60;

%********************************************************************
%计算中间矩阵
%********************************************************************
tic;
param5=param1.*(param3.^2-2.*param4.^2); %lambda=rou*(Vp^2-2*Vs^2)
param6=param1.*param4.^2; %mu=rou*Vs^2

sig_max=bj.*(m+1)./(150.*pi.*dx).*sqrt(max(max(param3)));

disp('Computing middle paramters for ...');
disp('Vx'); %Vx nx+1,ny
rou=single(zeros(nx+1,ny));
rou(2:nx,:)=1/2.*(param1(1:nx-1,:)+param1(2:nx,:));
rou([1 nx+1],:)=rou([2 nx],:);
MP_Vx=dt./rou;
%PML
sigx=sig_max.*repmat(([npml-1:-1:1 1:npml-1]./npml)'.^m,[1,ny]); %pml层数少1
alphax=repmat(alpha_max.*(([npml-1:-1:1 1:npml-1]./npml)'.^m),[1,ny]);
kappax=repmat(1+(kappa_max-1).*(([npml-1:-1:1 1:npml-1]./npml)'.^m),[1,ny]);
MP_Vx_xb=exp(-((sigx./kappax)+alphax));
MP_Vx_xa=sigx./(sigx.*kappax+kappax.^2.*alphax).*(MP_Vx_xb-1);

sigy=sig_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m,[nx-1,1]);
kappay=1+(kappa_max-1).*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml).^m,[nx-1,1]);
alphay=alpha_max.*repmat(([npml-0.5:-1:0.5 0.5:npml-0.5]./npml),[nx-1,1]);
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

sigx=sig_max.*repmat(([npml-1:-1:1 1:npml-1]./npml)'.^m,[1,ny-1]); %pml层数少1
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
clear sig_max param1 param2 param3 param4 param5 param6

%********************************************************************
%初始化
%********************************************************************

Vx  =zeros(nx+1,ny);   Vy  =zeros(nx,ny+1);   
Txy =zeros(nx+1,ny+1); 
Txx =zeros(nx,ny);     Tyy =zeros(nx,ny); 
% Update Variables in PML
UV_Txy_x =zeros(2*npml-2,ny-1);   UV_Txy_y =zeros(nx-1,2*npml-2);
UV_Txx_x =zeros(2*npml,ny);       UV_Txx_y =zeros(nx,2*npml);
UV_Tyy_x =zeros(2*npml,ny);       UV_Tyy_y =zeros(nx,2*npml);
UV_Vx_x =zeros(2*npml-2,ny);      UV_Vx_y =zeros(nx-1,2*npml);
UV_Vy_x =zeros(2*npml,ny-1);      UV_Vy_y =zeros(nx,2*npml-2);

ob_data_Vx=zeros(length(recv_x),length(recv_y),nt);
ob_data_Vy=zeros(length(recv_x),length(recv_y),nt);

%********************************************************************
%computing forward wavefield
%********************************************************************
myavi=VideoWriter('2D_forward_PSV_waves&Vy.avi','Motion JPEG AVI');
myavi.FrameRate=5;
open(myavi)

for j=single(1):single(nt)
%     eta=0; %eta 为地震波传播方向与隧道轴线的夹角，逆时针为正
% 	Vx(fsx,fsy)=Vx(fsx,fsy) + sin(eta)*source(j);
%     Vy(fsx,fsy)=Vy(fsx,fsy) + cos(eta)*source(j);
% 	Txy(round(nx/2),round(ny/2))=source(j);
%     Vx(round(nx/2),round(ny/2))=Vx(round(nx/2),round(ny/2)) + source(j);

    Txx(fsx,fsy)=Txx(fsx,fsy) + source(j);
    Tyy(fsx,fsy)=Tyy(fsx,fsy) + source(j);
%     Txy(fsx,fsy)=Txy(fsx,fsy) + source(j);
	
    Vx(2:nx,:)=Vx(2:nx,:)+ MP_Vx(2:nx,:).* ( (Txx(2:nx,:)-Txx(1:nx-1,:))./dx + (Txy(2:nx,2:ny+1)-Txy(2:nx,1:ny))./dy );	
    UV_Vx_x=MP_Vx_xb.*UV_Vx_x+MP_Vx_xa.*(Txx([2:npml nx-npml+2:nx],:)-Txx([1:npml-1 nx-npml+1:nx-1],:))./dx;
    Vx([2:npml nx-npml+2:nx],:)=Vx([2:npml nx-npml+2:nx],:) + MP_Vx([2:npml nx-npml+2:nx],:).*UV_Vx_x; %npml-1层pml,pml边界层Vx未衰减
    UV_Vx_y=MP_Vx_yb.*UV_Vx_y+MP_Vx_ya.*(Txy(2:nx,[2:npml+1 ny-npml+2:ny+1])-Txy(2:nx,[1:npml ny-npml+1:ny]))./dy;
    Vx(2:nx,[1:npml ny-npml+1:ny])=Vx(2:nx,[1:npml ny-npml+1:ny]) + MP_Vx(2:nx,[1:npml ny-npml+1:ny]).*UV_Vx_y;
    Vx([1 nx+1],:)=Vx([2 nx],:);
    
    Vy(:,2:ny)=Vy(:,2:ny)+MP_Vy(:,2:ny).* ( (Txy(2:nx+1,2:ny)-Txy(1:nx,2:ny))./dx + (Tyy(:,2:ny)-Tyy(:,1:ny-1))./dy );
    UV_Vy_x=MP_Vy_xb.*UV_Vy_x+MP_Vy_xa.*(Txy([2:npml+1 nx-npml+2:nx+1],2:ny)-Txy([1:npml nx-npml+1:nx],2:ny))./dx;
    Vy([1:npml nx-npml+1:nx],2:ny)=Vy([1:npml nx-npml+1:nx],2:ny) + MP_Vy([1:npml nx-npml+1:nx],2:ny).*UV_Vy_x;    
    UV_Vy_y=MP_Vy_yb.*UV_Vy_y+MP_Vy_ya.*(Tyy(:,[2:npml ny-npml+2:ny])-Tyy(:,[1:npml-1 ny-npml+1:ny-1]))./dy;
    Vy(:,[2:npml ny-npml+2:ny])=Vy(:,[2:npml ny-npml+2:ny]) + MP_Vy(:,[2:npml ny-npml+2:ny]).*UV_Vy_y;
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
	
    ob_data_Vx(:,:,j)=Vx(recv_x,recv_y);
    ob_data_Vy(:,:,j)=Vy(recv_x,recv_y);
    
	if mod(j,50)==0
		figure(3);
		set(gcf,'outerposition',get(0,'screensize')); %全屏
        subplot(2,1,1);
        imagesc(squeeze(Vx(npml+1:nx+1-npml,1+npml:ny-npml)));
%         line([1,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx-tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 隧道上界面
%         line([1,tunnel_length-npml],[round((nx+tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 隧道下界面
%         line([tunnel_length-npml,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 掌子面
%         line([301-npml,301-npml],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 垂直界面        
%         line([326,235],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 60度倾斜界面
%         line([360,141],[35,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 30度倾斜界面
%         pos = [237 12 40 40]; rectangle('Position',pos,'Curvature',[1 1]);
%         caxis([-0.001,0.001]);
        caxis([-10^-8,10^-8]);
        colorbar;
        shading interp;
        grid on;
        axis equal tight;
%         axis equal
        set(gca,'fontsize',15,'fontweight', 'Bold');
        title(['Vx at step = ',num2str(j), ',  time = ',num2str(j*dt),' sec']);
        
        subplot(2,1,2);
        imagesc(squeeze(Vy(npml+1:nx-npml,1+npml:ny+1-npml)));
%         line([1,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx-tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 隧道上界面
%         line([1,tunnel_length-npml],[round((nx+tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 隧道下界面
%         line([tunnel_length-npml,tunnel_length-npml],[round((nx-tunnel_width)/2-npml),round((nx+tunnel_width)/2-npml)],'LineWidth',2,'Color',[1 1 1]); % 掌子面
%         line([301-npml,301-npml],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 垂直界面        
%         line([326,235],[1,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 60度倾斜界面
%         line([360,141],[35,nx-2*npml],'LineWidth',2,'Color',[1 1 1]); % 30度倾斜界面
%         pos = [237 12 40 40]; rectangle('Position',pos,'Curvature',[1 1]);
%         colorbar('horiz');
%         caxis([-0.001,0.001]);
        caxis([-10^-8,10^-8]);
        colorbar;
        shading interp;
        grid on;
        axis equal tight;
%         axis equal
        set(gca,'fontsize',15,'fontweight', 'Bold');
        title(['Vy at step = ',num2str(j), ',  time = ',num2str(j*dt),' sec']);
        
        h=suptitle('Velocity component of PSV wave');
        set(h,'fontsize',18,'fontweight', 'Bold');
		writeVideo(myavi,getframe(gcf))
    end
%     if j==1400
%         break
%     end
end
close(myavi);
toc;
save ob_data.mat ob_data_Vx ob_data_Vy
% close all
