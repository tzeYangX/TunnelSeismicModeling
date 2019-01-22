clear;
clc;

%AVO����
Vp=zeros(1,4);Vs=zeros(1,4);density=zeros(1,4);
Vp(1)=2450;Vp(2)=3000;Vp(3)=3500;Vp(4)=4000;
Vs(1)=1500;Vs(2)=2000;Vs(3)=2500;Vs(4)=3000;
density(1)=2;density(2)=2.2;density(3)=2.4;density(4)=2.6;
%{
for i=1:2
    Vp(i)=2450;
    Vs(i)=1500;
    density(i)=2;
end
for i=3:4
    Vp(i)=3000;
    Vs(i)=1700;
    density(i)=2.2;
end
%}
RShuey=zeros(length(Vp)-1,40);
RFatti=zeros(length(Vp)-1,40);
%e,f,h �ֱ�Ϊ��P���迹��S���迹���ܶȱ仯�ʵ�ϵ����
%e1,f1,h1 �ֱ�Ϊ��P���ٶȣ�S���ٶȺ��ܶȱ仯�ʵ�ϵ��,��(Vp2-Vp1)/(Vp2+Vp1)��
e=zeros(length(Vp)-1,40);
f=zeros(length(Vp)-1,40);
h=zeros(length(Vp)-1,40);
e1=zeros(length(Vp)-1,40);
f1=zeros(length(Vp)-1,40);
h1=zeros(length(Vp)-1,40);
for i=1:length(Vp)-1
    VpAv=(Vp(i+1)+Vp(i))/2;
    VpDiff=Vp(i+1)-Vp(i);
    VsAv=(Vs(i+1)+Vs(i))/2;
    VsDiff=Vs(i+1)-Vs(i);    
    densityAv=(density(i+1)+density(i))/2;
    densityDiff=density(i+1)-density(i);
    VsDVpSqu=(Vs(i+1).^2/Vp(i+1).^2+Vs(i).^2/Vp(i).^2)/2;
    IpDiff=Vp(i+1)*density(i+1)-Vp(i)*density(i);
    IsDiff=Vs(i+1)*density(i+1)-Vs(i)*density(i);
    IpAv=(Vp(i+1)*density(i+1)+Vp(i)*density(i))/2;
    IsAv=(Vs(i+1)*density(i+1)+Vs(i)*density(i))/2;
    A=0.5*(VpDiff/VpAv+densityDiff/densityAv);
    B=0.5*VpDiff/VpAv-4*VsDVpSqu*VsDiff/VsAv-2*VsDVpSqu*densityDiff/densityAv;
    C=0.5*VpDiff/VpAv;
    for j=1:40
        RShuey(i,j)=A+B*sin((j-1)*2*pi/180).^2+C*sin((j-1)*2*pi/180).^2*tan((j-1)*2*pi/180).^2;
        RFatti(i,j)=0.5*(IpDiff/IpAv)*(1+tan((j-1)*2*pi/180).^2)-4*VsDVpSqu*(IsDiff/IsAv)*sin((j-1)*2*pi/180).^2;
        RAki(i,j)=0.5*(1-4*VsDVpSqu*sin((j-1)*2*pi/180).^2)*densityDiff/densityAv+0.5/cos((j-1)*2*pi/180).^2*VpDiff/VpAv...
                  -4*VsDVpSqu*sin((j-1)*2*pi/180).^2*VsDiff/VsAv;
        e(i,j)=1+tan((j-1)*2*pi/180).^2;
        f(i,j)=-8*VsDVpSqu*sin((j-1)*2*pi/180).^2;
        h(i,j)=4*VsDVpSqu*sin((j-1)*2*pi/180).^2-tan((j-1)*2*pi/180).^2;
        e1(i,j)=1/cos((j-1)*2*pi/180).^2;
        f1(i,j)=-8*VsDVpSqu*sin((j-1)*2*pi/180).^2;
        h1(i,j)=(1-4*VsDVpSqu*sin((j-1)*2*pi/180).^2);                  
    end
end
figure(1);
title('����ϵ��������Ǳ仯���ߣ�Ӧ��Aki-Richards���ƹ�ʽ��');
xlabel('�����'),ylabel('����ϵ��');
hold on;
for i=1:40
    %plot(i,RShuey(2,i));
    %plot(i,RFatti(2,i));   
    plot(i,RAki(2,i),'b*');
end

%AVO����
%��ȡd����
G=zeros(3*40,3*3);
G1=zeros(3*40,3);
G2=zeros(3*40,3);
G3=zeros(3*40,3);
for i=1:3*40
    temp1=0;temp2=0;
    if (i/3-fix(i/3))==0
        temp2=fix(i/3);
    else
        temp2=fix(i/3)+1;
    end
    if mod(i,3)==0
        temp1=3;
    else
        temp1=mod(i,3);
    end
    R(i,1)=RAki(temp1,temp2);
end

%��ȡG����
for j=1:40
    for i=1:3
         G1(3*(j-1)+i,i)=e1(i,j);
         G2(3*(j-1)+i,i)=f1(i,j);
         G3(3*(j-1)+i,i)=h1(i,j);
    end
end

G=[G1 G2 G3];
%����ֵ�ֽ�
[U S V]=svd(G);
%d=Gm,��R(�൱��d)��G��Rpsd(�൱��m)
Rpsd=V*pinv(S)*U'*R;

%ʹ�õ�����Ʒ�������P���ٶȱ仯�ʣ�S���ٶȱ仯�ʺ��ܶȱ仯�����P���ٶȣ�S���ٶȺ��ܶȡ�
invertedVp(1)=Vp(1);
invertedVs(1)=Vs(1);
inverteddensity(1)=density(1);
for i=1:3
    temp=1;
    temp1=1;
    temp2=1;
    for j=1:i
        temp=temp*(1+Rpsd(j))/(1-Rpsd(j));
        temp1=temp1*(1+Rpsd(j+3))/(1-Rpsd(j+3));
        temp2=temp2*(1+Rpsd(j+6))/(1-Rpsd(j+6));
    end
    invertedVp(i+1)=invertedVp(1)*temp;
    invertedVs(i+1)=invertedVs(1)*temp1;
    inverteddensity(i+1)=inverteddensity(1)*temp2;
end
figure(2);
subplot(231);
xlabel('�������'),ylabel('ԭʼ�ݲ��ٶ�');
hold on;
stem(Vp);
subplot(234);
xlabel('�������'),ylabel('���ݺ��ݲ��ٶ�');
hold on;
stem(invertedVp);
subplot(232);
xlabel('�������'),ylabel('ԭʼ�Შ�ٶ�');
hold on;
stem(Vs);
subplot(235);
xlabel('�������'),ylabel('���ݺ�Შ�ٶ�');
hold on;
stem(invertedVs);
subplot(233);
xlabel('�������'),ylabel('ԭʼ�ܶ�');
hold on;
stem(density);
subplot(236);
xlabel('�������'),ylabel('���ݺ��ܶ�');
hold on;
stem(inverteddensity);
gtext('����ǰ���ٶȡ��ܶȶԱ�');



    





        
        

    




    