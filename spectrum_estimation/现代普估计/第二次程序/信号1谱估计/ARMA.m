clear all;clc;close all;
%ARMA模型参数估计法(两步最小二乘法)；p=12;q=8
%可改变点数time,和阶数pp，qq
qq=8;                       %分子阶数
pp=10;                        %分母阶数
fs=1;                                     %采样频率
time=128;                               %持续时间256s,信号周期为100s;
t=0:1/fs:time-1/fs;                      %时间刻度
nfft1=129;
w=0:fs/2/nfft1:fs/2-fs/2/nfft1;    %频率刻度
N=time*fs;                           %信号的采样点数        
f1=0.21;                                 %设置信号的频率
f2=0.23;
snr1=10;                                        %信噪比
snr2=15;
b1=sqrt(2*10^(snr1/10));
b2=sqrt(2*10^(snr2/10));                  %求出信号的幅度大小
a1=2*pi*rand(1,50);          %产生0-2pi范围内的50个均匀随机数
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1=cell(50,1);                   %用于保存带噪信号的数据
p=cell(50,1);                      %用于保存功率谱估值
L=27;                       %L=k+m
k=20;                           

 for i=1:50
     R=normrnd(0,1,N,1)';         %产生均值为0，方差为1的高斯噪声
    y1{i}=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    A=arburg(y1{i},k);
    e=filter(A,1,y1{i});          %得到白噪声
   z=y1{i}(L+1:time)';              
   e1=e(L+1:time)';
   Z=zeros(time-L,qq+pp);
   for j=1:pp
      Z( : ,j)=y1{i}(L-j+1:time-j)'; 
   end
   for j=1:qq
      Z( : ,pp+j)=-e(L-j+1:time-j)'; 
   end
   o=-inv(conj(Z')*Z)*(conj(Z')*z);
   Aa=[1 o(1:pp)'];                          %分母系数
   Bb=[1 o(pp+1:pp+qq)'];                       %分子系数
   e2=Z*o+z;
   fc=1/(time-L)*e2'*e2;                %噪声功率
   [H,w1]=freqz(Bb,Aa,time);             
    H=abs(H);
   p{i}=(H.*H)*fc;           %理想信号功率谱
%     plot(w,p{i});
 end
p1=zeros(N,1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                      %50次AR谱估值的均值
figure(1)
subplot 211
plot(w1/pi/2,10*log10(p2));
xlabel('归一化频率/Hz');
ylabel('功率/db');
title(['p=',num2str(pp),';q=',num2str(qq),';N=',num2str(N),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';ARMA法功率谱估计均值']);
hold on;
zhenzhi=zeros(1,nfft1);
plot(w,zhenzhi,'r');
grid on;
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('估计值','真值','Location','NorthWest');  
hold off;
%%%%%%%%%%%%求方差
armafc=zeros(1,N);
dianzhi=zeros(1,50);
for i=1:N
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    armafc(i)=var(dianzhi);
end
subplot 212
plot(w1/pi/2,10*log10(armafc));
grid on;
title(['p=',num2str(pp),';q=',num2str(qq),';N=',num2str(N),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';ARMA法功率谱估计方差']);
xlabel('归一化频率/Hz');
ylabel('方差/db');