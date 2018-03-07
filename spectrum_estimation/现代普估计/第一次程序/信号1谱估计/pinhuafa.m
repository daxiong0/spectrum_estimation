clear all;clc;close all;
%周期图平滑法
%参数设置
fs=1;
time=512;                 %持续时间
t=0:1/fs:time-1/fs;                   
N=time*fs;
f1=0.21;
f2=0.23;
snr1=10;
snr2=15;
b1=sqrt(2*10^(snr1/10));
b2=sqrt(2*10^(snr2/10));
%设置窗函数
num=3;                          %1代表bartlett窗，2代表turkey窗，3代表Parzen窗
a1=2*pi*rand(1,50);          %产生0-2pi范围内的50个均匀随机数
a2=2*pi*rand(1,50); 
%%%平滑法功率谱估计

 p=cell(50,1);
nfft=N;
nfft1=nfft/2+1;
w=0:fs/2/nfft1:fs/2-fs/2/nfft1;
if num==1
    g=bartlett(N)';  %%%产生bartlett窗
    chuang='bartlett窗';
elseif(num==2)    
    g=tukeywin(N,0.5)';  %%%产生turkey窗
    chuang='turkey窗';
else
    g=parzenwin(N)';%Parzen窗
    chuang='Parzen窗';
end

for i=1:50
     R=normrnd(0,1,N,1)';         %产生均值为0，方差为1的高斯噪声
    y=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    x=fft(y);
       x=(abs(x).^2)/N;  %功率谱:
       y=(ifft(x));     %自相关函数
       y=y.*g;           %自相关函数乘以时滞窗
       F=fft(y);      %再对修正之后的自相关函数做傅立叶变换
       p{i}=abs(F(1:nfft1));
end
p1=zeros(1,nfft1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('归一化频率/Hz');
ylabel('功率/db');
title([chuang,';N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50次平滑法功率谱估计均值']);
grid on;
hold on;
zhenzhi=zeros(1,nfft1);
plot(w,zhenzhi,'r');
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('估计均值','真值','Location','NorthEast');  
hold off;
%%%%%%%%%%%%求方差
pinhuafc=zeros(1,nfft1);
dianzhi=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    pinhuafc(i)=var(dianzhi);
end
subplot 212
plot(w,10*log10(pinhuafc));
title([chuang,';N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50次平滑法功率谱估计方差']);
xlabel('归一化频率/Hz');
ylabel('方差/db');
grid on;
