clear all;clc;close all;
%music估计法
%可改变数据点数，信噪比
fs=1;                                     %采样频率
time=4096;                               %持续时间;
t=0:1/fs:time-1/fs;                      %时间刻度
NFFT=256;
nfft=NFFT/2+1;
w=0:fs/2/nfft:fs/2-fs/2/nfft;    %频率刻度
N=time*fs;                           %信号的采样点数        
f1=0.21;                                 %设置信号的频率
f2=0.23;
snr1=10;                                        %信噪比
snr2=15;
b1=sqrt(2*10^(snr1/10));
b2=sqrt(2*10^(snr2/10));                  %求出信号的幅度大小
a1=2*pi*rand(1,50);          %产生0-2pi范围内的50个均匀随机数
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%用music算法进行psd估计
y1=cell(50,1);                   %用于保存带噪信号的数据
p=cell(50,1);                      %用于保存功率谱估值

 for i=1:50
     R=normrnd(0,1,N,1)';         %产生均值为0，方差为1的高斯噪声
    y1{i}=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    X=corrmtx(y1{i},12,'mod');   %生成自相关矩阵
       p{i}=pmusic(X,4,'half');      % Uses the default NFFT of 256.
%     plot(w,10*log10(p{i}));
 end
p1=zeros(nfft,1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                      %50次music谱估值的均值
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('归一化频率/Hz');
ylabel('功率/db');
grid on;
title(['N=',num2str(N),':snr1=',num2str(snr1),';snr2=',num2str(snr2),';50次music法功率谱估计均值']);
hold on;
zhenzhi=zeros(1,nfft);
plot(w,zhenzhi,'r');
grid on;
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('估计值','真值','Location','NorthWest');  
hold off;
%%%%%%%%%%%%求方差
burgfc=zeros(1,nfft);
dianzhi=zeros(1,50);
for i=1:nfft
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    burgfc(i)=var(dianzhi);
end
subplot 212
plot(w,10*log10(burgfc));
grid on;
title(['N=',num2str(N),':snr1=',num2str(snr1),';snr2=',num2str(snr2),';50次music法功率谱估计方差']);
xlabel('归一化频率/Hz');
ylabel('方差/db');