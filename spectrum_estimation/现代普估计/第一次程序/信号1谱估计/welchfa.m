clear all;clc;close all; 
%改进的周期图法（welch法）
%参数设置
fs=1;
time=512;                 %持续时间
t=0:1/fs:time-1/fs;       % 时间刻度              
N=time*fs;                 %数据点数
f1=0.21;      
f2=0.23;
snr1=10;                  %信号1信噪比
snr2=15;   
b1=sqrt(2*10^(snr1/10));           %正弦波幅度
b2=sqrt(2*10^(snr2/10));
a1=2*pi*rand(1,50);          %产生0-2pi范围内的50个均匀随机数
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%psd估计
 p=cell(50,1);                     %保存psd
nfft=N/4;                           %每一段的长度
noverlap=nfft/2;                      %  两段重叠的数据长度
inc=nfft-noverlap;                   
k=(N-nfft)/inc+1;                   %数据分成7段
nfft1=nfft/2+1;                
window=boxcar(nfft)';         
w=0:fs/2/nfft1:fs/2-fs/2/nfft1;            %频率刻度
for i=1:50
     R=normrnd(0,1,	N,1)';         %产生均值为0，方差为1的高斯噪声
    y=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    p{i}=pwelch(y,window,noverlap,nfft);
end
p1=zeros(nfft1,1); 
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                  %谱估计均值
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('归一化频率/Hz');
ylabel('功率/db');
title(['k=',num2str(k),';N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50次welch法功率谱估计均值']);
grid on;
hold on;
zhenzhi=zeros(1,nfft1);
plot(w,zhenzhi,'r');
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('估计均值','真值','Location','northeast');  
hold off;
%%%%%%%%%%%%求方差
welfc=zeros(1,nfft1);
dianzhi=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    welfc(i)=var(dianzhi);
end
subplot 212
plot(w,10*log10(welfc));
title(['k=',num2str(k),';N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50次welch法功率谱估值方差']);
xlabel('归一化频率/Hz');
ylabel('方差/db');
grid on;