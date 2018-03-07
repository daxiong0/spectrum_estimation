clear all;clc;close all;
%周期图法
%参数设置%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=1;
time=512;                 %持续时间
t=0:1/fs:time-1/fs;         %时间刻度               
N=time*fs;                %采样点数
f1=0.21;                    %归一化频率
f2=0.23;
snr1=10;                    %信噪比参数
snr2=15;
b1=sqrt(2*10^(snr1/10));        %正弦信号幅度
b2=sqrt(2*10^(snr2/10));          
a1=2*pi*rand(1,50);          %产生0-2pi范围内的50个均匀随机数
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%周期图法功率谱估计
 p=cell(50,1);                %保存每一次周期图法计算的功率谱                   
nfft1=N/2+1;
w=0:fs/2/nfft1:fs/2-fs/2/nfft1;      %频率刻度
for i=1:50
    R=normrnd(0,1,N,1)';         %产生均值为0，方差为1的高斯噪声
    y=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;       %临时保存生成的白噪声+正弦信号
    F=fft(y);                         %临时保存信号的傅里叶变换
    p{i}=F(1:nfft1);                   
    p{i}=(abs(p{i}).^2)/N;             %保存功率谱
end
p1=zeros(1,nfft1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                   %功率谱均值
%%%%%%%%%%%%%%%%%%%%%%%%%%作图
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('归一化频率/Hz');
ylabel('功率/db');
title(['N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50次周期图法功率谱估计均值']);
grid on;
hold on;
zhenzhi=zeros(1,nfft1);
plot(w,zhenzhi,'r');
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('估计均值','真值','Location','northeast');  
hold off;
%%%%%%%%%%%%求方差
zqtfc=zeros(1,nfft1);
dianzhi=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    zqtfc(i)=var(dianzhi);
end
%作图
subplot 212
plot(w,10*log10(zqtfc));
title(['N=',num2str(N),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50次周期图法功率谱估值方差']);
xlabel('归一化频率/Hz');
ylabel('方差/db');
grid on;