clear all;close all;clc;
fs=pi;                                     %采样频率
nfft1=129;
nfft2=128;
w1=0:fs/nfft1:fs-fs/nfft1;    %频率刻度
w2=0:fs/nfft2:fs-fs/nfft2;    %频率刻度

A=[1,-1.3817,1.5632,-0.8843,0.4096];
B=[1,0.3544,0.3508,0.1736,0.2401];
[H,w] = freqz(B,A,256);          %理想信号
 H=abs(H);
 psd=10*log10((H.*H));           %理想信号功率谱
 
load ld_m;
ld=p2;
load burg_m;
burg=p2;
load arma_m;
arma=ESY;
figure(1);
subplot 211
plot(w1,10*log10(ld));
hold on;grid on;
plot(w1,10*log10(burg),'g');
plot(w2,10*log10(arma));
plot(w,psd,'r');
legend('L-D法','Burg法','ARMA法','真值'); 
xlabel('角频率/rad/s');
ylabel('功率/db')
title('不同参数谱估计在相同条件下的偏倚对比')
hold off;

load ldfc;
load burgfc;
load armafc;
subplot 212
plot(w1,10*log10(ldfc));
hold on;grid on;
plot(w1,10*log10(burgfc),'g');
plot(w2,10*log10(armafc));
legend('L-D法','Burg法','ARMA法');  
title('不同参数谱估计在相同条件下的方差对比')
hold off;