clear all;clc;close all;
fs=1;                                     %采样频率
nfft1=129;
nfft2=128;
w1=0:fs/2/nfft1:fs/2-fs/2/nfft1;    %频率刻度
w2=0:fs/2/nfft2:fs/2-fs/2/nfft2;    %频率刻度
load ld_mean;
ld=p2;
load burg_mean;
burg=p2;
load arma_mean;
arma=p2;
figure(1);
subplot 211
plot(w1,10*log10(ld));
hold on;grid on;
plot(w1,10*log10(burg),'g');
plot(w2,10*log10(arma));
zhenzhi=zeros(1,nfft1);
plot(w1,zhenzhi,'r');
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('L-D法','Burg法','ARMA法','真值'); 
xlabel('归一化频率/Hz');
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