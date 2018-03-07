clear all;clc;close all;
fs=1;
N=512;
nfft1=N/2+1;
w1=0:fs/2/nfft1:fs/2-fs/2/nfft1;      %周期图法、平滑法频率刻度
w2=0:fs/2/N:fs/2-fs/2/N;            %BT法频率刻度

nfft=N/4;                           %每一段的长度                
nfft1=nfft/2+1;                      
w3=0:fs/2/nfft1:fs/2-fs/2/nfft1;            %welch法频率刻度

load zhouqitu_mean
zqt=p2;
load BT_junzhi
BT=p2;
load welch_junzhi
welch=p2';
load pinhua_junzhi
pinhua=p2;
figure (1)
subplot 211
plot(w1,10*log10(zqt));
hold on;
plot(w2,10*log10(BT),'r');
plot(w3,10*log10(welch),'b');
plot(w1,10*log10(pinhua),'g');
hold off;
grid on;
legend('周期图法','BT法','welch法','平滑法');
xlabel('归一化频率/Hz');
ylabel('功率/db');
title('N=512;相同条件下不同经典谱估计50次均值');
load zhouqitu_fc
load BT_fc
load welch_fc
load pinhua_fc
subplot 212
plot(w1,10*log10(zqtfc));
hold on;
plot(w2,10*log10(BTfc),'r');
plot(w3,10*log10(welfc),'b');
plot(w1,10*log10(pinhuafc),'g');
hold off;
grid on;
legend('周期图法','BT法','welch法','平滑法');
xlabel('归一化频率/Hz');
ylabel('方差/db');
title('N=512;相同条件下不同经典谱估计50次方差');


