clear all;close all;clc;
fs=pi;                                     %����Ƶ��
nfft1=129;
nfft2=128;
w1=0:fs/nfft1:fs-fs/nfft1;    %Ƶ�ʿ̶�
w2=0:fs/nfft2:fs-fs/nfft2;    %Ƶ�ʿ̶�

A=[1,-1.3817,1.5632,-0.8843,0.4096];
B=[1,0.3544,0.3508,0.1736,0.2401];
[H,w] = freqz(B,A,256);          %�����ź�
 H=abs(H);
 psd=10*log10((H.*H));           %�����źŹ�����
 
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
legend('L-D��','Burg��','ARMA��','��ֵ'); 
xlabel('��Ƶ��/rad/s');
ylabel('����/db')
title('��ͬ�����׹�������ͬ�����µ�ƫ�жԱ�')
hold off;

load ldfc;
load burgfc;
load armafc;
subplot 212
plot(w1,10*log10(ldfc));
hold on;grid on;
plot(w1,10*log10(burgfc),'g');
plot(w2,10*log10(armafc));
legend('L-D��','Burg��','ARMA��');  
title('��ͬ�����׹�������ͬ�����µķ���Ա�')
hold off;