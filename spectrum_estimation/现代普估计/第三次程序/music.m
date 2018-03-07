clear all;clc;close all;
%music���Ʒ�
%�ɸı����ݵ����������
fs=1;                                     %����Ƶ��
time=4096;                               %����ʱ��;
t=0:1/fs:time-1/fs;                      %ʱ��̶�
NFFT=256;
nfft=NFFT/2+1;
w=0:fs/2/nfft:fs/2-fs/2/nfft;    %Ƶ�ʿ̶�
N=time*fs;                           %�źŵĲ�������        
f1=0.21;                                 %�����źŵ�Ƶ��
f2=0.23;
snr1=10;                                        %�����
snr2=15;
b1=sqrt(2*10^(snr1/10));
b2=sqrt(2*10^(snr2/10));                  %����źŵķ��ȴ�С
a1=2*pi*rand(1,50);          %����0-2pi��Χ�ڵ�50�����������
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��music�㷨����psd����
y1=cell(50,1);                   %���ڱ�������źŵ�����
p=cell(50,1);                      %���ڱ��湦���׹�ֵ

 for i=1:50
     R=normrnd(0,1,N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
    y1{i}=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    X=corrmtx(y1{i},12,'mod');   %��������ؾ���
       p{i}=pmusic(X,4,'half');      % Uses the default NFFT of 256.
%     plot(w,10*log10(p{i}));
 end
p1=zeros(nfft,1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                      %50��music�׹�ֵ�ľ�ֵ
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
grid on;
title(['N=',num2str(N),':snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��music�������׹��ƾ�ֵ']);
hold on;
zhenzhi=zeros(1,nfft);
plot(w,zhenzhi,'r');
grid on;
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('����ֵ','��ֵ','Location','NorthWest');  
hold off;
%%%%%%%%%%%%�󷽲�
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
title(['N=',num2str(N),':snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��music�������׹��Ʒ���']);
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');