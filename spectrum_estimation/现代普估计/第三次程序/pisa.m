clear all;clc;close all;
%pisarenko���Ʒ�
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
    xt=y1{i}';
    sin_num=2;
    N=length(xt);
    rxx=xcorr(xt,'biased');
    rxx=rxx(N:(2*sin_num)+N);
    %Frequencies estimation
    Rxx=toeplitz(rxx);
    ev=eig(Rxx);
    [S k]=min(ev);
    [V D]=eig(Rxx);
    a=V(:,k);
    rts=roots(a);
    w_est=[];
    for j=1:4
        w_est(j)= angle(rts(j));
    end
    p{i}=(w_est/(2*pi))';
    
end
p1=zeros(4,1);                      %���ڱ��湦���׹�ֵ
for i=1:50
    p1=p1+p{i};
end
F=p1'./50;
j=1;
for i=1:4
    if F(i)<0
        continue;
    else
      F1(j)=F(i);
      j=j+1;
    end
end
figure(1)

line([F1(1) F1(1)],[0 15],'color','b'); 
hold on;
line([F1(2) F1(2)],[0 10],'color','b'); 
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
grid on;
title(['N=',num2str(N),':snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��pisarenko�������׹��ƾ�ֵ']);

zhenzhi=zeros(1,nfft);
plot(w,zhenzhi,'r');
grid on;
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');

hold off;

