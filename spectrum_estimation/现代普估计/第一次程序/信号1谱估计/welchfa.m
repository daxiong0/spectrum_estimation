clear all;clc;close all; 
%�Ľ�������ͼ����welch����
%��������
fs=1;
time=512;                 %����ʱ��
t=0:1/fs:time-1/fs;       % ʱ��̶�              
N=time*fs;                 %���ݵ���
f1=0.21;      
f2=0.23;
snr1=10;                  %�ź�1�����
snr2=15;   
b1=sqrt(2*10^(snr1/10));           %���Ҳ�����
b2=sqrt(2*10^(snr2/10));
a1=2*pi*rand(1,50);          %����0-2pi��Χ�ڵ�50�����������
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%psd����
 p=cell(50,1);                     %����psd
nfft=N/4;                           %ÿһ�εĳ���
noverlap=nfft/2;                      %  �����ص������ݳ���
inc=nfft-noverlap;                   
k=(N-nfft)/inc+1;                   %���ݷֳ�7��
nfft1=nfft/2+1;                
window=boxcar(nfft)';         
w=0:fs/2/nfft1:fs/2-fs/2/nfft1;            %Ƶ�ʿ̶�
for i=1:50
     R=normrnd(0,1,	N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
    y=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    p{i}=pwelch(y,window,noverlap,nfft);
end
p1=zeros(nfft1,1); 
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                  %�׹��ƾ�ֵ
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
title(['k=',num2str(k),';N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��welch�������׹��ƾ�ֵ']);
grid on;
hold on;
zhenzhi=zeros(1,nfft1);
plot(w,zhenzhi,'r');
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('���ƾ�ֵ','��ֵ','Location','northeast');  
hold off;
%%%%%%%%%%%%�󷽲�
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
title(['k=',num2str(k),';N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��welch�������׹�ֵ����']);
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
grid on;