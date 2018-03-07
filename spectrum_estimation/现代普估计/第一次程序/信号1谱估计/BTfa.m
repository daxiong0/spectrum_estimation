clear all;clc;close all;
%BT����
%��������
fs=1;
time=512;                 %����ʱ��200s
t=0:1/fs:time-1/fs;                   
N=time*fs;
f1=0.21;
f2=0.23;
snr1=10;
snr2=15;
b1=sqrt(2*10^(snr1/10));
b2=sqrt(2*10^(snr2/10));
a1=2*pi*rand(1,50);          %����0-2pi��Χ�ڵ�50�����������
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BT��
p=cell(50,1);
w=0:fs/2/N:fs/2-fs/2/N;
for i=1:50
    R=normrnd(0,1,N,1)';               %������ֵΪ0������Ϊ1�ĸ�˹����
    y=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    y1=xcorr(y,'biased');             %������غ���
    F=fft(y1);                         %������غ�������Ҷ�任 
    p{i}=abs(F(1:N));                   %����һ��Ĺ�����
end
p1=zeros(1,N);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;
figure(1)
subplot 211
plot(w,10*log10(p2));
grid on;
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
title(['N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��BT�������׹��ƾ�ֵ']);
hold on;
zhenzhi=zeros(1,N);
plot(w,zhenzhi,'r');
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('���ƾ�ֵ','��ֵ','Location','NorthEast');  
hold off;
%%%%%%%%%%%%�󷽲�
BTfc=zeros(1,N);
dianzhi=zeros(1,50);
for i=1:N
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    BTfc(i)=var(dianzhi);
end
subplot 212
plot(w,10*log10(BTfc));
title(['N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��BT�������׹��Ʒ���']);
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
grid on;