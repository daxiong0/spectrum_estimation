clear all;clc;close all;
%����ͼ��
%��������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=1;
time=512;                 %����ʱ��
t=0:1/fs:time-1/fs;         %ʱ��̶�               
N=time*fs;                %��������
f1=0.21;                    %��һ��Ƶ��
f2=0.23;
snr1=10;                    %����Ȳ���
snr2=15;
b1=sqrt(2*10^(snr1/10));        %�����źŷ���
b2=sqrt(2*10^(snr2/10));          
a1=2*pi*rand(1,50);          %����0-2pi��Χ�ڵ�50�����������
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ͼ�������׹���
 p=cell(50,1);                %����ÿһ������ͼ������Ĺ�����                   
nfft1=N/2+1;
w=0:fs/2/nfft1:fs/2-fs/2/nfft1;      %Ƶ�ʿ̶�
for i=1:50
    R=normrnd(0,1,N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
    y=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;       %��ʱ�������ɵİ�����+�����ź�
    F=fft(y);                         %��ʱ�����źŵĸ���Ҷ�任
    p{i}=F(1:nfft1);                   
    p{i}=(abs(p{i}).^2)/N;             %���湦����
end
p1=zeros(1,nfft1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                   %�����׾�ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%��ͼ
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
title(['N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50������ͼ�������׹��ƾ�ֵ']);
grid on;
hold on;
zhenzhi=zeros(1,nfft1);
plot(w,zhenzhi,'r');
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('���ƾ�ֵ','��ֵ','Location','northeast');  
hold off;
%%%%%%%%%%%%�󷽲�
zqtfc=zeros(1,nfft1);
dianzhi=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    zqtfc(i)=var(dianzhi);
end
%��ͼ
subplot 212
plot(w,10*log10(zqtfc));
title(['N=',num2str(N),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50������ͼ�������׹�ֵ����']);
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
grid on;