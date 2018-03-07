clear all;clc;close all;
%����ͼƽ����
%��������
fs=1;
time=512;                 %����ʱ��
t=0:1/fs:time-1/fs;                   
N=time*fs;
f1=0.21;
f2=0.23;
snr1=10;
snr2=15;
b1=sqrt(2*10^(snr1/10));
b2=sqrt(2*10^(snr2/10));
%���ô�����
num=3;                          %1����bartlett����2����turkey����3����Parzen��
a1=2*pi*rand(1,50);          %����0-2pi��Χ�ڵ�50�����������
a2=2*pi*rand(1,50); 
%%%ƽ���������׹���

 p=cell(50,1);
nfft=N;
nfft1=nfft/2+1;
w=0:fs/2/nfft1:fs/2-fs/2/nfft1;
if num==1
    g=bartlett(N)';  %%%����bartlett��
    chuang='bartlett��';
elseif(num==2)    
    g=tukeywin(N,0.5)';  %%%����turkey��
    chuang='turkey��';
else
    g=parzenwin(N)';%Parzen��
    chuang='Parzen��';
end

for i=1:50
     R=normrnd(0,1,N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
    y=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    x=fft(y);
       x=(abs(x).^2)/N;  %������:
       y=(ifft(x));     %����غ���
       y=y.*g;           %����غ�������ʱ�ʹ�
       F=fft(y);      %�ٶ�����֮�������غ���������Ҷ�任
       p{i}=abs(F(1:nfft1));
end
p1=zeros(1,nfft1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
title([chuang,';N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��ƽ���������׹��ƾ�ֵ']);
grid on;
hold on;
zhenzhi=zeros(1,nfft1);
plot(w,zhenzhi,'r');
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('���ƾ�ֵ','��ֵ','Location','NorthEast');  
hold off;
%%%%%%%%%%%%�󷽲�
pinhuafc=zeros(1,nfft1);
dianzhi=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    pinhuafc(i)=var(dianzhi);
end
subplot 212
plot(w,10*log10(pinhuafc));
title([chuang,';N=',num2str(time),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��ƽ���������׹��Ʒ���']);
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
grid on;
