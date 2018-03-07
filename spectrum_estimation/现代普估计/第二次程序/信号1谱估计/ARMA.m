clear all;clc;close all;
%ARMAģ�Ͳ������Ʒ�(������С���˷�)��p=12;q=8
%�ɸı����time,�ͽ���pp��qq
qq=8;                       %���ӽ���
pp=10;                        %��ĸ����
fs=1;                                     %����Ƶ��
time=128;                               %����ʱ��256s,�ź�����Ϊ100s;
t=0:1/fs:time-1/fs;                      %ʱ��̶�
nfft1=129;
w=0:fs/2/nfft1:fs/2-fs/2/nfft1;    %Ƶ�ʿ̶�
N=time*fs;                           %�źŵĲ�������        
f1=0.21;                                 %�����źŵ�Ƶ��
f2=0.23;
snr1=10;                                        %�����
snr2=15;
b1=sqrt(2*10^(snr1/10));
b2=sqrt(2*10^(snr2/10));                  %����źŵķ��ȴ�С
a1=2*pi*rand(1,50);          %����0-2pi��Χ�ڵ�50�����������
a2=2*pi*rand(1,50); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1=cell(50,1);                   %���ڱ�������źŵ�����
p=cell(50,1);                      %���ڱ��湦���׹�ֵ
L=27;                       %L=k+m
k=20;                           

 for i=1:50
     R=normrnd(0,1,N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
    y1{i}=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    A=arburg(y1{i},k);
    e=filter(A,1,y1{i});          %�õ�������
   z=y1{i}(L+1:time)';              
   e1=e(L+1:time)';
   Z=zeros(time-L,qq+pp);
   for j=1:pp
      Z( : ,j)=y1{i}(L-j+1:time-j)'; 
   end
   for j=1:qq
      Z( : ,pp+j)=-e(L-j+1:time-j)'; 
   end
   o=-inv(conj(Z')*Z)*(conj(Z')*z);
   Aa=[1 o(1:pp)'];                          %��ĸϵ��
   Bb=[1 o(pp+1:pp+qq)'];                       %����ϵ��
   e2=Z*o+z;
   fc=1/(time-L)*e2'*e2;                %��������
   [H,w1]=freqz(Bb,Aa,time);             
    H=abs(H);
   p{i}=(H.*H)*fc;           %�����źŹ�����
%     plot(w,p{i});
 end
p1=zeros(N,1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                      %50��AR�׹�ֵ�ľ�ֵ
figure(1)
subplot 211
plot(w1/pi/2,10*log10(p2));
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
title(['p=',num2str(pp),';q=',num2str(qq),';N=',num2str(N),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';ARMA�������׹��ƾ�ֵ']);
hold on;
zhenzhi=zeros(1,nfft1);
plot(w,zhenzhi,'r');
grid on;
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('����ֵ','��ֵ','Location','NorthWest');  
hold off;
%%%%%%%%%%%%�󷽲�
armafc=zeros(1,N);
dianzhi=zeros(1,50);
for i=1:N
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    armafc(i)=var(dianzhi);
end
subplot 212
plot(w1/pi/2,10*log10(armafc));
grid on;
title(['p=',num2str(pp),';q=',num2str(qq),';N=',num2str(N),';snr1=',num2str(snr1),';snr2=',num2str(snr2),';ARMA�������׹��Ʒ���']);
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');