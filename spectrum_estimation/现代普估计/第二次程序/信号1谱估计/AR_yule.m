clear all;clc;close all;
%ARģ�Ͳ������Ʒ�;p=20;
%�ɸı����ݵ���������ȣ�����
fs=1;                                     %����Ƶ��
time=128;                               %����ʱ��;
t=0:1/fs:time-1/fs;                      %ʱ��̶�
nfft=129;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
% dingjie=b1*sin(2*pi*f1*t+a1(1))+b2*sin(2*pi*f2*t+a2(1));%���ö����ź�
% 
% objectfun=zeros(1,100);              %�洢ÿ�εľ������
% for i=1:100                          %��1��100��
%     [a,E]=arburg(dingjie,i);              
%     objectfun(i) =((t1+i+1)/(t1-i-1))*E;  %�洢ÿ�εľ������
%     if i==1
%          orderpredict = i;                
%     elseif (objectfun(i) > objectfun(i-1) )      
%         orderpredict = i-1;                   %���յĽ���
%         break;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����������Ľ�����burg�㷨����psd����
y1=cell(50,1);                   %���ڱ�������źŵ�����
p=cell(50,1);                      %���ڱ��湦���׹�ֵ
orderpredict=20;                     %ģ�ͽ���
 for i=1:50
     R=normrnd(0,1,N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
    y1{i}=b1*sin(2*pi*f1*t+a1(i))+b2*sin(2*pi*f2*t+a2(i))+R;
    p{i}=pyulear(y1{i},orderpredict);                 %��burg�㷨���ƹ�����
%     plot(w,p{i});
 end
p1=zeros(nfft,1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;                      %50��AR�׹�ֵ�ľ�ֵ
figure(1)
subplot 211
plot(w,10*log10(p2));
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');
grid on;
title(['p=',num2str(orderpredict),';N=',num2str(N),':snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��Yulewalker�������׹��ƾ�ֵ']);
hold on;
zhenzhi=zeros(1,nfft);
plot(w,zhenzhi,'r');
grid on;
line([0.21 0.21],[0 10],'color','r');
line([0.23 0.23],[0 15],'color','r');
legend('����ֵ','��ֵ','Location','NorthWest');  
hold off;
%%%%%%%%%%%%�󷽲�
ldfc=zeros(1,nfft);
dianzhi=zeros(1,50);
for i=1:nfft
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    ldfc(i)=var(dianzhi);
end
subplot 212
plot(w,10*log10(ldfc));
grid on;
title(['p=',num2str(orderpredict),';N=',num2str(N),':snr1=',num2str(snr1),';snr2=',num2str(snr2),';50��Yulewalker�������׹��Ʒ���']);
xlabel('��һ��Ƶ��/Hz');
ylabel('����/db');