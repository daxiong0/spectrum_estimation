%ARģ�Ͳ������Ʒ���p=4;
clear all;clc;close all;
pp=40;     %����
N=128      %256����  
nfft1=129;
w1=0:pi/nfft1:pi-pi/nfft1;
noise=cell(50,1);
for i=1:50
    noise{i}=normrnd(0,1,N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
end
A=[1 -1.3817 1.5632 -0.8843 0.4096];
B=1;
[H,w] = freqz(B,A,256);          %�����ź�
 H=abs(H);
 psd=10*log10((H.*H));           %�����źŹ�����
 %%%%%%%%%%%%
%  dingjie=filter(B,A,noise{1});           %����   
% objectfun=zeros(1,100);              %�洢ÿ�εľ������
% for i=1:100                          %��1��100��
%     [a,E]=arburg(dingjie,i);              
%     objectfun(i) =((N+i+1)/(N-i-1))*E;  %�洢ÿ�εľ������
%     if i==1
%          orderpredict = i;                
%     elseif (objectfun(i) > objectfun(i-1) )      
%         orderpredict = i-1;                   %���յĽ���
%         break;
%     end
% end
y1=cell(50,1);                   %���ڱ�������źŵ�����
p=cell(50,1);                      %���ڱ��湦���׹�ֵ

 for i=1:50
%      R=normrnd(0,1,N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
    y1{i}=filter(B,A,noise{i}); 
    p{i}=pburg(y1{i},pp);                 %��burg�㷨���ƹ�����
%     plot(w,p{i});
 end
p1=zeros(nfft1,1);
for i=1:50
   p1=p1+p{i};
end
p2=p1/50;  
figure(1)
subplot 211
plot(w1,10*log10(p2));
grid on;
xlabel('��Ƶ��/rad/s');
ylabel('����/db');
title(['p=',num2str(pp),';N=',num2str(N),';50��Burg�������׹��ƾ�ֵ']);
hold on;
plot(w,psd,'r');
legend('���ƾ�ֵ','��ֵ','Location','NorthEast');  
hold off;
%%%%%%%%%%%%�󷽲�
burgfc=zeros(1,nfft1);
dianzhi=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi(j)=p{j}(i);
    end
    burgfc(i)=var(dianzhi);
end
subplot 212
plot(w1,10*log10(burgfc));
grid on;
title(['p=',num2str(pp),';N=',num2str(N),';50��Burg�������׹��Ʒ���']);
xlabel('��Ƶ��/rad/s');
ylabel('����/db');


