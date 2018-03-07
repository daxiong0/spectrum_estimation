%AR模型参数估计法；p=4;
clear all;clc;close all;
pp=40;     %阶数
N=128      %256个点  
nfft1=129;
w1=0:pi/nfft1:pi-pi/nfft1;
noise=cell(50,1);
for i=1:50
    noise{i}=normrnd(0,1,N,1)';         %产生均值为0，方差为1的高斯噪声
end
A=[1 -1.3817 1.5632 -0.8843 0.4096];
B=1;
[H,w] = freqz(B,A,256);          %理想信号
 H=abs(H);
 psd=10*log10((H.*H));           %理想信号功率谱
 %%%%%%%%%%%%
%  dingjie=filter(B,A,noise{1});           %定阶   
% objectfun=zeros(1,100);              %存储每次的均方误差
% for i=1:100                          %从1：100阶
%     [a,E]=arburg(dingjie,i);              
%     objectfun(i) =((N+i+1)/(N-i-1))*E;  %存储每次的均方误差
%     if i==1
%          orderpredict = i;                
%     elseif (objectfun(i) > objectfun(i-1) )      
%         orderpredict = i-1;                   %最终的阶数
%         break;
%     end
% end
y1=cell(50,1);                   %用于保存带噪信号的数据
p=cell(50,1);                      %用于保存功率谱估值

 for i=1:50
%      R=normrnd(0,1,N,1)';         %产生均值为0，方差为1的高斯噪声
    y1{i}=filter(B,A,noise{i}); 
    p{i}=pburg(y1{i},pp);                 %用burg算法估计功率谱
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
xlabel('角频率/rad/s');
ylabel('功率/db');
title(['p=',num2str(pp),';N=',num2str(N),';50次Burg法功率谱估计均值']);
hold on;
plot(w,psd,'r');
legend('估计均值','真值','Location','NorthEast');  
hold off;
%%%%%%%%%%%%求方差
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
title(['p=',num2str(pp),';N=',num2str(N),';50次Burg法功率谱估计方差']);
xlabel('角频率/rad/s');
ylabel('方差/db');


