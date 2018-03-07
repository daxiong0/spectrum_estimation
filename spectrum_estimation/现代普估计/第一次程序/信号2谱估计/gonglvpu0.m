%周期图法；做50次功率谱；BT法；平均法；平滑法
clear all;clc;close all;
N=512;                                %数据点数
e=cell(50,1);                         
for i=1:50
    e{i}=normrnd(0,1,N,1)';         %产生均值为0，方差为1的高斯噪声
end
A=[1,-1.3817,1.5632,-0.8843,0.4096];     %系数
B=[1,0.3544,0.3508,0.1736,0.2401];
[H,w] = freqz(B,A,256);          %理想信号，w真值频率刻度
 H=abs(H);   
 psd=10*log10((H.*H));
 %%%%%%%%%%%%周期图法
nfft1=N/2+1;
w1=0:pi/nfft1:pi-pi/nfft1;                   %估值频率刻度
psd1=cell(50,1);
for i=1:50
    x=filter(B,A,e{i});                  %已知信号

    f=fft(x);
    F=f(1:nfft1);
    psd1{i}=(abs(F).^2)/N;
end
psd2=zeros(1,nfft1);
for i=1:50
    psd2=psd2+psd1{i};
end
figure(1)
subplot 221
plot(w1,10*log10(psd2/50));
grid on;
hold on;
plot(w,psd);
title(['N=',num2str(N),';50次周期图法功率谱估计均值']);
xlabel('角频率/rad/s');
ylabel('功率/db');
legend('估计均值','真值');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%welch法
nfft2=N/4;
noverlap=nfft2/2;                                   %将数据分为7段
inc=nfft2-noverlap;
k=(N-nfft2)/inc+1;                   %数据段数
nfft3=nfft2/2+1;
window=boxcar(nfft2)';
w2=0:pi/nfft3:pi-pi/nfft3;       %welch法估值频率刻度
psd3=cell(50,1);
for i=1:50
    x=filter(B,A,e{i});                  %已知信号
    psd3{i}=pwelch(x,window,noverlap,nfft2);
end
psd4=zeros(nfft3,1);
for i=1:50
    psd4=psd4+psd3{i};
end
subplot 222
plot(w2,10*log10(psd4/50));
grid on;
hold on;
plot(w,psd);
title(['k=',num2str(k),';N=',num2str(N),';50次welch法功率谱估计均值']);
xlabel('角频率/rad/s');
ylabel('功率/db');
legend('估计均值','真值');
hold off;
%%%%%%%%%%%%%%BT法
% w1=0:pi/nfft1:pi-pi/nfft1;
w3=0:pi/N:pi-pi/N;            %BT法估值频率刻度
psd5=cell(50,1);
for i=1:50
    x=filter(B,A,e{i});                  %已知信号
    y=xcorr(x,'biased');
    F1=fft(y);
    psd5{i}=abs(F1(1:N));
end
psd6=zeros(1,N);
for i=1:50
   psd6=psd6+psd5{i};
end
subplot 223
plot(w3,10*log10(psd6/50));
grid on;
hold on;
plot(w,psd);
title(['N=',num2str(N),';50次BT法功率谱估计均值']);
xlabel('角频率/rad/s');
ylabel('功率/db');
legend('估计均值','真值');
hold off;
%%%%%%%%%%%%%%5%%%%%%%%%%%平滑法
num=2;     %1代表bartlett窗，2代表turkey窗，3代表Parzen窗
if num==1
    g=bartlett(N)';  %%%产生bartlett窗
    chuang='bartlett窗';
elseif(num==2)
    g=tukeywin(N,0.5)';  %%%产生turkey窗
    chuang='turkey窗';
else
    g=parzenwin(N)';
    chuang='Parzen窗';
end
%%%%%%%%%%%%%%%%%%%%%%%%
psd7=cell(50,1);
for i=1:50
    x1=filter(B,A,e{i});                  %已知信号
    x2=fft(x1);
       x3=(abs(x2).^2)/N;  %功率谱:
       x4=(ifft(x3));     %自相关函数
       x4=x4.*g;           %自相关函数乘以时滞窗
     F2=fft(x4);      %再对修正之后的自相关函数做傅立叶变换
       psd7{i}=abs(F2(1:nfft1));
end
psd8=zeros(1,nfft1);
for i=1:50
   psd8=psd8+psd7{i};
end
subplot 224
plot(w1,10*log10(psd8/50));
grid on;
hold on;
plot(w,psd);
title([chuang,';N=',num2str(N),';50次平滑法功率谱估计均值']);
xlabel('角频率/rad/s');
ylabel('功率/db');
legend('估计均值','真值');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%方差
figure(2)
%%%%%%%%%%%%%%%%%%%%%%周期图法方差
zqtfc=zeros(1,nfft1);
dianzhi0=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi0(j)=psd1{j}(i);
    end
    zqtfc(i)=var(dianzhi0);
end
plot(w1,10*log10(zqtfc));
grid on;
hold on;
title(['N=',num2str(N),';50次功率谱估值方差']);
xlabel('角频率/rad/s');
ylabel('方差/db');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%welch法方差
welfc=zeros(1,nfft3);
dianzhi1=zeros(1,50);
for i=1:nfft3
    for j=1:50
        dianzhi1(j)=psd3{j}(i);
    end
    welfc(i)=var(dianzhi1);
end
plot(w2,10*log10(welfc),'linewidth',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BT法方差
BTfc=zeros(1,N);
dianzhi2=zeros(1,50);
for i=1:N
    for j=1:50
        dianzhi2(j)=psd5{j}(i);
    end
    BTfc(i)=var(dianzhi2);
end
plot(w3,10*log10(BTfc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%平滑法方差
pinhuafc=zeros(1,nfft1);
dianzhi3=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi3(j)=psd7{j}(i);
    end
    pinhuafc(i)=var(dianzhi3);
end
plot(w1,10*log10(pinhuafc));
legend('周期图法','welch法','BT法','平滑法');
hold off;