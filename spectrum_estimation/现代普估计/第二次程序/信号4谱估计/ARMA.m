%ARMA模型参数估值;P=4;q=0
%可改变点数N,和阶数pp，qq
clear all;clc;close all;
N=128;      %256个点  
pp=4;    %分母阶数
qq=0;       %分子阶数
w1=0:pi/N:pi-pi/N;
Aa=[1 -1.3817 1.5632 -0.8843 0.4096];
Bb=1;
[H1,w] = freqz(Bb,Aa,256);          %理想信号
 H1=abs(H1);
 psd=10*log10((H1.*H1));           %理想信号功率谱
 figure (1)
 subplot 211
 plot(w,psd);
 xlabel('角频率/rad/s');
 ylabel('功率/db');
 hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Ts=1;

for m=1:50
	K=idpoly(Aa,Bb,Ts);            %用确定的参数构造一个多项式模型
	e=iddata([],idinput(N,'rgs'));   %iddata Create a data object to encapsulate the input/output data and
                                        % their properties.
                                           %idinput Generates input signals for identification.
                                        %'RGS': Generates a Random, Gaussian Signal
                                        %只有输入没有输出
	data=sim(K,e);                   %用K模型仿真e数据,输出simulink仿真数据
	M=[pp,qq];                       %阶数
	armadat=armax(data,M);            %armax  Estimate armax polynomial model using time domain data.
                                        %data is the time-domain estimation data given as an IDDATA object
                                        %armadat contains the estimated values for A, B, and C polynomials along with
                                        % their covariances and structure information.
	B(m,:)=armadat.C;                   %B系数
	A(m,:)=armadat.A;                    %A系数

	[H(m,:)]=freqz(B(m,:),A(m,:),N);        %频率响应
end
ESY=mean(abs(H));                      %求均值
plot(w1,20*log10(ESY)');
legend('真值','ARMA估计均值');
title(['p=',num2str(pp),';q=',num2str(qq),';N=',num2str(N),';ARMA法功率谱估计均值']);
grid on;
hold off;
%%%%%%%%%%%%求方差
armafc=zeros(1,N);
for i=1:N            
    armafc(i)=var(abs(H( : ,i)));
end
subplot 212
plot(w1,10*log10(armafc));
grid on;
title(['p=',num2str(pp),';q=',num2str(qq),';N=',num2str(N),';ARMA法功率谱估计方差']);
xlabel('角频率/rad/s');
ylabel('方差/db');