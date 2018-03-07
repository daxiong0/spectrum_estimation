%ARMAģ�Ͳ�����ֵ;P=4;q=0
%�ɸı����N,�ͽ���pp��qq
clear all;clc;close all;
N=128;      %256����  
pp=4;    %��ĸ����
qq=0;       %���ӽ���
w1=0:pi/N:pi-pi/N;
Aa=[1 -1.3817 1.5632 -0.8843 0.4096];
Bb=1;
[H1,w] = freqz(Bb,Aa,256);          %�����ź�
 H1=abs(H1);
 psd=10*log10((H1.*H1));           %�����źŹ�����
 figure (1)
 subplot 211
 plot(w,psd);
 xlabel('��Ƶ��/rad/s');
 ylabel('����/db');
 hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Ts=1;

for m=1:50
	K=idpoly(Aa,Bb,Ts);            %��ȷ���Ĳ�������һ������ʽģ��
	e=iddata([],idinput(N,'rgs'));   %iddata Create a data object to encapsulate the input/output data and
                                        % their properties.
                                           %idinput Generates input signals for identification.
                                        %'RGS': Generates a Random, Gaussian Signal
                                        %ֻ������û�����
	data=sim(K,e);                   %��Kģ�ͷ���e����,���simulink��������
	M=[pp,qq];                       %����
	armadat=armax(data,M);            %armax  Estimate armax polynomial model using time domain data.
                                        %data is the time-domain estimation data given as an IDDATA object
                                        %armadat contains the estimated values for A, B, and C polynomials along with
                                        % their covariances and structure information.
	B(m,:)=armadat.C;                   %Bϵ��
	A(m,:)=armadat.A;                    %Aϵ��

	[H(m,:)]=freqz(B(m,:),A(m,:),N);        %Ƶ����Ӧ
end
ESY=mean(abs(H));                      %���ֵ
plot(w1,20*log10(ESY)');
legend('��ֵ','ARMA���ƾ�ֵ');
title(['p=',num2str(pp),';q=',num2str(qq),';N=',num2str(N),';ARMA�������׹��ƾ�ֵ']);
grid on;
hold off;
%%%%%%%%%%%%�󷽲�
armafc=zeros(1,N);
for i=1:N            
    armafc(i)=var(abs(H( : ,i)));
end
subplot 212
plot(w1,10*log10(armafc));
grid on;
title(['p=',num2str(pp),';q=',num2str(qq),';N=',num2str(N),';ARMA�������׹��Ʒ���']);
xlabel('��Ƶ��/rad/s');
ylabel('����/db');