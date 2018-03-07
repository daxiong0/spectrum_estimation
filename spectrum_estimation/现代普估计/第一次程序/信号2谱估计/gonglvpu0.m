%����ͼ������50�ι����ף�BT����ƽ������ƽ����
clear all;clc;close all;
N=512;                                %���ݵ���
e=cell(50,1);                         
for i=1:50
    e{i}=normrnd(0,1,N,1)';         %������ֵΪ0������Ϊ1�ĸ�˹����
end
A=[1,-1.3817,1.5632,-0.8843,0.4096];     %ϵ��
B=[1,0.3544,0.3508,0.1736,0.2401];
[H,w] = freqz(B,A,256);          %�����źţ�w��ֵƵ�ʿ̶�
 H=abs(H);   
 psd=10*log10((H.*H));
 %%%%%%%%%%%%����ͼ��
nfft1=N/2+1;
w1=0:pi/nfft1:pi-pi/nfft1;                   %��ֵƵ�ʿ̶�
psd1=cell(50,1);
for i=1:50
    x=filter(B,A,e{i});                  %��֪�ź�

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
title(['N=',num2str(N),';50������ͼ�������׹��ƾ�ֵ']);
xlabel('��Ƶ��/rad/s');
ylabel('����/db');
legend('���ƾ�ֵ','��ֵ');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%welch��
nfft2=N/4;
noverlap=nfft2/2;                                   %�����ݷ�Ϊ7��
inc=nfft2-noverlap;
k=(N-nfft2)/inc+1;                   %���ݶ���
nfft3=nfft2/2+1;
window=boxcar(nfft2)';
w2=0:pi/nfft3:pi-pi/nfft3;       %welch����ֵƵ�ʿ̶�
psd3=cell(50,1);
for i=1:50
    x=filter(B,A,e{i});                  %��֪�ź�
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
title(['k=',num2str(k),';N=',num2str(N),';50��welch�������׹��ƾ�ֵ']);
xlabel('��Ƶ��/rad/s');
ylabel('����/db');
legend('���ƾ�ֵ','��ֵ');
hold off;
%%%%%%%%%%%%%%BT��
% w1=0:pi/nfft1:pi-pi/nfft1;
w3=0:pi/N:pi-pi/N;            %BT����ֵƵ�ʿ̶�
psd5=cell(50,1);
for i=1:50
    x=filter(B,A,e{i});                  %��֪�ź�
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
title(['N=',num2str(N),';50��BT�������׹��ƾ�ֵ']);
xlabel('��Ƶ��/rad/s');
ylabel('����/db');
legend('���ƾ�ֵ','��ֵ');
hold off;
%%%%%%%%%%%%%%5%%%%%%%%%%%ƽ����
num=2;     %1����bartlett����2����turkey����3����Parzen��
if num==1
    g=bartlett(N)';  %%%����bartlett��
    chuang='bartlett��';
elseif(num==2)
    g=tukeywin(N,0.5)';  %%%����turkey��
    chuang='turkey��';
else
    g=parzenwin(N)';
    chuang='Parzen��';
end
%%%%%%%%%%%%%%%%%%%%%%%%
psd7=cell(50,1);
for i=1:50
    x1=filter(B,A,e{i});                  %��֪�ź�
    x2=fft(x1);
       x3=(abs(x2).^2)/N;  %������:
       x4=(ifft(x3));     %����غ���
       x4=x4.*g;           %����غ�������ʱ�ʹ�
     F2=fft(x4);      %�ٶ�����֮�������غ���������Ҷ�任
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
title([chuang,';N=',num2str(N),';50��ƽ���������׹��ƾ�ֵ']);
xlabel('��Ƶ��/rad/s');
ylabel('����/db');
legend('���ƾ�ֵ','��ֵ');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
figure(2)
%%%%%%%%%%%%%%%%%%%%%%����ͼ������
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
title(['N=',num2str(N),';50�ι����׹�ֵ����']);
xlabel('��Ƶ��/rad/s');
ylabel('����/db');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%welch������
welfc=zeros(1,nfft3);
dianzhi1=zeros(1,50);
for i=1:nfft3
    for j=1:50
        dianzhi1(j)=psd3{j}(i);
    end
    welfc(i)=var(dianzhi1);
end
plot(w2,10*log10(welfc),'linewidth',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BT������
BTfc=zeros(1,N);
dianzhi2=zeros(1,50);
for i=1:N
    for j=1:50
        dianzhi2(j)=psd5{j}(i);
    end
    BTfc(i)=var(dianzhi2);
end
plot(w3,10*log10(BTfc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ƽ��������
pinhuafc=zeros(1,nfft1);
dianzhi3=zeros(1,50);
for i=1:nfft1
    for j=1:50
        dianzhi3(j)=psd7{j}(i);
    end
    pinhuafc(i)=var(dianzhi3);
end
plot(w1,10*log10(pinhuafc));
legend('����ͼ��','welch��','BT��','ƽ����');
hold off;