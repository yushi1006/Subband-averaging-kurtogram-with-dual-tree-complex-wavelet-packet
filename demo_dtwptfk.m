clc;clear;close all;

%��з����ź�---------------------------------------------------------------
fs = 10000;                  % ����Ƶ��
load xf2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������
% figure()
% [freq_s,sig_n]=envelope(x,N,fs);
% plot(freq_s(2:end/2),sig_n(2:end/2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Ҷ��
% figure()
% a=abs(fft(x))/(N/2);
% plot(f(2:end),a(2:length(f)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure,plot(t,x,'b'),title('(c) Signal with Hidden Periodic Transients'),xlabel('Time (s)'),ylabel('Amplitude');
% set(gca,'XLim',[0 3.2768]);
% % set(gca,'XTick',-2:1:2);
% % set(gca,'YLim',[-1 1]);
% % set(gca,'YTick',-1:0.5:1);
% set(gcf,'position',[961.0000 238.6000 560.0000 224.8000])
% set(gca,'FontName','Times New Roman','FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlevel = 7;     % number of decomposition levels



tic;
DTWPT_MFK(x,8,nlevel,fs);
toc;

tic;
c = Fast_kurtogram(x,nlevel,fs);
toc;



