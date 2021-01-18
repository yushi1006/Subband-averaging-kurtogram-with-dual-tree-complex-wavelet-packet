% function y=envelope(signal,N,Fs)
function [freq_s,sig_n]=envelope(signal,N,Fs)

analy=hilbert(signal);
y=abs(analy);
% figure();
T=N/Fs;
sig_f=abs(fft(y(1:N)',N));
sig_n=sig_f/(norm(sig_f));
freq_s=(0:N-1)/T;
% plot(freq_s(2:round(end/2)),sig_n(2:round(end/2)));title('Envelope Detection : Hilbert Transform')
% plot(freq_s(2:500),sig_n(2:500));title('Envelope Detection : Hilbert Transform')


