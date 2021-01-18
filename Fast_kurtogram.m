function c = Fast_kurtogram(x,nlevel,Fs)
% Fast_Kurtogram(x,nlevel,Fs)
% Computes the fast kurtogram of signal x up to level 'nlevel' via a fast decimated filterbank tree.
% Maximum number of decomposition levels is log2(length(x)), but it is 
% recommended to stay by a factor 1/8 below this.
% Fs = sampling frequency of signal x (default is Fs = 1)
%
% --------------------------
% Reference: J. Antoni, Fast Computation of the Kurtogram for the Detection of Transient Faults, 
% Mechanical Systems and Signal Processing, Volume 21, Issue 1, 2007, pp.108-124.
% --------------------------
% Author: J. Antoni
% Last Revision: 12-2014
% --------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(x);
N2 = log2(N) - 7;
if nlevel > N2
   error('Please enter a smaller number of decomposition levels');
end

if nargin < 3
    Fs = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAST COMPUTATION OF THE KURTOGRAM (by means of wavelet packets or STFT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x - mean(x);

% Analytic generating filters
N = 16;			fc = .4;                        % a short filter is just good enough!
h = fir1(N,fc).*exp(2i*pi*(0:N)*.125);
n = 2:N+1;
g = h(1+mod(1-n,N)).*(-1).^(1-n);
%
N = fix(3/2*N);
h1 = fir1(N,2/3*fc).*exp(2i*pi*(0:N)*.25/3);
h2 = h1.*exp(2i*pi*(0:N)/6);
h3 = h1.*exp(2i*pi*(0:N)/3);
%
Kwav = K_wpQ(x,h,g,h1,h2,h3,nlevel,'kurt2');		% kurtosis of the complex envelope
Kwav = Kwav.*(Kwav>0);							% keep positive values only!
% size(Kwav)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICAL DISPLAY OF RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
Level_w = 1:nlevel;	Level_w = [Level_w;Level_w+log2(3)-1];	Level_w = Level_w(:); Level_w = [0 Level_w(1:2*nlevel-1)'];
freq_w = Fs*((0:3*2^nlevel-1)/(3*2^(nlevel+1)) + 1/(3*2^(2+nlevel)));
imagesc(freq_w,1:2*nlevel,Kwav),colorbar,[I,J,M] = max_IJ(Kwav);
xlabel('Frequency (Hz)'),set(gca,'ytick',1:2*nlevel,'yticklabel',round(Level_w*10)/10),ylabel('Level k')
fi = (J-1)/3/2^(nlevel+1);   fi = fi + 2^(-2-Level_w(I));
title(['(a) K_{max}=',num2str(round(10*M)/10),' @ Level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, fc=',num2str(Fs*fi),'Hz'])
set(gca,'FontName','Times New Roman','FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL FILTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lev=Level_w(I);
[c,Bw,fc] = Find_wav_kurt(x,h,g,h1,h2,h3,lev,fi,'kurt2',Fs);



% END OF THE MAIN ROUTINE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIST OF SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = kurt(x,opt)
% Computes the kurtosis of signal x
if strcmp(opt,'kurt2')
   if all(x == 0), K = 0;return;end
   x = x - mean(x);
   E = mean(abs(x).^2);
   if E < eps, K = 0; return;end
   K = mean(abs(x).^4)/E^2;
   if all(isreal(x))
      K = K - 3;							% real signal
   else
      K = K - 2;
   end
elseif strcmp(opt,'kurt1')
   if all(x == 0), K = 0;return;end
   x = x - mean(x);
   E = mean(abs(x));
   if E < eps, K = 0; return;end
   K = mean(abs(x).^2)/E^2;
   if all(isreal(x))
      K = K - 1.57;							% real signal
   else
      K = K - 1.27;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,J,M] = max_IJ(X)
% Returns the row and column indices of the maximum in matrix X.

[temp,tempI] = max(X);
[M,J] = max(temp);
I = tempI(J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = binary(i,k)
% Returns the coefficients of the binary expansion of i: 
% i = a(1)*2^(k-1) + a(2)*2^(k-2) + ... + a(k)

if i>=2^k
   error('i must be such that i < 2^k !!')
end
a = zeros(1,k);
temp = i;
for l = k-1:-1:0
   a(k-l) = fix(temp/2^l);
   temp = temp - a(k-l)*2^l;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = raylinv(p,b)
% Inverse of the Rayleigh cumulative distribution function with parameter B at the probabilities in P.

x = zeros(size(p));

k1 = find(b <= 0| p < 0 | p > 1);
if any(k1) 
    tmp   = NaN;
    x(k1) = tmp(ones(size(k1)));
end

k = find(p == 1);
if any(k)
    tmp  = Inf;
    x(k) = tmp(ones(size(k))); 
end

k = find(b > 0 & p > 0  &  p < 1);
if any(k)
    pk = p(k);
    bk = b(k);
    x(k) = sqrt((-2*bk .^ 2) .* log(1 - pk));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIST OF SUBROUTINES FOR THE WAVELET PACKET KURTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = K_wpQ(x,h,g,h1,h2,h3,nlevel,opt,level)
% K = K_wpQ(x,h,g,h1,h2,h3,nlevel)
% Computes the kurtosis K of the complete "binary-ternary" wavelet packet transform w of signal x, 
% up to nlevel, using the lowpass and highpass filters h and g, respectively. 
% The values in K are sorted according to the frequency decomposition.

L = floor(log2(length(x)));
if nargin == 8
   if nlevel >= L
      error('nlevel must be smaller !!');
   end
   level = nlevel;
end
x = x(:);										 % shapes the signal as column vector if necessary

[KD,KQ] = K_wpQ_local(x,h,g,h1,h2,h3,nlevel,opt,level);

K = zeros(2*nlevel,3*2^nlevel);
K(1,:) = KD(1,:);
for i = 1:nlevel-1
   K(2*i,:) = KD(i+1,:);
   K(2*i+1,:) = KQ(i,:);
end
K(2*nlevel,:) = KD(nlevel+1,:);


%--------------------------------------------------------------------------
function [K,KQ] = K_wpQ_local(x,h,g,h1,h2,h3,nlevel,opt,level)
% (Subroutine of K_wpQ)

% performs one analysis level into the analysis tree
[a,d] = DBFB(x,h,g);                    

N = length(a);                       
d = d.*(-1).^(1:N)';
Lh = length(h);
Lg = length(g);

K1 = kurt(a(Lh:end),opt);
K2 = kurt(d(Lg:end),opt);

if level > 2    % ternary decomposition
   [a1,a2,a3] = TBFB(a,h1,h2,h3);
   [d1,d2,d3] = TBFB(d,h1,h2,h3);
   
   Ka1 = kurt(a1(Lh:end),opt);
   Ka2 = kurt(a2(Lh:end),opt);
   Ka3 = kurt(a3(Lh:end),opt);
   Kd1 = kurt(d1(Lh:end),opt);
   Kd2 = kurt(d2(Lh:end),opt);
   Kd3 = kurt(d3(Lh:end),opt);
else
   Ka1 = 0;
   Ka2 = 0;
   Ka3 = 0;
   Kd1 = 0;
   Kd2 = 0;
   Kd3 = 0;
end

if level == 1
   K =[K1*ones(1,3),K2*ones(1,3)];
   KQ = [Ka1 Ka2 Ka3 Kd1 Kd2 Kd3];
end

if level > 1
   [Ka,KaQ] = K_wpQ_local(a,h,g,h1,h2,h3,nlevel,opt,level-1);
   [Kd,KdQ] = K_wpQ_local(d,h,g,h1,h2,h3,nlevel,opt,level-1);
   
   K1 = K1*ones(1,length(Ka));
   K2 = K2*ones(1,length(Kd));
   K = [K1 K2; Ka Kd];
   
   Long = 2/6*length(KaQ);
   Ka1 = Ka1*ones(1,Long);
   Ka2 = Ka2*ones(1,Long);
   Ka3 = Ka3*ones(1,Long);
   Kd1 = Kd1*ones(1,Long);
   Kd2 = Kd2*ones(1,Long);
   Kd3 = Kd3*ones(1,Long);
   KQ = [Ka1 Ka2 Ka3 Kd1 Kd2 Kd3; KaQ KdQ];
end

if level == nlevel    
   K1 = kurt(x,opt);        % kurtosis of the raw signal is computed here
   K = [K1*ones(1,length(K));K];
   
   [a1,a2,a3] = TBFB(x,h1,h2,h3);
   Ka1 = kurt(a1(Lh:end),opt);
   Ka2 = kurt(a2(Lh:end),opt);
   Ka3 = kurt(a3(Lh:end),opt);
   Long = 1/3*length(KQ);
   Ka1 = Ka1*ones(1,Long);
   Ka2 = Ka2*ones(1,Long);
   Ka3 = Ka3*ones(1,Long);   
   KQ = [Ka1 Ka2 Ka3; KQ(1:end-2,:)];
end

% ------------------------------------------------------------------------
function [a,d] = DBFB(x,h,g)
% Double-band filter-bank (binary decomposition).
%   [a,d] = DBFB(x,h,g) computes the approximation
%   coefficients vector a and detail coefficients vector d,
%   obtained by passing signal x though a two-band analysis filter-bank.
%   h is the decomposition low-pass filter and
%   g is the decomposition high-pass filter.
% (Subroutine of 'K_wpQ_local')

N = length(x);

% lowpass filter
a = filter(h,1,x);
a = a(2:2:N);
a = a(:);

% highpass filter
d = filter(g,1,x);
d = d(2:2:N);
d = d(:);

% ------------------------------------------------------------------------
function [a1,a2,a3] = TBFB(x,h1,h2,h3)
% Trible-band filter-bank (ternary decomposition).
% (Subroutine of 'K_wpQ_local')

N = length(x);

% lowpass filter
a1 = filter(h1,1,x);
a1 = a1(3:3:N);
a1 = a1(:);

% passband filter
a2 = filter(h2,1,x);
a2 = a2(3:3:N);
a2 = a2(:);

% highpass filter
a3 = filter(h3,1,x);
a3 = a3(3:3:N);
a3 = a3(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,Bw,fc,i] = Find_wav_kurt(x,h,g,h1,h2,h3,Sc,Fr,opt,Fs)
% [c,Bw,fc,i] = Find_wav_kurt(x,h,g,h1,h2,h3,Sc,Fr,opt2)
% Sc = -log2(Bw)-1 with Bw the bandwidth of the filter
% Fr is in [0 .5]

if nargin < 10
   Fs = 1;
end

level = fix(Sc) + (rem(Sc,1)>=0.5)*(log2(3)-1);
Bw = 2^(-level-1);
freq_w = (0:2^level-1)/(2^(level+1))+Bw/2;
[temp,J] = min(abs(freq_w-Fr));
fc = freq_w(J);
i = round((fc/Bw-1/2));

if rem(level,1) == 0
   acoeff = binary(i,level);
   bcoeff = [];
   temp_level = level;
else
   i2 = fix(i/3);
   temp_level = fix(level)-1;
   acoeff = binary(i2,temp_level);
   bcoeff = i-i2*3;
end
acoeff = acoeff(end:-1:1);

c = K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,temp_level);

kx = kurt(c,opt);
sig = median(abs(c))/sqrt(pi/2);
threshold = sig*raylinv(.999,1);

figure()
t = (0:length(x)-1)/Fs;
tc = linspace(t(1),t(end),length(c));
plot(tc,real(c),'b'),
title('(b) Filtered signal');
xlabel('Time (s)'),xlim([tc(1) tc(end)]);
ylabel('Amplitude');
set(gcf,'position',[413.0000 147.4000 560.0000 233.6000]);
set(gca,'FontName','Times New Roman','FontSize',12);

nfft = 2*ceil(length(c)/2);
env = abs(c).^2;
S = abs(fft((env(:)-mean(env)).*hanning(length(env))/length(env),nfft));
f = linspace(0,.5*Fs/2^level,nfft/2);
figure()
plot(f,S(1:nfft/2),'b'),title('(c) SES'),xlabel('Frequency (Hz)');
ylabel('Amplitude'),xlim([f(1) f(end)]);
set(gcf,'position',[413.0000 147.4000 560.0000 233.6000]);
set(gca,'FontName','Times New Roman','FontSize',12);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,level)
% c = K_wpQ_filt(x,h,g,h1,h2,h3,acoeff,bcoeff,level)
% Calculates the complete "binary-ternary" wavelet packet transform w of signal x, 
% using the lowpass and highpass filters h and g, respectively. 
% The WP coefficients are sorted according to the frequency decomposition.
% This version handles both real and analytical filters, but does not yield WP coefficients
% suitable for signal synthesis.
% -------------------------------
% (Subroutine of 'Find_wav_kurt')

nlevel = length(acoeff);
L = floor(log2(length(x)));
if nargin == 8
   if nlevel >= L
      error('nlevel must be smaller !!');
   end
   level = nlevel;
end
x = x(:);									% shapes the signal as column vector if necessary

if nlevel == 0
   if isempty(bcoeff)
      c = x;
   else
      [c1,c2,c3] = TBFB(x,h1,h2,h3);
      if bcoeff == 0;
         c = c1(length(h1):end);
      elseif bcoeff == 1;
         c = c2(length(h2):end);
      elseif bcoeff == 2;
         c = c3(length(h3):end);
      end
   end
else
   c = K_wpQ_filt_local(x,h,g,h1,h2,h3,acoeff,bcoeff,level);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = K_wpQ_filt_local(x,h,g,h1,h2,h3,acoeff,bcoeff,level)
% (Subroutine of function 'K_wpQ_filt')

[a,d] = DBFB(x,h,g);                % perform one analysis level into the analysis tree

N = length(a);                       
d = d.*(-1).^(1:N)';

if level == 1
   if isempty(bcoeff)
      if acoeff(level) == 0
         c = a(length(h):end);
      else
         c = d(length(g):end);
      end
   else
      if acoeff(level) == 0
         [c1,c2,c3] = TBFB(a,h1,h2,h3);
      else
         [c1,c2,c3] = TBFB(d,h1,h2,h3);
      end
      if bcoeff == 0;
         c = c1(length(h1):end);
      elseif bcoeff == 1;
         c = c2(length(h2):end);
      elseif bcoeff == 2;
         c = c3(length(h3):end);
      end     
   end
end

if level > 1
   if acoeff(level) == 0
      c = K_wpQ_filt_local(a,h,g,h1,h2,h3,acoeff,bcoeff,level-1);
   else
      c = K_wpQ_filt_local(d,h,g,h1,h2,h3,acoeff,bcoeff,level-1);
   end 
end





