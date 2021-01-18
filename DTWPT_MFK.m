% Reference: Wang L, Liu Z, Cao H, et al. Subband averaging kurtogram with dual-tree complex wavelet packet 
% transform for rotating machinery fault diagnosis[J]. Mechanical Systems and Signal Processing, 2020, 142: 106755.
% --------------------------
% Author: Wang Lei
% Email: wang_llei@163.com
function DTWPT_MFK(x,K,nlevel,Fs)
load dtcwpt_filters_long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输入：x是输入信号；K是平均的次数；nlevel是分解层数；Fs是采样频率。
%检查输入数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(x);
if n==1
    if m>1
        x=x';
    else
        error('输入必须为向量');
    end
end
    
N = length(x);
if mod(N,K)~=0
    error('输入必须是K的倍数');
end

if nargin < 4
    Fs = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算子带的峭度值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x - mean(x);
L=N/K;
% s=ones(K,L);
Kwav = zeros(nlevel+1,2^nlevel);
for i=1:K
    s=x(1+(i-1)*L:i*L);
    KQ = K_dtwpQ(s,nlevel);		% 子带的峭度值
    KQ = KQ.*(KQ>0);
    Kwav = Kwav+KQ/K;
end
Kwav = Kwav.*(Kwav>0);
freq_w=Fs*linspace(0,.5,2^nlevel);
level=0:nlevel;
figure()
imagesc(freq_w,level,Kwav),colorbar,[I,J,M] = max_IJ(Kwav);
xlabel('Frequency (Hz)'),set(gca,'ytick',level),ylabel('Level k')
fi = (ceil(J/2^(nlevel-I+1))-.5)*2^(-I+1)*0.5*Fs;
title(['(a) K_{max}=',num2str(round(10*M)/10),' @ Level ',num2str(I-1),', Bw= ',num2str(Fs*2^(-I)),'Hz, fc=',num2str(fi),'Hz'])
set(gca,'FontName','Times New Roman','FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%双树复小波包变换
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1=SRDTWPT(x,first_1,h,f,I-1);
y2=SRDTWPT(x,first_2,g,f,I-1);
for j=1:2^(I-1)
    if j~=ceil(J/2^(nlevel-I+1))
        y1{j}=zeros(size(y1{j}));
        y2{j}=zeros(size(y2{j}));
%     else
%         figure()
%         plot(abs(y1{j}+1i*y2{j}));
    end
end
x1 = ISRDTWPT(y1,first_1(:,end:-1:1),h(:,end:-1:1),f(:,end:-1:1));
x2 = ISRDTWPT(y2,first_2(:,end:-1:1),g(:,end:-1:1),f(:,end:-1:1));
sx=(x1+x2)/2;



env = abs(hilbert(sx)).^2;
S = abs(fft((env(:)-mean(env)).*hanning(length(env))/length(env),N));
f = linspace(0,.5*Fs,N/2);
t = [0:length(sx)-1]/Fs;
figure()
plot(t,sx,'b');title('(b) Filtered signal'),xlabel('Time (s)'),ylabel('Amplitude');
set(gca,'XLim',[0 t(end)]);
set(gca,'YLim',[-1 1]);
set(gca,'YTick',-1:0.5:1);
set(gcf,'position',[413.0000 147.4000 560.0000 233.6000]);
set(gca,'FontName','Times New Roman','FontSize',12);
figure()
plot(f,S(1:N/2),'b'),title('(c) SES'),xlabel('Frequency (Hz)'),xlim([f(1) f(end)]),ylabel('Amplitude');
set(gca,'XLim',[0 300]);
set(gcf,'position',[413.0000 147.4000 560.0000 233.6000]);
set(gca,'FontName','Times New Roman','FontSize',12);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 子程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = kurt(x) %计算输入序列的峭度值
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,J,M] = max_IJ(X)
% 输出矩阵X最大值及其位置.

[temp,tempI] = max(X);
[M,J] = max(temp);
I = tempI(J);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function KQ = K_dtwpQ(x,nlevel)		% 子带的峭度值
    load dtcwpt_filters_long;
    KQ=ones(nlevel+1,2^nlevel);
    y1=cell(nlevel,2^nlevel);
    y2=y1;
    z=y1;
    
    KQ(1,:)=kurt(hilbert(x))*KQ(1,:);
    for i=1:nlevel
        if i==1
            %==============================================================
            fil0=first_1(1,:);
            fil1=first_1(2,:);
            yy=afb(x,fil0,fil1);
            y1{1,1}=yy(1,:);
            y1{1,2}=yy(2,:);
            %第一层实部分解
            %==============================================================
            fil0=first_2(1,:);
            fil1=first_2(2,:);
            yy=afb(x,fil0,fil1);
            y2{1,1}=yy(1,:);
            y2{1,2}=yy(2,:);
            %第一层虚部分解
            %==============================================================
            z{1,1}=y1{1,1}+1i*y2{1,1};
            z{1,2}=y1{1,2}+1i*y2{1,2};
            KQ(i+1,1:2^(nlevel-i))=kurt(z{1,1})*KQ(i+1,1:2^(nlevel-i));
            KQ(i+1,2^(nlevel-i)+1:2^nlevel)=kurt(z{1,2})*KQ(i+1,2^(nlevel-i)+1:2^nlevel);
        elseif i==2
            %==============================================================
            fil0=h(1,:);
            fil1=h(2,:);
            yy=afb(y1{1,1},fil0,fil1);
            y1{2,1}=yy(1,:);
            y1{2,2}=yy(2,:);
            %低频第二层实部分解
            %==============================================================
            fil0=g(1,:);
            fil1=g(2,:);
            yy=afb(y2{1,1},fil0,fil1);
            y2{2,1}=yy(1,:);
            y2{2,2}=yy(2,:);
            %低频第二层虚部分解
            %==============================================================
            fil0=h(1,:);
            fil1=h(2,:);
            yy=afb(y1{1,2},fil0,fil1);
            y1{2,3}=yy(1,:);
            y1{2,4}=yy(2,:);
            y1=rearrange(y1,2,3,4);
            %高频第二层实部分解
            %==============================================================
            fil0=g(1,:);
            fil1=g(2,:);
            yy=afb(y2{1,2},fil0,fil1);
            y2{2,3}=yy(1,:);
            y2{2,4}=yy(2,:);
            y2=rearrange(y2,2,3,4);
            %高频第二层虚部分解
            %==============================================================
            z{2,1}=y1{2,1}+1i*y2{2,1};
            z{2,2}=y1{2,2}+1i*y2{2,2};
            z{2,3}=y1{2,3}+1i*y2{2,3};
            z{2,4}=y1{2,4}+1i*y2{2,4};
            KQ(i+1,1:2^(nlevel-i))=kurt(z{2,1})*KQ(i+1,1:2^(nlevel-i));
            KQ(i+1,2^(nlevel-i)+1:2^(nlevel-i)*2)=kurt(z{2,2})*KQ(i+1,2^(nlevel-i)+1:2^(nlevel-i)*2);
            KQ(i+1,2^(nlevel-i)*2+1:2^(nlevel-i)*3)=kurt(z{2,3})*KQ(i+1,2^(nlevel-i)*2+1:2^(nlevel-i)*3);
            KQ(i+1,2^(nlevel-i)*3+1:2^(nlevel-i)*4)=kurt(z{2,4})*KQ(i+1,2^(nlevel-i)*3+1:2^(nlevel-i)*4);
        elseif i>2
            %==============================================================
            for k=1:2^(i-1)
                if k==1 || k==2^(i-1)
                    fil0=h(1,:);
                    fil1=h(2,:);
                else
                    fil0=f(1,:);
                    fil1=f(2,:);
                end
                yy=afb(y1{i-1,k},fil0,fil1);
                y1{i,2*k-1}=yy(1,:);
                y1{i,2*k}=yy(2,:);
%                 N = binary_number(k-1);
%                 if mod(N,2)~=0
                if mod(k,2)==0
                    y1=rearrange(y1,i,2*k-1,2*k);
                end
                if k==1 || k==2^(i-1)
                    fil0=g(1,:);
                    fil1=g(2,:);
                else
                    fil0=f(1,:);
                    fil1=f(2,:);
                end
                yy=afb(y2{i-1,k},fil0,fil1);
                y2{i,2*k-1}=yy(1,:);
                y2{i,2*k}=yy(2,:);
%                 N = binary_number(k-1);
%                 if mod(N,2)~=0
                if mod(k,2)==0
                    y2=rearrange(y2,i,2*k-1,2*k);
                end
                z{i,2*k-1}=y1{i,2*k-1}+1i*y2{i,2*k-1};
                z{i,2*k}=y1{i,2*k}+1i*y2{i,2*k};
                KQ(i+1,2^(nlevel-i)*(2*k-2)+1:2^(nlevel-i)*(2*k-1))=kurt(z{i,2*k-1})*KQ(i+1,2^(nlevel-i)*(2*k-2)+1:2^(nlevel-i)*(2*k-1));
                KQ(i+1,2^(nlevel-i)*(2*k-1)+1:2^(nlevel-i)*2*k)=kurt(z{i,2*k})*KQ(i+1,2^(nlevel-i)*(2*k-1)+1:2^(nlevel-i)*2*k);
            end
        end
    end
end
function [x] = ISRDTWPT(y,h_first,h,f)
%inverse subband rearranged dual-tree wavelet packet transform
%author:wang lei, Xi’an Jiaotong University.
%e-mail：wang_llei@163.com
%the original DTWOT was developed by I. Bayram.
%full packet inverse transform
%y : the cell array arranged as in DTWPT
%h_first : first stage synthesis filters([h0_first;h1_first])
%h : dual-tree synthesis filters([h0;h1])
%f : the 'same' synthesis filters([f0;f1])

max_level = log2(size(y,2));
% k = size(y,2)/2;
% n = max_level;
% for k=1:2^(n-1)
%     N = binary_number(k-1);
%     if mod(N,2)~=0
%         y=rearrange(y,1,2*k-1,2*k);
%     end
% end


xx = y(1,:);
xx=r2w(xx);
if max_level>2
    for n = max_level:-1:3
        for k=1:2^(n-1)
            if mod(k,2^(n-2))==1
                fil0=h(1,:);
                fil1=h(2,:);
            else
                fil0=f(1,:);
                fil1=f(2,:);
            end
%             N = binary_number(k-1);
%             if mod(N,2)~=0
%                 xx=rearrange(xx,1,2*k-1,2*k);
%             end
            x2{k} = sfb([xx{2*k-1};xx{2*k}],fil0,fil1);
        end
        xx=x2;
    end
end

if max_level>1
    %second stage
    fil0=h(1,:);
    fil1=h(2,:);
%     xx=rearrange(xx,1,3,4);
    x2{1}=sfb([xx{1};xx{2}],fil0,fil1);
    x2{2}=sfb([xx{3};xx{4}],fil0,fil1);
    xx=x2;
end

%first stage
fil0=h_first(1,:);
fil1=h_first(2,:);
x=sfb([xx{1};xx{2}],fil0,fil1);
end

function [y] = SRDTWPT(x,h_first,h,f,max_level)
%subband rearranged dual-tree wavelet packet transform
%author:wang lei, Xi’an Jiaotong University.
%e-mail：wang_llei@163.com
%the original DTWOT was developed by I. Bayram.
%x : input
%h_first : first stage filters([h0_first;h1_first])
%h : dual-tree filters([h0;h1])
%f : the 'same' filters([f0;f1])
%max_level : maximum level 
%y : output cell array containing all of the branches)-- the last row
%gives the full packet. The subbands are in order of decreasing frequency. 
y=cell(max_level,2^max_level);
%first stage
fil0=h_first(1,:);
fil1=h_first(2,:);
yy=afb(x,fil0,fil1);
y{1,1}=yy(1,:);
y{1,2}=yy(2,:);

if max_level>1
    %second stage
    fil0=h(1,:);
    fil1=h(2,:);
    yy=afb(y{1,1},fil0,fil1);
    y{2,1}=yy(1,:);
    y{2,2}=yy(2,:);

    yy=afb(y{1,2},fil0,fil1);
    y{2,3}=yy(1,:);
    y{2,4}=yy(2,:);
%     y=rearrange(y,2,3,4);
end
if max_level>2
    for n=3:max_level
        for k=1:2^(n-1)
            if mod(k,2^(n-2))==1
                fil0=h(1,:);
                fil1=h(2,:);
            else
                fil0=f(1,:);
                fil1=f(2,:);
            end
            yy=afb(y{n-1,k},fil0,fil1);
            y{n,2*k-1}=yy(1,:);
            y{n,2*k}=yy(2,:);
%             N = binary_number(k-1);
%             if mod(N,2)~=0
%                 y=rearrange(y,n,2*k-1,2*k);
%             end
        end
    end
end

y = y(max_level,:);
y=w2r(y);
end

function [y] = afb(x,h0,h1)
%x : input
%h0, h1 : analysis filters
%y : output -> [lowpass_channel ; highpass_channel] 

fil=h0;
temp=conv(fil,x);
temp(1:length(temp)-length(x))=temp(1:length(temp)-length(x))+temp(length(x)+1:end);
take=temp(1:length(x));
y0=take(1:2:length(x));

fil=h1;
temp=conv(fil,x);
temp(1:length(temp)-length(x))=temp(1:length(temp)-length(x))+temp(length(x)+1:end);
take=temp(1:length(x));
y1=take(1:2:length(x));

y=[y0;y1];
end

function x = sfb(y,g0,g1)

% y : output from 'afb.m'
% g0,g1 : synthesis filters
% x : reconstructed input

x0 = zeros(1,length(y)*2);
x1 = zeros(1,length(y)*2);

x0(1:2:end) = y(1,1:end);
x1(1:2:end) = y(2,1:end);

x = x0;
fil = g0;
temp = conv(fil,x);
temp(1:length(temp)-length(x)) = temp(1:length(temp)-length(x))+temp(length(x)+1:end);
take = temp(1:length(x));
shift = length(fil);
take0 = take(mod((0:length(x)-1)+shift-1,length(x))+1);

x = x1;
fil = g1;
temp = conv(fil,x);
temp(1:length(temp)-length(x)) = temp(1:length(temp)-length(x))+temp(length(x)+1:end);
take = temp(1:length(x));
shift = length(fil);
take = take(mod((0:length(x)-1)+shift-1,length(x))+1);

x = take+take0;
end
function x=w2r(y)
%y=cell(1,8);for i=1:8,y(i)={i};end;x=w2r(y);
if iscell(y)
    x=y;
    n=length(y);
    k=log2(n);
    for i=1:n
        a=binary(i-1,k);
        b=zeros(size(a));
        b(1)=a(1);
        r=1+a(1)*2^(k-1);
        for j=2:k
            b(j)=mod(a(j-1)+a(j),2);
            r=r+b(j)*2^(k-j);
        end
        x(i)=y(r);
    end
else
    error('y must be cell!');
end
end

function x=r2w(y)
%y=cell(1,8);for i=1:8,y(i)={i};end;x=w2r(y);
if iscell(y)
    x=y;
    n=length(y);
    k=log2(n);
    for i=1:n
        a=binary(i-1,k);
        b=zeros(1,k);
        c=zeros(1,k);
        b(1)=a(1);
        c(1)=a(1);
        w=1+a(1)*2^(k-1);
        for j=2:k
            c(j)=c(j-1)+a(j);
            b(j)=mod(c(j),2);
            w=w+b(j)*2^(k-j);
        end
%         a
%         c
%         b
        x(i)=y(w);
    end
else
    error('y must be cell!');
end
end

function x=rearrange(y,k,i,j)
if iscell(y)
    x=y;
    x(k,i)=y(k,j);
    x(k,j)=y(k,i);
else
    error('y must be cell!');
end
end
function N = binary_number(i)
%author:wang lei, Xi’an Jiaotong University.
%e-mail：wang_llei@163.com
% Returns the number of nonzero coefficients of the binary expansion of i: 
% i = a(1)*2^(k-1) + a(2)*2^(k-2) + ... + a(k)

k=fix(log2(i))+1;
a = zeros(1,k);
temp = i;
for l = k-1:-1:0
   a(k-l) = fix(temp/2^l);
   temp = temp - a(k-l)*2^l;
end
% a
N=length(find(a~=0));
end

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
end
