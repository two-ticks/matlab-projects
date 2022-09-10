clc;
clear;

M = 5;                     % number of multipaths
N = 10^5;                  % number of samples to generate
Ts = 0.0001;               % sampling period in seconds
freq = 5*10^9;             % hz, maximum frequency
vr = 50;                   % m/s, speed of reciever 

fd = (vr * freq) / (3*10^8);  % maximum doppler spread in hertz

x=(0:N-1)*Ts;
initialSignal = 12*cos(20*x+pi*x) + 20*sin(0.2*x).*cos(50*x);

h = rayleighFading(M, N, fd, Ts); % transfer function of Rayleigh Channel
h_re = real(h); h_im = imag(h);

% output signal from channel
filteredSignal = ifft(fft(initialSignal).*fft(h));

% channel capacity 
% noise at room temperature
% 4·10^(-21) Watts/Hz

% channel_snr = snr(real(filteredSignal), 1/Ts);


function [h]=rayleighFading(M,N,fd,Ts)
% function to generate Rayleigh Fading samples based on Clarke's model
% M = number of multi-paths in the channel
% N = number of samples to generate
% fd = maximum Doppler frequency
% Ts = sampling period
a=0;
b=2*pi;
alpha=a+(b-a)*rand(1,M); %uniformly distributed from 0 to 2 pi
beta=a+(b-a)*rand(1,M);  %uniformly distributed from 0 to 2 pi
theta=a+(b-a)*rand(1,M); %uniformly distributed from 0 to 2 pi
m=1:M;
for n=1:N;
x=cos(((2.*m-1)*pi+theta)/(4*M));
h_re(n)=1/sqrt(M)*sum(cos(2*pi*fd*x*n'*Ts+alpha));
h_im(n)=1/sqrt(M)*sum(sin(2*pi*fd*x*n'*Ts+beta));
end
h=h_re+j*h_im;
end

