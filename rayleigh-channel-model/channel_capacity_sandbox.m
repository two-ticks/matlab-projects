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

snrdB=-10:0.5:20; %Range of SNRs to simulate

% h = (randn(1,100) + 1i*randn(1,100) )/sqrt(2); % Rayleigh flat channel
sigma_z=1; %Noise power - assumed to be unity


snr = 10.^(snrdB/10); %SNRs in linear scale
P=(sigma_z^2)*snr./(mean(abs(h).^2)); %Calculate corresponding values for P

C_erg_awgn= (log2(1+ mean(abs(h).^2).*P/(sigma_z^2))); %AWGN channel capacity (Bound)
C_erg = mean((log2(1+ ((abs(h).^2).')*P/(sigma_z^2)))); %ergodic capacity for Fading channel

plot(snrdB,C_erg_awgn,'b'); hold on;
plot(snrdB,C_erg,'r'); grid on;
legend('AWGN channel capacity','Fading channel Ergodic capacity');
title('flat fading channel - Ergodic capacity');
xlabel('SNR (dB)');ylabel('Capacity (bps/Hz)');

% outage probability 

snrth = 1; 
s = zeros(1, N);
snrDB = 2:2:12;

p_sim = zeros(1, length(snrDB));
p1_ana = zeros(1, length(snrDB));

for jj = 1:length(snrDB)
    count = 0;
    snr = 10.^(snrDB(jj)/10);
    for i = 1:N 
        s(i) = (abs(h(i).^2)).*snr;
        if s(i) < snrth
            count = count + 1;
        end
    end
    p_sim(jj) = count/N;
p1_ana(jj) = 1 - exp(-snrth/snr).*(snrth./snr+1);    
end
semilogy(snrDB, p_sim,'b-', LineWidth=1); hold on; grid on;
semilogy(snrDB, p1_ana,'r o', LineWidth=1); hold on; grid on;

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

