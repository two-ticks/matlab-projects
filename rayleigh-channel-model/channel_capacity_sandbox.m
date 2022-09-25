%% v1 19 September, 2022  
% TODO : 
% [ ] check for corner/test cases of SNR
%     - [x] truncated capacity > inversion
%     - [x] water filling > other capacities 
% [ ] title and labels for graphs 

% clear 
clc;
clear;
close all;

% parameters
M             = 5;                                         % number of multipaths
N             = 10^5;                                      % number of samples to generate
Ts            = 0.0001;                                    % sampling period in seconds
freq          = 5*10^9;                                    % hz, maximum frequency
vr            = 50;                                        % m/s, speed of reciever 
B             = 2*freq;
fd            = (vr * freq) / (3*10^8);                    % maximum doppler spread in hertz

% rayleigh channel 
h             = rayleighFading(M, N, fd, Ts);              % transfer function of rayleigh channel
h_re          = real(h); h_im = imag(h);

% channel capacity 
snrdB         = -1:2:20;                                  % range of SNRs to simulate
sigma_z       = 1;                                         % noise power - assumed to be unity
snr           = 10.^(snrdB/10);                            % SNRs in linear scale
P             = (sigma_z^2)*snr./(mean(abs(h).^2));        % calculate corresponding values for Power


% normalized capacity: C/B
C_erg_awgn    =...
          (log2(1+ mean(abs(h).^2).*P/(sigma_z^2)));       % AWGN channel capacity (Bound)
C_erg         =...
    mean((log2(1+ ((abs(h).^2).')*P/(sigma_z^2))));        % ergodic capacity for Fading channel

figure
plot(snrdB, C_erg_awgn,'b'); hold on;
plot(snrdB, C_erg,'r'); grid on;

legend('AWGN channel capacity','Fading channel Ergodic capacity');
title('flat fading channel - Ergodic capacity');
xlabel('SNR (dB)');ylabel('Capacity (bps/Hz)');
hold off 

% outage probability 
snrTh = 1; 
s = zeros(1, N);                                           % initialize instanteneous snr

p_sim = zeros(1, length(snrdB));
p_ana = zeros(1, length(snrdB));

for jj = 1:length(snrdB)                                   % loop over all snr
    count = 0; 
    snr = 10.^(snrdB(jj)/10);
    for i = 1:N 
        s(i) = (abs(h(i).^2)).*snr;                        % instanteneous snr  
        if s(i) < snrTh
            count = count + 1;                             % counting instanteneous snr < threshold
        end
    end
p_sim(jj) = count/N;                                       % simulated probabilty
p_ana(jj) = 1 - exp(-snrTh/snr).*(snrTh./snr+1);           % analytical outage probability 
end
figure 
semilogy(snrdB, p_sim, 'b-', LineWidth=1); hold on; grid on; 
semilogy(snrdB, p_ana, 'r +', LineWidth=1); grid on;
hold off
C_out_sim = sum((p_sim).*(log2(1+snr)));
C_out_ana = sum((p_ana).*(log2(1+snr)));

% inversion channel

inv_snr_expectation = sum((p_sim)./snr);
C_inv = (log2(1+ 1/inv_snr_expectation));



% water filling 

% dependent parameters --------------------------------------------------
    % convert from dB to scalar
    signalToNoiseRatio      = 10.^(0.1*snrdB);
    % number of elements in SNR vector
    nOfsignalToNoiseRatio   = length(signalToNoiseRatio);
    % allocate memory for capacity 
    CapacityVec             = ...
        zeros(nOfsignalToNoiseRatio,N);
    noisePower = 1;

for j = 1 : nOfsignalToNoiseRatio
            % load SNR for this iteration
            snr = signalToNoiseRatio(j);    % SNR at this iteration
            % calaculate transmit power
            txPower = noisePower*snr;       % we know SNR = txPower/noisePower
                                            % with the assumption that
                                            % there is no pathloss. It
                                            % means here txPower = rxPower;
            for i = 1 : N     % loop over number of iterations
                % generate channel coefficient
                % find allocated power based on waterfilling algorithm
                allocatedPower = waterFilling(snr,txPower);
                % calculate the capacity
                capacity = sum(log2(1+allocatedPower.*snr));
                % store the value
                CapacityVec(j,i) = capacity;
            end
end

mean_capacity = zeros(1:nOfsignalToNoiseRatio);
for ii = 1:nOfsignalToNoiseRatio
    mean_capacity(ii) = mean(CapacityVec(ii,:));
end

% figure
% x = (1:nOfsignalToNoiseRatio);
% y = CapacityVec(x,1);
% plot((1:nOfsignalToNoiseRatio), y ,'b'); hold on;
% 
% % plot(snrdB, 1/5,'b'); hold on;
% hold off

figure 
indices = find(signalToNoiseRatio<snrTh);
signalToNoiseRatio(indices) = [];
curve1 = 1./signalToNoiseRatio;


plot(signalToNoiseRatio, curve1 ,'r'); hold on;
curve2 =  1/snrTh*ones(1,length(signalToNoiseRatio));
plot(signalToNoiseRatio, curve2, 'b'); hold on;
x2 = [signalToNoiseRatio, fliplr(signalToNoiseRatio)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'g');
axis([1 8 0 2])
legend('1/SNR_{o}','1/SNR');
title('1/SNR vs SNR');
xlabel('SNR');ylabel('1/SNR');

% truncated capacity
%snr                         = signalToNoiseRatio;
snr_trunc_index             = find(snr <= snrTh);
snr_trunc                   = snr;
snr_trunc(snr_trunc_index)  = [];
p_sim_trunc                 = p_sim; 
p_sim_trunc(snr_trunc_index)= [];
inv_snr_expectation_trunc = sum((p_sim_trunc)./snr_trunc);
C_inv_trunc = (log2(1+ 1/inv_snr_expectation));


% water filling algo
function P = waterFilling(SNR,Pmax)
    % initial power allocation
    initP = (Pmax + sum(1./SNR)) ./ ( length(SNR) ) - 1./SNR;

    % waterfilling algorithm
    while any( initP < 0 )
        negIndex        = initP <= 0;
        posIndex        = initP >  0;
        NkRem           = nnz(posIndex); % # of non-zero elements in posIndex 
        SNRRem          = SNR(posIndex); 
        powAllcTemp     = (Pmax + sum(1./SNRRem)) ./ (NkRem) - 1./SNRRem;
        initP(negIndex) = 0;
        initP(posIndex) = powAllcTemp;
    end
    P              = initP;
end

function [h]=rayleighFading(M,N,fd,Ts)
% function to generate Rayleigh Fading samples based on Clarke's model
% M = number of multi-paths in the channel
% N = number of samples to generate
% fd = maximum Doppler frequency
% Ts = sampling period
a=0;
h_re = zeros(1,N);
h_im = zeros(1, N);
b=2*pi;
alpha=a+(b-a)*rand(1,M); %uniformly distributed from 0 to 2 pi
beta=a+(b-a)*rand(1,M);  %uniformly distributed from 0 to 2 pi
theta=a+(b-a)*rand(1,M); %uniformly distributed from 0 to 2 pi
m=1:M;
for n=1:N
    x=cos(((2.*m-1)*pi+theta)/(4*M));
    h_re(n)=1/sqrt(M)*sum(cos(2*pi*fd*x*n'*Ts+alpha));
    h_im(n)=1/sqrt(M)*sum(sin(2*pi*fd*x*n'*Ts+beta));
end
h=h_re+1i*h_im;
end

