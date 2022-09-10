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
figure;
plot((0:N-1)*Ts, initialSignal);
title('transmitted signal: s(t)');

h = rayleighFading(M, N, fd, Ts); % transfer function of Rayleigh Channel
h_re = real(h); h_im = imag(h);

% output signal from channel
filteredSignal = ifft(fft(initialSignal).*fft(h));
figure;
plot((0:(N-1))*Ts, abs(filteredSignal));
title('recieved signal: |s(t) \ast h(t)|');
xlabel('time(s)'); ylabel('amplitude');


% find: 1) the number of fades in the initial signal; and 2) where these fades occur
% find the number of fades in the initial signal, by first quantising the signal around the mean
quantised = abs(filteredSignal) / mean(abs(filteredSignal)) < 1;
crossingsUp  = find(diff(quantised) == -1);      % index of when signal crosses above threshold of the mean, suggesting the end of a fade
crossingsDown  = find(diff(quantised) == 1);     % index of when signal crosses below threshold of the mean, suggesting the beginning of a fade

% remove stray crossings at start & end of data
if(length(crossingsUp) > 1 && length(crossingsDown) > 1)
    if(crossingsUp(1) < crossingsDown(1))
        crossingsUp(1) = [];
    end
    if(crossingsDown(end) > crossingsUp(end))
        crossingsDown(end) = [];
    end
end

% total fades 
fprintf('total fades: %d \n', length(crossingsUp));

% fade duration
differenceCrossing = crossingsUp - crossingsDown;
fprintf('%f sec\n', mean(differenceCrossing)/10000);
figure;
plot((0:length(differenceCrossing)-1)*Ts, differenceCrossing);
title('fade duration');
hold on
yline(mean(differenceCrossing), 'r');
legend('fade duration','average fade duration');
hold off

% autocorrelation
acf = abs(xcorr(filteredSignal, filteredSignal));
figure;
plot((0:length(acf)-1)*Ts, acf);
title('Auto Correlation Function');
figure;

% comparing the PDF of overal response of the channel against the PDF of Rayleigh distribution
[val,bin] = hist(abs(h),1000); % pdf of generated Raleigh Fading samples
plot(bin, val/trapz(bin,val)); % normalizing the PDF to match theoretical result
% trapz function gives the total area under the PDF curve. It is used as the normalization factor

hold on;
z=0:0.1:3; sigma=1;
y=2*z/(sigma^2).*exp(-z.^2/(sigma^2)); % theoretical Raleigh pdf
plot(z,y,'r');
title('Probability Density Function');
legend('Simulated pdf','Theoretical Rayleigh pdf');

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

