% Find: 1) the number of fades in the initial signal; and 2) where these
% fades occur.
%--------------------------------------------------------------------------

clc;
clear;

% doppler fade

freq = 5*10^9; %Hz
v_rx = 50;     % speed of reciever 
f_doppler = (v_rx * freq) / (3*10^8);  %doppler spread 

% Find the number of fades in the initial signal, by first quantising the
% signal around the mean

M=15; %number of multipaths
N=10^5; %number of samples to generate
fd= f_doppler; % Maximum doppler spread in hertz
Ts=0.0001; % Sampling period in seconds
x=(0:N-1)*Ts;
initial_signal = 12*sin(5*x+2)+20*cos(5*10^3*x);
figure;
plot((0:N-1)*Ts, initial_signal);
title('initial signal');

h=rayleighFading(M,N,fd,Ts);

h_re=real(h);
h_im=imag(h);

% output signal from channel
filtered_signal = ifft(fft(initial_signal).*fft(h));

figure;
% subplot(2,1,1);
plot((0:(N-1))*Ts, abs(filtered_signal));
title('abs(filtered_signal)');
xlabel('time(s)');ylabel('Amplitude |hI(t)|');


quantised = abs(filtered_signal) / mean(abs(filtered_signal)) < 1;
crossings_up  = find(diff(quantised) == -1); % Index of when signal crosses above threshold of the mean, suggesting the end of a fade
crossings_down  = find(diff(quantised) == 1); % Index of when signal crosses below threshold of the mean, suggesting the beginning of a fade

% Remove stray crossings at start & end of data
if(length(crossings_up) > 1 && length(crossings_down) > 1)
    if(crossings_up(1) < crossings_down(1))
        crossings_up(1) = [];
    end
    if(crossings_down(end) > crossings_up(end))
        crossings_down(end) = [];
    end
end

% Hence the number of fades in the initial signal is...

% plot(x, initial_signal);

nfades = length(crossings_up);

% fade duration

dif_crossing = crossings_up - crossings_down;
mean_dif_crossing = mean(dif_crossing);
fprintf('%f sec\n',mean_dif_crossing/10000);
figure;
% subplot(2,1,1);
plot((0:length(dif_crossing)-1)*Ts, dif_crossing);
title('fade duration');


% autocorrelation

acf = abs(xcorr(filtered_signal, filtered_signal));

figure;
plot((0:length(acf)-1)*Ts, acf);
title('auto correlation function');
figure;
%comparing the PDF of overal response of the channel against the PDF of Rayleigh distribution
[val,bin]=hist(abs(h),1000); % pdf of generated Raleigh Fading samples
plot(bin,val/trapz(bin,val)); %Normalizing the PDF to match theoretical result
%Trapz function gives the total area under the PDF curve. It is used as the normalization factor
hold on;
z=0:0.1:3; sigma=1;
y=2*z/(sigma^2).*exp(-z.^2/(sigma^2)); % theoretical Raleigh pdf
plot(z,y,'r');
title('Probability density function');
legend('Simulated pdf','Theoretical Rayleigh pdf');

