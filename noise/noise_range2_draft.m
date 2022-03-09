%% write a code for white noise (1khz-16khz) and one narrow band sound (1khz-2khz)

%cf1 = 8500;         % central fq
% bw = 0.25;            % bandwidth in octave, log2(f2/f1) 
low_f1 = 1000; %cf1 / 2 ^ (bw/2);      % lower limit of the fq range
high_f1 = 16000; %cf1 * 2 ^ (bw/2);     % upper limit of the fq range
fs = 33000;          % sampling frequency
time = 2;           % in seconds



% geberate a gaussian white signal
n = round(time*fs);  % number of samples
x = randn(n ,1);
% band pass the white noise to get the narrow band 
y = bandpass(x, [low_f1 high_f1],fs);
[p, f] = pwelch(y, 1024, 768, 1024, fs);

% read/genrate sound
[signal,fsig] = audioread('sound_input.mp3');
% f_sig1 = 1000;
% f_sig2 = 2000;
% t   = linspace(0, time, fs*len);                 % Time Vector
% signal = sin(2*pi*f_sig1*t); 

% band pass the sound to get the narrow band 
sig_bp = bandpass(signal, [f_sig1 f_sig2],fs);
[sig_nb, f_sig_nb] = pwelch(sig_bp, 1024, 768, 1024, fsig);
% sound(y, fs) 
% sound(signal, fs)
figure;
subplot(2,1,1);
plot(f, 10*log10(p))

subplot(2,1,2);
figure(2)
plot(f_sig_nb, 10*log10(sig_nb))

% write audio file - output.wav
% signal_final = sig_bp + y;
audiowrite('output.wav',signal_final,fs)