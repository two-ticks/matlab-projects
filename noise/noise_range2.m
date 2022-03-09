%% write a code for white noise (1khz-16khz) and one narrow band sound (1khz-2khz)

low_f1 = 1000;     % lower limit of the fq range
high_f1 = 16000;   % upper limit of the fq range
fs = 33000;        % sampling frequency
time = 15;         % in seconds
len  = 15;         % in seconds

%% geberate a gaussian white signal
n = round(time*fs);  % number of samples
x = randn(n ,1);

%% band pass the white noise to get the narrow band 
y = bandpass(x, [low_f1 high_f1],fs);
[p, f] = pwelch(y, 1024, 768, 1024, fs);

%% read/generate sound
f_sig1 = 1000; % lower limit of the fq range
f_sig2 = 2000; % upper limit of the fq range
t   = linspace(0, len, fs*time);              % time Vector
signal = sin(2*pi*f_sig1*t) + t.*cos(2*pi*f_sig2*t) + 2*sin(2*pi*f_sig1*1.5*t); 

% reading audio from sound_input.mp3
% [signal,fsig] = audioread('sound_input.mp3'); % reading signal

%% band pass the sound to get the narrow band 
sig_bp = bandpass(signal, [f_sig1 f_sig2],fs);
[sig_nb, f_sig_nb] = pwelch(sig_bp, 1024, 768, 1024, fs);

%% plotting PSD
subplot(2,1,1);
plot(f, 10*log10(p))

subplot(2,1,2);
plot(f_sig_nb, 10*log10(sig_nb))

%% write audio file - output.wav
signal_final =  y' + sig_bp(1:size(y,1));
sound(signal_final,fs)
audiowrite('output.wav',signal_final,fs)

