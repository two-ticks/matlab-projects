f1=100;
f2=300;
time=2;                   % length of time
Ts=1/1000;                % time interval between samples
x=randn(1,time/Ts);              % generate noise signal
freqs=[ 0.0  0.2  0.3  0.5  0.6  1];
amps=[0 0 1 1 0 0];
b=firpm(100,freqs,amps);         % BP filter
ylp=filter(b,1,x);               % do the filtering
% figure(1),plotspec(ylp,Ts)       % plot the output spectrum
figure(1)
freqz(b,1,1024,1/Ts)        % Filter Bode Plot