% UNIT IMPULSE FUNCTION
n=-20:20;
delta_n = n==0;
stem(n,delta_n);
xlabel('Time Samples');
ylabel('Amplitude');