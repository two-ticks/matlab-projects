clc;
clear all;
close all;

X = 4;                                % group 4
K = 1;
a = 0.2 + 0.025 * X * 5;
b = 100 + 0.5 * X * 5;
yfinal = 1;                           % by final value theorem  

sim('IE474_O1_Task1_simulink_4.mdl')  % simulate
S = stepinfo(yout,tout,yfinal);       % step response characteristics

% plotting
figure, plot(tout, yout)
grid
axis([-0.01, 0.09, 1.002, 1.014]);

fprintf("Percentage Overshoot %d \n",S.Overshoot);
fprintf("Settling Time %d \n",S.SettlingTime);

% Routh Hurwitz
% s^3 + 7570K s^2 + (837999K - 3920.0121)s + 582890K = 0 
% for stability Inner Product > Outer Product in third order CE
% therefore,  7570K(837999K - 3920.0121) > 582890 
% => K > 0.0048 for stable system
