X = 5; 
K = 1;
a = 0.2 + 0.025 * X * 5;
b = 100 + 0.5 * X * 5;
yfinal = (K*7570*a*b)/(K*7570*a*b - (62.61)^2); % by final value theorem  
sim('IE474_O1_Task1_simulink_X.mdl')
S = stepinfo(yout,tout,yfinal);
figure, plot(tout,yout)
grid
disp(S.Overshoot);
disp(S.SettlingTime);