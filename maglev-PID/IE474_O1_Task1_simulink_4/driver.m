X = 4;                                % group 4
K = 1;
a = 0.2 + 0.025 * X * 5;
b = 100 + 0.5 * X * 5;
yfinal = 1;                           % by final value theorem  
sim('IE474_O1_Task1_simulink_4.mdl')
S = stepinfo(yout,tout,yfinal);
figure, plot(tout, yout)
grid
disp(S.Overshoot);
disp(S.SettlingTime);