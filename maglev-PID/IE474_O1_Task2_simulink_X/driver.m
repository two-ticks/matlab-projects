close all;
clc
X = 5;
a_vec = 0.2 + 0.025 * X * (1:10); % vector including all possible a values 
b_vec = 100 + 0.5 * X * (1:10);   % vector including all possible a values 
yfinal = 1.0056;                  % final value that you have calculated in (a) using the final value theorem. 
K_vec = (0:.01:1);
for ii=1:length(K_vec)
 for jj=1:length(K_vec)
 a = a_vec(ii);                    %
 b = b_vec(ii);                    %
 K = K_vec(ii);
 sim('IE474_O1_Task2_simulink_X.mdl') % You simulate vector including the
 S(ii,jj) = stepinfo(yout, tout, yfinal);
 PercentOveshoot_matrix(ii,jj) = S(ii,jj).Overshoot;
 SettlingTime_matrix(ii,jj) = S(ii,jj).SettlingTime;
 
 end
end