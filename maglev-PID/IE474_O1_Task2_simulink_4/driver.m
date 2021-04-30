close all;                         % k > 0.0048
clc;
X      = 4;
a_vec  = 0.2 + 0.025 * X * (1:10);   % vector including all possible a values 
b_vec  = 100 + 0.5 * X * (1:10);     % vector including all possible a values 
yfinal = 0.1;                        % final value that you have calculated in (a) using the final value theorem. 
K_vec  = 1+0.0048*(1:10);      
for ii=1:length(K_vec)
    for jj=1:length(K_vec)
     a = a_vec(jj);                  
     b = b_vec(ii);
     K = K_vec(jj);
     sim('IE474_O1_Task2_simulink_4.mdl') % You simulate vector including the
     S = stepinfo(yout, tout, yfinal);
     PercentOveshoot_matrix(ii,jj) = S.Overshoot;
     SettlingTime_matrix(ii,jj) = S.SettlingTime;
    end
end

% find max overshoot and minimum settling time
figure
%mesh(a_vec, b_vec,PercentOveshoot_matrix(1:10,1:10));
mesh(a_vec, b_vec,SettlingTime_matrix(1:10,1:10));

