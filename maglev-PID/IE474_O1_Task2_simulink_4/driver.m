clear all;
close all;                           % k > 0.0048 
clc;
X      = 4;
a_vec  = 0.2 + 0.025 * X * (1:10);   % vector including all possible a values 
b_vec  = 100 + 0.5 * X * (1:10);     % vector including all possible a values 
yfinal = 0.1;                        % final value that you have calculated in (a) using the final value theorem. 
K_vec  = ones(1:10);      
for ii=1:length(K_vec)
    for jj=1:length(K_vec)
         a = a_vec(jj);                  
         b = b_vec(ii);
         K = K_vec(jj);
         sim('IE474_O1_Task2_simulink_4.mdl') % simulate
         S = stepinfo(yout, tout, yfinal);
         PercentOveshoot_matrix(ii,jj) = S.Overshoot;
         SettlingTime_matrix(ii,jj) = S.SettlingTime;
    end
end

% plotting
subplot(2,1,1);
mesh(a_vec, b_vec,PercentOveshoot_matrix);
title('Percent Oveshoot');

subplot(2,1,2); 
mesh(a_vec, b_vec,SettlingTime_matrix);
title('Settling Time');

% find max overshoot and minimum settling time
max_matrix = (PercentOveshoot_matrix-SettlingTime_matrix); % finding element having least settling time and maximum percentage oveshoot
[row, col] = find(max_matrix == max(max_matrix(:)));

fprintf('Optimum a value = %d\n', a_vec(col));
fprintf('Optimum b value = %d\n', b_vec(row));
