N = 1000;
e7 = 0;
e127 = 0;
e511 = 0;
p_vec = 0.05:0.05:0.4;
n_vec = [7 13];

e = zeros(length(p_vec), length(n_vec));

n_ptr = 0;
for n = [7 13]
n_ptr = n_ptr + 1; 
p_ptr = 0;
    for p = 0.05:0.05:0.4
        p_ptr = p_ptr + 1;
        G = [ones(n) randi([0, 1], [n,2^n-n-1])]; % [I : P] 
            for itr = 1: 1000
                x = randi([0 1], [1, n]);
                c = x*G;
                pd = makedist('Binomial','N', 1, 'p', p);
                z = random(pd, 1, 2^n-1);
                % s = c + n
                s = mod(c + z, 2); 
                % MAP 
                if ((x~=s)) 
                    e(p_ptr, n_ptr) = e(p_ptr, n_ptr) + 1;
                end
            end
    
    end
end

BLER = e./N;


figure
plot(p_vec, BLER(:,1)); hold on
plot(p_vec, BLER(:,2)); hold on
plot(p_vec, BLER(:,3)); 




