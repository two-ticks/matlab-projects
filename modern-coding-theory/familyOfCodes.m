% clear
clc;
clear;

% parameters
N = 1000; % number of trials
e7 = 0; e127 = 0; e511 = 0;

p_vec = 0.05:0.05:0.4; % p vector
m_vec = [3 5 7];       % m vector

e = zeros(length(p_vec), length(m_vec)); % matrix containing number of errors 

n_ptr = 0;  % intialize index of n = 0 
for m = m_vec
    n_ptr = n_ptr + 1; % index of n vector
    G = [ones(m) randi([0, 1], [m,2^m-m-1])]; % [I : P]
    M = (dec2bin(0:2^m-1)' - '0')';
    C = mod(M*G, 2); % 2^m codewords of 'm' length in binary field
    C_red =  sum(C,2); % Hamming weight of each codeword 
    
        p_ptr = 0; % intialize index of p = 0 
        for p = 0.05:0.05:0.4
            p_ptr = p_ptr + 1; % index of p vector
                for itr = 1: N
                    x = randi([0 1], [1, m]);
                    c = mod(x*G,2); % codeword 
                    pd = makedist('Binomial','N', 1, 'p', p);
                    z = random(pd, 1, 2^m-1); % generated noise with given bit flip probability distribution 
                    r = mod(c + z, 2); % r = c + m, received vector 
                    
                    % MAP decoding  
                    % summary: 
                    % step 1: finding hamming distance
                    % step 2: comparing recieved vector with codeword
                    % step 3: adding 1 to # of error if r is not equal to estimated
                    
                    [min_value, codeword_index] = min(sum(repmat(r,2^m,1).*C)); % hamming distance matching
    
                    if (C(codeword_index)~=c)
                        e(p_ptr, n_ptr) = e(p_ptr, n_ptr) + 1;
                    end
                end
        end
end

BLER = e./N;

e7 = BLER(:,1);
e127 = BLER(:,2);
e511 = BLER(:,3);

figure
plot(p_vec, BLER(:,1), Color=	"#0072BD"); hold on
plot(p_vec, BLER(:,2), Color=	"#D95319"); hold on
plot(p_vec, BLER(:,3), Color=	"#EDB120"); hold on 

title('BLER vs p')
legend('simulated m1','simulated  m2','simulated  m3','Location','northwest')

