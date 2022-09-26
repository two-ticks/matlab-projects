% Question A
N = 1000;
e7 = 0;
e127 = 0;
e511 = 0;
p_vec = 0.05:0.05:0.4;
n_vec = [7 127 511];

e = zeros(length(p_vec), length(n_vec));

n_ptr = 0;
for n = [7 127 511]
n_ptr = n_ptr + 1; 
p_ptr = 0;
    for p = 0.05:0.05:0.4
        p_ptr = p_ptr + 1;
        G = ones(1,n); % [I : P] 
            for itr = 1: N
                x = randi([0 1]);
                c = x*G;
                pd = makedist('Binomial','N', 1,'p', p);
                z = random(pd, 1, n);
                % s = c + n
                s = mod(c + z, 2); 
                % MAP  
                if ((x == 1) && (sum(s(:)) < 0.5*n)) 
                    e(p_ptr, n_ptr) = e(p_ptr, n_ptr) + 1;
                elseif ((x == 0) && (sum(s(:)) > 0.5*n)) 
                    e(p_ptr, n_ptr) = e(p_ptr,n_ptr) + 1;
                end
            end
    
    end
end

BLER = e./N;

% theoretical 

% Question B 

% probability of error in repetition code increases with flip probability
% of BSC

% pe = ( ([n/2]) choose (i=0) p^(n-i)(1-p)^i ) 

% commented for debugging
pe_matrix = zeros(length(p_vec), length(n_vec));

for nn = 1: length(n_vec)
    nnn = n_vec(nn);
    for pp = 1:length(p_vec)
        pe = 0;
            for ii = 0:floor(nnn/2)
                 % pe = ( ([n/2]) choose (i=0) p^(n-i)(1-p)^i )
                 pe = pe + nchoosek(sym(nnn),ii)*((p_vec(pp)^(nnn-ii))*((1-(p_vec(pp)))^(ii))); 
            end
        pe_matrix(pp,nn) = pe;
    end
end


figure 

plot(p_vec, BLER(:,1), Color=	"#0072BD"); hold on
plot(p_vec, BLER(:,2), Color=	"#D95319"); hold on
plot(p_vec, BLER(:,3), Color=	"#EDB120"); hold on 

plot(p_vec, pe_matrix(:,1),  LineStyle = "--", Color=	"#7E2F8E"); hold on
plot(p_vec, pe_matrix(:,2), LineStyle = "--", Color=	"#77AC30"); hold on
plot(p_vec, pe_matrix(:,3), LineStyle = "--",  Color=	"#A2142F"); 

title('BLER vs p')
legend('simulated e7','simulated  e127','simulated  e511','theoretical e7','theoretical e127','theoretical e511','Location','northwest')