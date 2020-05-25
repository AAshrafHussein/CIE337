%{
Communication Theory and Systems
CIE 337 - Spring 20 20
By: Ahmed Ashraf Hussein
Email: s-ahmedashraf@zewailcity.edu.eg
Final Assessment
PART2_THEORETICAL_ERROR_PROBABILITY
%}
%{
Prompt:
It's required to plot the average error probability vs Eb/N0 for
p = 0.25, 0.5 and 0.75 all on the same figure.
Note that ? is a function of p
%}

clear;
clc;
T_b = 1;
A = .1:.1:20;
E_b = A.^2*T_b;
N_o = 25;
dB = 10*log10(E_b/N_o);
p = [0.25, 0.5, 0.75];

for j = 1:length(p)
    pe_all = [];
    for i = 1:length(A)
        pe = (1-p(j))*qfunc((((N_o/(16*A(i)*T_b))*log((1-p(j))/p(j)))+A(i))*sqrt(2)/sqrt(N_o/T_b))+(p(j))*qfunc((-1*((N_o/(16*A(i)*T_b))*log((1-p(j))/p(j)))+A(i))*sqrt(2)/sqrt(N_o/T_b));
        pe_all = [pe_all pe];
    end
    if p(j) == 0.25
        semilogy(dB, pe_all, 'o');
        hold on;
        continue;
    end        
    semilogy(dB, pe_all);
    hold on;
    
end
axis([0 15 10^-10 10^-2])
legend('p=0.25', 'p=0.5', 'p=0.75')
xlabel('E_b/N_o, dB')
ylabel('Probability of Error, P_e')
title('Average Error Probabilities vs E_b/N_o for different transmission probabilities')