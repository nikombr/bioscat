clear; close all; clc; 

d = linspace(0,3,300);

kSE = @(tau, ell, d) tau^2*exp(-d.^2/(2*ell^2));
% 
% tau = [0.25, 0.5, 1, 10]
% 
% for t = tau
%     plot(d,kSE(t,1,d))
%     hold on
% end

ps = [0 1 3]

for p = ps
    plot(d,kM(p,1,d),'label',sprintf('$p=%d$',p))
    hold on
end




function res = kM(p, ell, d)

nu = p + 1/2;

sum = 0;

for i = 0:p
    sum = sum + factorial(p+i)/(factorial(i)*factorial(p-i))*(d*sqrt(8*nu)/ell).^(p-i);
end

res = exp(-d*sqrt(2*nu)/ell)*factorial(p)/factorial(2*p).*sum;




end