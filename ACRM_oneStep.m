function [ x ] = ACRM_oneStep(C,n,m,xk)
%% ACRM_oneStep.m 
%  This file is to implemention of one step CRM with approximat projection
%   Reference: R., Behling, etc. (2022). 
%                     Circumcentering approximate reflections for solving the convex feasibility problem. 
%                     Fixed Point Theory and Algorithms for Sciences and Engineering, 2022(1), 1.
% 
%
%

wk=zeros(n,1);
s=zeros(n,1);
e =sqrt(eps);

for i=1:m  %the cyclic loop
    f = xk'*C{1,i}*xk + 2*xk'*C{2,i} -C{3,i};
    g= 2*C{1,i}*xk+2*C{2,i};
%     x = xk -(max(0,f))/(norm(g)^2)*g;
%     vik = xk - Approximate_Projection(C{1,i},C{2,i},C{3,i},xk);    
    vik = (max(0,f))/(norm(g)^2)*g;    
    wk =wk +vik;
    s(i) = norm(vik)^2;
end
wk =wk/m;

% calculate the alpha_k
if norm(wk)<=e
    alpha_k = 0;
else
    alpha_k = (sum(s))/(m*norm(wk)^2);
end
xk= xk -alpha_k*wk;

x = xk;
end

