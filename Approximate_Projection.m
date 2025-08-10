function x = Approximate_Projection(A,b,alpha,x0)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

%% Intitialzation
% y = x0; 
f = x0'*A*x0 + 2*x0'*b -alpha;
g= 2*A*x0+2*b;
x = x0 -(max(0,f))/(norm(g)^2)*g;



end

