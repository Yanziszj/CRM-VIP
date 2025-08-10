function [ F_xk ] = F_monotone(M,xk, type_problem)
%% F_MONOTONE 
%   to evaluate the F at point xk


%% Inputs
%  n: dimension
%  m: number of sets
%  A: the matrix
%  B:
%  C: 
%  c:
%  xk:
%  j:

%% Outputs
% F_xk

%% Initialization

[n,m] = size(M);

F_xk = zeros(n,1);   % 
c = zeros(n,1) + 20;
alpha = zeros(n,1)+3;
alpha  = alpha./sum(alpha); % alpha is the convex combinaition coefficients

switch type_problem
%% paramonote which could be a subgradient of a convex function    
%  F(x) = Mx + \alppha.*x^3
%  M: sysmetric and semipositive definitie matrix.

    case 1
        F_xk = M*xk;
        for i=1:n
            F_xk(i) = F_xk(i) + alpha(i)*xk(i)^3 + c(i);
        end
        
%% paramonote which is not any subgradient of a convex function      
%  F(x)  = Bx
    case 2  % monotone
        F_xk = M*xk + c;
        
%%  Linear operator which is monotone but not paramonotone      
    case 3  
        F_xk = M*xk + c;

end

