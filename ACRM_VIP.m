function [xk,time,error2,index_iteration,index_stopping] =...
    ACRM_VIP( C,M,e,Max_iteration,x0,type_problem)
%% ACRM_VIP
%   This file to solve the variation inequality with paramonotone operators
%   Feasible sets: intersections of ellipsoids 
%       E: = \{x \in R^n: x' A x+2b' x \leq \alpha \}
%

%% Inupts
%   C is a 3*n cell to store the information of all the ellipsoids
%       C{i,1}: the matrix A of i-th ellipsoid
%       C{i,2}: the vector b of i-th ellipsoid
%       C{i,3}: the constant alpha of i-th ellipsoid
%       e: error tolerance
%       Max_iteration: max number of iteration
%       x0: initial point
%       type_problem: 1.
%                               2.
%                               3.
%% Outputs



%% Initailizations
xk = x0;
beta_k = 1;   % should be summable
eta_k = 1;

e2= 10^(-5);

[n,m]=size(C);   % n is the dimension of the problem
[n,n] = size(C{1,1});

x_temp = zeros(n,1);
Max_iteration_projection = 500;

index_stopping = 'xk and xk+1 are too close';
index_iteration = 1;
p = 1;

t1 = cputime;
x_temp = 0;

%% Debug variables
Error = []; % store the whole error
X=[]; % to store the history information of xk

% while (norm(xk-x_temp)>e)

%% Iteration
while (1)
    x_temp = xk;
    
    beta_k = 1/(index_iteration^p);  % summable1
    F_xk = F_monotone(M,xk,type_problem);
    eta_k = max(1,norm(F_xk));
%     xk =xk- beta_k/eta_k*F_xk;
    % choose the type of epsilon
%     fprintf('beta/eta:%f\n',beta_k/eta_k);
    xk = ACRM_oneStep(C,n,m,xk- beta_k/eta_k*F_xk);
%     fprintf('Feasibility of first ellipsoid:%f\n',xk'*C{1,1}*xk +2*C{2,1}'*xk -C{3,1});
%     fprintf('Feasibility of second ellipsoid:%f\n',xk'*C{1,2}*xk +2*C{2,2}'*xk -C{3,2});
    error = norm(xk-x_temp);
%     error = norm(Projection(C{1,1},C{2,1},C{3,1},xk- beta_k/eta_k*F_xk)-xk);
%     error = norm(xk-Dykstra(C,n,m,e,Max_iteration_projection,xk -0.1*F_monotone(M,xk,type_problem)));
%     fprintf('Error is:%f\n',error);
    Error = [Error;error];
    if error <=e
        break;
    end
    
    ACRM_oneStep(C,n,m,x_temp- beta_k/eta_k*F_xk);
    
    % the max-iteration check
    if index_iteration>Max_iteration
        
        break;    % stop the cyclic loop
    end

    index_iteration = index_iteration+1;
%             X = [X xk];           

end

% These steps are to check the stopping criteria
% beta_k = 1/(index_iteration^2);  % summable1
% F_xk = F_monotone(M,xk,type_problem);
% eta_k = max(1,norm(F_xk));


% if norm(Projection_ellipsoids_cvx(C, n, m, xk- beta_k/eta_k*F_xk))<=e2 && index_iteration <= Max_iteration
%     %  to indicate that program stops since the interation exceed the max
%     index_stopping = 'succesful';
% elseif index_iteration > Max_iteration
%     index_stopping =  'the iteration exceed the max';
% else
%     index_stopping =  'converges to a point which is not solution';
% end
error1 = error;
t2 = cputime;
time = t2-t1;

error2 = norm(xk-Dykstra(C,n,m,e,Max_iteration_projection,xk - 0.1*F_monotone...
    (M,xk,type_problem)));

end



