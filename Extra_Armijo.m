function [ xk,time,error,index_iteration] = Extra_Armijo(C,M,e,Max_iteration,x0,type_problem)
% function [ xk,time,error,index_iteration] = Extragradient(C,M,e,Max_iteration,x0,type_problem)
%% EXTR 
%   This .m is to ultilize the Extragradient method for the variational
%   inequality with monotone operator

t1 = cputime; 

error = 1;

[n,m]=size(C);
[n,n] = size(C{1,1});
x_temp = zeros(n,1);
lambda = 0.05;


beta_k = 0.1;
delta = 0.5; % delta \in (0,1)

% Max_iteration = 50;
Max_iteration_projection = 1000;

xk=Dykstra(C,n,m,e,Max_iteration_projection,x0);

index_iteration = 0;

while (1)
    
    zk = xk - beta_k*F_monotone(M,xk,type_problem);
    
    P_C_zk = Dykstra(C,n,m,e,Max_iteration_projection,zk);
    
    %% Stopping criteria
    if norm(xk-P_C_zk)<e
        break;
    end
    
    if index_iteration > Max_iteration
        break;
    end
    
    %% Line search
    j=1;
    while (1) % linesearch
        a = 2^(-j);
        if (F_monotone(M,a*P_C_zk +(1-a)*xk,type_problem))'*(xk-P_C_zk)>=(delta/beta_k)*norm(xk-P_C_zk)^2
            break;
        end

        j = j +1;
    end
    %% main loop: Extragradient method with line search 
    alpha_k = 2^(-j);
    
    yk = alpha_k*P_C_zk + (1-alpha_k)*xk;
    
    gamma_k = (F_monotone(M,yk,type_problem))'*(xk-yk)/norm(F_monotone(M,yk,type_problem))^2;
    
    wk = xk-gamma_k*F_monotone(M,yk,type_problem);
    
    xk = Dykstra(C,n,m,e,Max_iteration_projection,wk);
    
    index_iteration = index_iteration + 1;
    
end

error = norm(xk-Dykstra(C,n,m,e,Max_iteration_projection,xk-...
    0.1*F_monotone(M,xk,type_problem)));

t2= cputime;
time = t2-t1;

end



