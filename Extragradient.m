function [ xk,time,error,index_iteration] = Extragradient(C,M,e,Max_iteration,x0,type_problem)
%% EXTR 
%   This .m is to ultilize the Extragradient method for the variational
%   inequality with monotone operator

t1 = cputime; 

error = 1;

[n,m]=size(C);
[n,n] = size(C{1,1});
x_temp = zeros(n,1);
lambda = 0.05;
xk=x0;


% Max_iteration = 50;
Max_iteration_projection = 1000;

index_iteration = 0;

while (1)

%     yk = Projection_ellipsoids_cvx(C,n,m,xk -lambda*F_monotone(M,xk,type_problem));
%     
%     xk = Projection_ellipsoids_cvx(C,n,m,xk -lambda*F_monotone(M,yk,type_problem));
    
%     
    yk = Dykstra(C,n,m,e,Max_iteration_projection,xk -lambda*F_monotone(M,xk,type_problem));
    
    xk = Dykstra(C,n,m,e,Max_iteration_projection,xk -lambda*F_monotone(M,yk,type_problem));
    
    index_iteration = index_iteration + 1;
    
    
    %% Stopping criteria
    error = norm(xk-yk);
    if error <= e
        break;
    end
    
    if index_iteration > Max_iteration
        break;
    end
    
end

error = norm(xk-Dykstra(C,n,m,e,Max_iteration_projection,xk-...
    0.1*F_monotone(M,xk,type_problem)));

t2= cputime;
time = t2-t1;

end

