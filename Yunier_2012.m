function [yk,time,error,index_iteration_outer,G_yk ] =...
    Yunier_2012(C,M,e,Max_iteration,x0, Slater_point,type_problem)
%YUNIER_2012 
%   
%ECM_VIP 此处显示有关此函数的摘要
%   此处显示详细说明

% beta_k = 1;
% eta_k = 1;

%% Inputs
%
% 


%% variables for debug
G_yk = [];
Error = [];
S = [];
Yk = [];


%%


t1 = cputime;

[n,m]=size(C);
[n,n] = size(C{1,1});

zk = x0;
xk=zeros(n,1);

Max_iteration_projection = 1000;

index_iteration_outer = 0;
error = 1;

p = 1;
theta =2;
sigma_k = 0;

% zk = ACRM_oneStep(C,n,m,x0);
xk = x0;

while (1) % outer loop is to update x^k =\sum_{i=1}^{k-1} y^i
    

    beta_k = 1/(index_iteration_outer+1)^p;
    
    jk = Find_max_function_value(C,n,m,zk);
    if zk'*C{1,jk}*zk +2*C{2,jk}'*zk - C{3,jk} <=0
        yk=zk;
    else
        yk = zk;
        index_iteration_inner = 0;
        while (1) % inner loop to update y^k            
            for i=1:m
                F(i) = xk'*C{1,i}*xk +2*C{2,i}'*xk - C{3,i};
            end
            [value_max,index_max]=max(F);
            yk  = Approximate_Projection(C{1,index_max},C{2,index_max},C{3,index_max},yk);
            
%             yk = ACRM_oneStep(C,n,m,yk);
            jk = Find_max_function_value(C,n,m,yk);

            g_yk = yk'*C{1,jk}*yk +2*C{2,jk}'*yk - C{3,jk};
            G_yk = [G_yk g_yk];

            g_Slater_point = Slater_point'*C{1,jk}*Slater_point +2*C{2,jk}'*Slater_point - C{3,jk};
            f = (g_yk* norm(yk - Slater_point))/(g_yk - g_Slater_point);

            if f<=theta*beta_k
                break;
            end
            
            index_iteration_inner = index_iteration_inner + 1;
            if index_iteration_inner > 10
                break;
            end
            
        end
    end
    eta_k = max(1,norm(F_monotone(M,yk,type_problem)));
    for i=1:m
        F(i) = xk'*C{1,i}*xk +2*C{2,i}'*xk - C{3,i};
    end
    [value_max,index_max]=max(F);
    zk  = Approximate_Projection(C{1,index_max},C{2,index_max},C{3,index_max},yk - beta_k/eta_k*F_monotone(M,yk,type_problem));
%     zk = ACRM_oneStep(C,n,m,yk - beta_k/eta_k*F_monotone(M,yk,type_problem));
    
    error = norm (zk -yk); % stopping criteria
    Error = [Error;error];
    if error <e 
        break;
    end
    
    sigma_k = sigma_k + beta_k/eta_k;
    xk = (1 - beta_k/(eta_k*sigma_k))*xk + beta_k/(eta_k*sigma_k)*yk;
    S = [S;(1 - beta_k/(eta_k*sigma_k))];
    
    
    index_iteration_outer = index_iteration_outer + 1;
    
    if index_iteration_outer > Max_iteration
        break;  
    end
         
end



t2 = cputime;
time = t2 - t1;

error = norm(yk-Dykstra(C,n,m,e,Max_iteration_projection,yk-...
    0.1*F_monotone(M,yk,type_problem)));

end




