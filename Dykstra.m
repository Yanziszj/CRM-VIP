function [x,time,index_iteration,error]  = Dykstra(C,n,m,e,Max_iteration,x0)
%DYKSTRA 
%%  This Dykstra.m is to uitilize the Dktras's method for the best
%            approximation with m sets. The Algorithm is in the form of the
%            algorithm in Birgin, E. G., & Raydan, M. (2005). 
%            Robust stopping criteria for Dykstra's algorithm. 
%            SIAM Journal on Scientific Computing, 26(4), 1405-1414.
%            Author: Di Liu (di.liu@impa.br)
%
%

%% Stopping Criteria
%
%
% 

%% Inputs
%    C is cell to restor the ellipsoid
%    e is the error tolerance
%    x0 is the initial point or the projected point

%% Outputs
%
% 
%

t1 = cputime;
ck_I = 1;

xk_1 = zeros(n,m+1); % the matirx to store all the vertorx_i^{k-1} to each C_i as Column vector
xk =zeros(n,m+1);
yk_1 = zeros(n,m+1);
yk = zeros(n,m+1);

x_linshi = zeros(n,m+1);
y_linshi = zeros(n,m+1);

%% The first iteration
k=1;
for i=1:m
    xk_1(:,i) = 0;
    yk_1(:,i) = 0;
end

xk_1(:,m+1) = x0;
yk_1(:,m+1) = 0;


%% Later iterations
while ck_I>e
    ck_I=0;
    xk(:,1) = xk_1(:,m+1);
    yk(:,1) = 0;
    for i=2:m+1
        y_linshi = yk;
        x_linshi = xk;
        xk(:,i) = Projection(C{1,i-1},C{2,i-1},C{3,i-1},xk(:,i-1) - yk_1(:,i));
        yk(:,i) = xk(:,i) - (xk(:,i-1) - yk_1(:,i));
        xk_1 = x_linshi;
        yk_1 = y_linshi;
        ck_I =ck_I + norm(yk(:,i)-yk_1(:,i))^2;
    end
    
    k=k+1;
    
    if k>Max_iteration
        break;
    end
    
end
index_iteration = k;
x = xk(:,m+1);
error =  ck_I;

t2 = cputime;

time  = t2-t1;
end

