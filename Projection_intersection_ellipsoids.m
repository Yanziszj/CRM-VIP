function [ x, time, index_iteration,error ] = Projection_intersection_ellipsoids(C, n, m,e,Max_iteration, x0)
%% PROJECTION_INTER_ELLIPSOID 
%  This file is to Implement an algorithm 8 for problem (6.2) in paper
%    Jia, Zehui, Xingju Cai, and Deren Han. "Comparison of several fast algorithms for projection onto an ellipsoid." 
%    Journal of Computational and Applied Mathematics 319 (2017): 320-337.
% Generally speaking it is to find the projection of some point onto
% intesections of m ellipsoid by ultilizing ADMM

t1 =cputime;
% Max_iteration = 1000;

R_history=[]; % to store each residual in each iteration
X =[];

B_block = cell(1,m);

b_bar_block = cell(1,m);
B=[];
b_bar=[];
r=zeros(m,1);
Y=zeros(n,m);
%% Intitialzation
% y = x0;
y =zeros(n,m);

y_bf = zeros(n*m,1);

% lambda = zeros(n,m);  
lambda= rand(n,m)*0.2;

lambda_bf = reshape (lambda,m*n,1); % \in R^{mn,1}

y=rand(n,m);


w=zeros(n,m);



%y = 1/(sqrt(upsilon))*x0;

I =eye(n);

for i=1:m % to define B = [B1;B2;...;Bm] \in R^{nm \times n}
    B_block{1,i} = chol(C{1,i});
    B=[B;B_block{1,i}];
%     b_bar_block{1,i} = -chol(inv(C{1,i}))*C{2,i};
    b_bar_block{1,i} = -inv(B_block{1,i})*C{2,i};
    b_bar=[b_bar;b_bar_block{1,i}];
end

% vartheta = norm(inv(B))/norm(B) ;
% vartheta =norm(inv(B*B'))/norm(B*B');
% vartheta =norm(inv(B'*B))/norm(B'*B);
% vartheta_m  = zeros(m,1);
% for i=1:m
%     vartheta_m(i) = norm(inv(B_block{1,i}))/norm(B_block{1,i});
% end
% vartheta = max(vartheta_m);
% vartheta = zeros(m,1);
% for i=1:m
%     vartheta(i) = norm(inv(B_block{1,i}))/norm(B_block{1,i});
% end

% vartheta = norm(inv(B*B')*B,2)/norm(B,2);
% vartheta = norm(inv(B),2)/norm(B,2);

%b_Bar = - inv(B)'*b;
%A_Bar = inv(I +vartheta*B'*B);
%r=sqrt(alpha+norm(b_Bar)^2);

for i=1:m
    r(i) = sqrt(C{3,i}+norm(b_bar_block{1,i})^2);
end

A_Bar =inv(I + vartheta*B'*B );   %\in R^{n \times n}

u = x0 + B'*lambda_bf +vartheta*B'*(y_bf+b_bar);  %\in R^{n \times n}

x = A_Bar*u; % nm


for i=1:m % to calculate w = [w1 w2 ... wm] \in R^{n \times m}
    
    lambda = reshape(lambda_bf,n,m);
    w(:,i) = B_block{1,i}*x - (1/vartheta)*lambda(:,i)-b_bar_block{1,i};
        if norm(w(:,i))<=r(i)
            y(:,i)=w(:,i);
        else
            y(:,i) = r(i)/(norm(w(:,i)))*w(:,i);
        end
end




y_bf = reshape(y,m*n,1); % mn

% for i=1:m % to calculate lambada = [lambda_1 lambda_2 ... lambda_m] \in R^{n\times m}
%     ambda(:,i) = lambda(:,i) - vartheta*(B*x - y - b_Bar);
% end

lambda_bf = lambda_bf - vartheta*(B*x - y_bf - b_bar); % to calculate lambada = [lambda_1 lambda_2 ... lambda_m] \in R^{n\times m}
     
%  x_bf=[]; 
%  x0_bf=[];
%  for i=1:m
%      x_bf=[x_bf;x];
%      x0_bf = [x0_bf;x0];
%  end

%% Following while is the main loop of Algorithm 8 

alpha=[];
for i=1:m
    alpha = [alpha;C{3,i}];
end

index_iteration=1;

while Residual_intersection(x, y_bf, lambda_bf, B, b_bar, x0,alpha,n,m)>e %
    
    
    u = x0 + B'*lambda_bf +vartheta*B'*(y_bf+b_bar); % Step 3 in Algorithm 8 
    
    x = A_Bar*u;
 
    for i = 1:m % Step 4 in Algorithm 8
        lambda = reshape(lambda_bf,n,m);
        w(:,i) = B_block{1,i}*x - (1/vartheta)*lambda(:,i)-b_bar_block{1,i};
        
        if norm(w(:,i))<=r(i)
            y(:,i)=w(:,i);
        else
            y(:,i) = r(i)/(norm(w(:,i)))*w(:,i);
        end
    end
      
    y_bf = reshape(y ,m*n,1);
%     y = r/(norm(w))*w;
    lambda_bf = lambda_bf - vartheta*(B*x - y_bf - b_bar);
    
    X=[X x];
    R_history=[R_history;Residual_intersection(x, y_bf, lambda_bf, B, b_bar, x0,alpha,n,m)];
    
    index_iteration=index_iteration+1;
    if index_iteration>Max_iteration
        break;
    end
    
end
error = Residual_intersection(x, y_bf, lambda_bf, B, b_bar, x0,alpha,n,m);
t2 = cputime;
time = t2-t1;

end

