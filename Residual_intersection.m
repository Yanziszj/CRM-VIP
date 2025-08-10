function r = Residual_intersection(x,y_bf,lambda_bf,B,b_Bar,x0,alpha,n,m)
%% RESIDUAL is here to check the optimal condition of projection onto intersections of ellipsoids
%   
%  x \in R^{n \times 1}
%  y \in R^{nm \times 1}
%  lambda_bf \in R^{nm \times 1}
%  B \in R^{nm \times m}
%  b_Bar \in R^{nm \times 1}
%  x0 \in R^{n \times 1}
%  alpha \in R^{m \times 1}:   store the constant term in each ellipsoid
%  n: dimension of the problem
%  m: number of the ellipsoid
%
%  The difference of Residual_intersection and Residual is, to slove the
%  projection onto intersection of ellipsoids, in the equalent form (6.1),
%  we have more constraint w.r.t. varaible y_bf i.e. we have one
%  constraints for each y_i. Consequently, in R(x,y,lambada) in Algorithm
%  8, we have m terms w.r.t. each y_i
%

y = reshape (y_bf,n,m);

lambda = reshape(lambda_bf,n,m);

b_Bar = reshape(b_Bar,n,m);
R_x = x - x0 -B'*lambda_bf;

R_y = zeros(n,m);  % residual of y_i for i =1,..,m
P_y = zeros(n,m);  % projection of y_i -lambda_i onto B_i

for i=1:m
    if norm(y(:,i) -lambda(:,i))^2 <=alpha(i)+norm(b_Bar(:,i))^2
        P_y(:,i) = y(:,i) -lambda(:,i);
    else
        P_y(:,i) = (sqrt(alpha(i)+norm(b_Bar(:,i))^2))/(norm(y(:,i) -lambda(:,i)))*(y(:,i) -lambda(:,i));
    end
    % P_y = (sqrt(alpha+norm(b_Bar)^2))/(norm(y -lambda))*(y -lambda);
    R_y(:,i) =y(:,i) - P_y(:,i);
end

R_y_bf=reshape(R_y,n*m,1);

b_Bar = reshape(b_Bar,n*m,1);

R_lambda_bf = B*x -y_bf -b_Bar;

R=[R_x;R_y_bf;R_lambda_bf];

r=norm(R);
end

