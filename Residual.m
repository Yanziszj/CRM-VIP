function r = Residual(x,y,lambda,B,b_Bar,x0,alpha)
%RESIDUAL 此处显示有关此函数的摘要
%   此处显示详细说明
R_x = x - x0 -B'*lambda;

if norm(y -lambda)^2 <=alpha+norm(b_Bar)^2
    P_y = y -lambda;
else
    P_y = (sqrt(alpha+norm(b_Bar)^2))/(norm(y -lambda))*(y -lambda);
end
% P_y = (sqrt(alpha+norm(b_Bar)^2))/(norm(y -lambda))*(y -lambda);
R_y =y - P_y;
R_lambda = B*x -y -b_Bar;

 R=[R_x;R_y;R_lambda];
r=norm(R);
end

