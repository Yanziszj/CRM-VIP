function [ C, Slater_point ] = GeneratingPolyhedral(n,m)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
%  each polyhydral is in the form of
%      f_i(x) = (a^i)^T x - b_i \leq 0
%      the Salter point is original point to garantee this, 
%             b_i \geq 0 for all i=1,2,...

C = cell(2,m);
lambda = 10;
mu = 5;

for i=1:m  % to generated m polyhedral
    C{1,i} = (rand(n,1)-0.5) * lambda;
    C{2,i} = rand(1,1)*mu;

end
Slater_point = zeros(n,1);


