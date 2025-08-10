function [ index_max ] = Find_max_function_value( C,n,m,xk)
%% FIND_MAX_FUNCTION_VALUE 
%   To find the index of max function value at xk 
%   functions are in related to ellipsoid stored in C

%% Inputs
%

%% Outputs
%

Rank_FunctionValues=zeros(1,m); % restore the function value of g_i at xk
for j=1:m % calculate the largest function value
    
    Rank_FunctionValues(j) = xk'*C{1,j}*xk +2*C{2,j}'*xk - C{3,j};
end
[value_max,index_max] =max(Rank_FunctionValues);
% max_value=max(Rank_FunctionValues);
% Max_index= find(max_value==Rank_FunctionValues); 
end

