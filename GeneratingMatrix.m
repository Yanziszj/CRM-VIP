function [ M ] = GeneratingMatrix(n,type_problem)
%% GENERATINGMATRIX 
%   This file to generate three different type of matrix which is to
%   constrast the paramonotoe operator and monotone operator
%

%% Inputs


%% Outputs
%  M


%% Initialinazation

M = zeros(n,n);

switch type_problem
    %% Paramonotone and is a subgradient of a convex function
    % F(x)  = Mx +alpha x^3
    case 1
        eigenvalues = rand(n,1)* 5;         % 
        eigenvalues(n) = 0 ;                     % the last eigenvalue is 0;
        Q = orth(randn(n));                      % Random Orthogonal Matrix
        M = Q * diag(eigenvalues) * Q';  % Ensure symmetry

%     case 2
%         M = triu(rand(n));              % Random upper triangular matrix
%         for i = 1:n
%             % Ensure diagonal dominance
%             M(i,i) = sum(abs(M(i,:))) + rand(); 
%         end
%         M = M + n * eye(n); 
%    

%%
    case 2
        n1 = floor (n/2);
        n2 = n-n1;
        
        M1 = triu(rand(n1));              % Random upper triangular matrix
        for i = 1:n1
            % Ensure diagonal dominance
            M1(i,i) = sum(abs(M1(i,:))) + rand(); 
        end
        M1 = M1 + n1 * eye(n1); 
        for i = 2:n1
            for j = 1:i-1
                M1(i,j) = -M1(j,i);  % Бо A(i,j) = -A(j,i)
            end
        end 
        
        % to construct M2 which is a sysmetric, postive semidefinite metrix
        % whith n2 dimension
        eigenvalues = rand(n2,1)*5;         % 
        eigenvalues(n2) = 0 ;                     % the last eigenvalue is 0;
        Q = orth(randn(n2));                      % Random Orthogonal Matrix
        M2 = Q * diag(eigenvalues) * Q';  % Ensure symmetry
        
        Z1 = zeros(n1,n2);
        Z2 = zeros(n2,n1);
   
        M= [M1 Z1;Z2 M2]; 
        
    case 3
        n1 = floor (n/(10/8));
        n2 = n-n1;
        
        M1 = triu((rand(n1)-0.5)*10);              % Random skew-symmetric matrix
        for i = 1:n1
            % Ensure diagonal is 0
%             M1(i,i) = sum(abs(M1(i,:))) + rand(); 
            M1(i,i) = 0;
        end
%         M1 = M1 + n1 * eye(n1); 
        for i = 2:n1
            for j = 1:i-1
                M1(i,j) = -M1(j,i);  % Бо A(i,j) = -A(j,i)
            end
        end 
        
        % to construct M2 which is a sysmetric, postive semidefinite metrix
        % whith n2 dimension
        eigenvalues = rand(n2,1)* 5;         % 
        eigenvalues(n2) = 0 ;                     % the last eigenvalue is 0;
        Q = orth(randn(n2));                      % Random Orthogonal Matrix
        M2 = Q * diag(eigenvalues) * Q';  % Ensure symmetry
        
        Z1 = zeros(n1,n2);
        Z2 = zeros(n2,n1); 
        
        M= [M1 Z1;Z2 M2]; 
end
        
end

