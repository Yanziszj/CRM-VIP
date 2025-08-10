% function [C,c2,p1,p,h] = GeneratingExamples(n,m)
function [C,Salter_ponit] = GeneratingEllipsoid2(n,m)
%% GENERATINGEXAMPLES is to generate the m sets
%  Input: n is the dimension, m is the number of convex closed sets
%  Output: C is a cell containing all the example
%                        C(1,i) is the positive matrix A, C(2,i) is the
%                        vector b and C(3,i) is the constant.
% we need m is odd.

%% Generate the cell to store the information about the m ellipsoids
C = cell(3,m);
lambda=0.7;
m2=m/2;

%% To generate the ellipsoid 
% first to generate the center of the ellipsoid
for i=1:m2
    center = rand*5+5;
    d = center*1.2;
    u= rand(n,1)*d*lambda;
    % to generate the matrix a
    A=eye(n);
    c=zeros(n,1);
    
    
    for j=1:n
        if j==i;
            c(j)=center;
            A(j,j) = d;
        else
            A(j,j)= u(j);
        end
 
    end
%     Q = rand(n,n);   
%     for k=1:n
%         Q(k,1) = d(k);
%     end
%     Q = GramSchmidt(Q);
%     A = Q*A*Q';
    A = inv(A'*A);
 
    C{1,2*i-1} = A;
    C{2,2*i-1} = -A*c;
    C{3,2*i-1} = -c'*A*c+1;

    C{1,2*i} = A;
    C{2,2*i} = -A*(-c);
    C{3,2*i} = -(-c)'*A*(-c)+1;

end

Salter_ponit = zeros(n,1);
end


