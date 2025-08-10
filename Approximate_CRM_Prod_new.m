function [X,xk,index_stopping,Error,time,index_iteration ] = Approximate_CRM_Prod_new( C,e,Max_iteration,xk )
%CRM_PROD 此处显示有关此函数的摘要
%   

[n,m]=size(C);
[n,n] = size(C{1,1});
x_temp = zeros(n,1);
% Max_iteration = 3000;
Error = zeros(m,Max_iteration)+1000; % store the whole error
index_stopping = 'xk and xk+1 are too close';
index_iteration = 1;
X=[];
t1 = cputime;
x_temp = 0;
% while (norm(xk-x_temp)>e)
while (1)
    x_temp = xk;
    
    % choose the type of epsilon

    wk=zeros(n,1);
    s=zeros(n,1);
    for i=1:m  %the cyclic loop
        vik = xk - Approximate_Projection(C{1,i},C{2,i},C{3,i},xk);
        wk =wk +vik;
        
        s(i) = norm(vik)^2;
    end
    wk =wk/m;
    
    % calculate the alpha_k
    alpha_k = (sum(s))/(m*norm(wk)^2);
    
    xk= xk -alpha_k*wk;
    
    for i=1:m
        tem=xk'*C{1,i}*xk+2*C{2,i}'*xk -C{3,i};
%         tem=xk'*C{1,i}*xk+2*C{2,i}'*xk -C{3,i};
        Error(i,index_iteration) = max([0,tem]);
    end                
    
%     if max(Error(:,index_iteration))*10^(index_iteration)<e   % to test the finite convergence
%     if max(Error(:,index_iteration))<e
    if max(Error(:,index_iteration))==0
        index_stopping = 'successful';
        break;
    end     
    
    % the max-iteration check
    if index_iteration>=Max_iteration
        index_stopping = 'the iteration exceed the max';  % to indicate that program stops since the interation exceed the max
        break;    % stop the cyclic loop
    end

    index_iteration = index_iteration+1;
%             X = [X xk];           
end
t2 = cputime;
time = t2-t1;
error=Error(:,index_iteration);
end



