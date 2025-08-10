function [] = Plot2DEllipsoid(C,m)
%% PLOT2DELLIPSOID.m is to plot the all the ellipsoid in xy-plane and the history of xk
%
%   n: dimension
%   m: number of the sets, 
%   X: is all the points for Cyclic cCRM
%   Y: is all the points for the Most Violatedconstraints in function value
%   

figure(1); % to show the trejectory of the X

for i=1:m
    A = C{1,i};
    b = C{2,i};
    c = C{3,i};
    f=@(x,y)(A(1,1)*x^2 + (A(1,2)+A(2,1))*x*y+A(2,2)*y^2 +2*b(1)*x +2*b(2)*y - c);
    ezplot(f,[-5 5 -5 5]);
    hold on;
end
% 
% xlim([-15 15]);
% ylim([-15 15]);

% plot(X(1,:),X(2,:),'r--');
% hold on;
% grid on;

% figure(2);
% for i=1:m
%     A = C{1,i};
%     b = C{2,i};
%     c = C{3,i};
%     f=@(x,y)(A(1,1)*x^2 + (A(1,2)+A(2,1))*x*y+A(2,2)*y^2 +2*b(1)*x +2*b(2)*y - c);
%     ezplot(f);
%     hold on;
% end
% 
% plot(Y(1,:),Y(2,:),'r--');
% hold on;
 grid on;

end

